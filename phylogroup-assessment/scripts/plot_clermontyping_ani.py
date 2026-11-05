#!/usr/bin/env python3
"""Plot a 300-DPI ANI clustermap with qualified two-of-three phylogroup voting.

Run: bash run_ani_phylgroup_prediction.sh (fixed seed 13), or
python3 plot_clermontyping_ani.py --seed 13 [--output figure.png]
Requires numpy, pandas, scipy, matplotlib, seaborn, and biopython.

Inputs are read from ../input and default outputs go to ../results, relative
to this script, not the working directory. Only genomes in
../input/EcoZUR_genomovar_metadadta.tsv are included.
ClermonTyping_results supplies phylogroups; ClermonTyping_metadata.tsv is used
to check those assignments against the original assembly-named records.
Also reads the corresponding ../input/ClermonTyping_output/*_mash_screen.tab and *.xml
files and the V2_gamma Newick tree specified by TREE_PATH below. Other genomes
in those inputs are excluded. Original input files/assignments are never edited.

For each unordered genome pair (including self matches), retain the ANI from
the row with the largest FIFTH column (total query fragments), as in the
requested example. Equal maxima are sampled uniformly at random; --seed makes
this reproducible. If omitted, a random seed is generated and recorded in the
log and evidence TSV. The selected ANI fills both symmetric matrix positions.
Missing comparisons or annotations are errors, never imputed values.

Both axes use the same complete-linkage clustering of distance = 100 - ANI,
with zero self-distances. Original self-match ANI values remain in the heatmap,
including small deviations from 100 reported by fastANI.
Genomovar colors come directly from the metadata. Phylogroups use a colorblind
palette, with Unknown in gray. Labels appear once per contiguous annotation
block, using whichever of black or white has the higher contrast.

Prediction rules (heuristics, not calibrated probabilities or definitive calls):
* References must have an original assigned PCR call agreeing with a unique
  best Mash phylogroup at identity >=0.95, and no conflicting assigned PCR calls
  in their genomovar. Mash screen identity is NOT fastANI or a probability.
* Exclude the focal genome from every comparison. Prefer trusted same-genomovar
  peers at ANI >=99.5, provided the genomovar has no conflicting original calls.
* Otherwise exclude the ENTIRE focal genomovar. ANI votes only if its best match
  exceeds 98.5%, at least ONE external genomovar supports it above 98.5%, and
  its advantage over the best competing phylogroup is >=0.1 percentage point.
* Tree independently votes for its closest reference group only if its distance
  is <=0.8 times that to the closest competing group. Its candidate is NOT
  restricted to the ANI candidate. Missing competitors/zero-distance ties abstain.
* Mash votes for its unique best reference group only at identity >=0.95.
  Each method has one vote; a failed quality check is an abstention (Unknown).
* Any TWO qualified methods agreeing determine consensus, irrespective of the
  third method. No individual method vetoes consensus; predicted labels are
  never recycled as references.
* Without consensus, a singleton with an assigned original call retains that
  call, flagged singleton_original_fallback. Otherwise the result is Unknown.
  The third bar displays final assignments, INCLUDING singleton fallbacks.

Writes a per-genome evidence TSV (default: ClermonTyping_phylogroup_predictions.tsv).
It includes original calls, candidate/prediction/status/reasons, neighbors,
fastANI, tree distances, Mash hits, and PCR/XML audit details.
manuscript_phylogroup uses the final assignment when available and otherwise
retains the original ClermonTyping call. The PNG still uses predicted_phylogroup,
including singleton fallbacks and Unknown for unresolved cases. Fallbacks do not
validate the original call or feed back into prediction. Votes, abstention
reasons, consensus, and assignment_basis distinguish support from fallback.
XML hit presence
alone is not a valid PCR product; the recorded PCR result is not recalculated.
Missing optional Mash/XML files are warned about and recorded, never silently
treated as negative hits. Missing Mash evidence makes a genome ineligible as a
reference. Malformed files or missing required tree tips are errors.

Method/file-format sources:
https://github.com/A-BN/ClermonTyping
https://github.com/A-BN/ClermonTyping/blob/master/bin/clermont.py
https://github.com/A-BN/ClermonTyping/blob/master/bin/add_mash_minimal.R
https://mash.readthedocs.io/en/latest/tutorials.html
"""

import argparse
import ast
import csv
import logging
import random
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET

from Bio import Phylo
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_hex, to_rgb
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
import seaborn as sns


BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent / "input"
RESULTS_DIR = BASE_DIR.parent / "results"
DPI = 300
LOGGER = logging.getLogger(__name__)
TREE_PATH = INPUT_DIR / "EcoZUR_fasttree.tre"
WITHIN_ANI = 99.5
EXTERNAL_ANI = 98.5
MIN_ANI_MARGIN = 0.1
MAX_TREE_RATIO = 0.8
MIN_EXTERNAL_GENOMOVARS = 1
MASH_MIN_IDENTITY = 0.95


def genome_id(name: str) -> str:
    return Path(name).name.split("_SPAdes_", 1)[0].replace("_", "-")


def phylogroup_label(value: str) -> str:
    value = value.strip()
    if value.lower() in {"unknown", "unassigned", "undetermined", "na", "n/a"}:
        return "Unknown"
    return value


def load_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = ("genome", "genomovar_complete", "genomovar_complete_colors")
    missing = set(required) - set(metadata.columns)
    if missing:
        raise ValueError(f"{path}: missing columns: {sorted(missing)}")
    if len(metadata) < 2:
        raise ValueError(f"{path}: at least two genomes are required")
    for column in required:
        if metadata[column].str.strip().eq("").any():
            raise ValueError(f"{path}: empty values in {column}")
    duplicates = metadata.loc[metadata["genome"].duplicated(), "genome"].tolist()
    if duplicates:
        raise ValueError(f"{path}: duplicate genome IDs: {duplicates}")
    for genome, color in zip(metadata["genome"], metadata[required[2]]):
        if not is_color_like(color) or color.lower() == "none":
            raise ValueError(f"{path}: invalid color {color!r} for {genome}")
    canonical_colors = metadata[required[2]].map(to_hex)
    color_counts = canonical_colors.groupby(metadata["genomovar_complete"]).nunique()
    if color_counts.gt(1).any():
        raise ValueError(f"{path}: a genomovar has inconsistent colors")
    return metadata.set_index("genome")


def load_phylogroups(
    path: Path, genomes: pd.Index, *, raw_metadata: bool = False
) -> pd.Series:
    assignments = {}
    allowed = set(genomes)
    with path.open() as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields:
                continue
            genome = fields[0]
            if raw_metadata:
                genome = genome_id(genome)
            if genome not in allowed:
                continue
            column = 4 if raw_metadata else 1
            if len(fields) <= column or not fields[column].strip():
                raise ValueError(f"{path}:{line_number}: missing phylogroup")
            if genome in assignments:
                raise ValueError(f"{path}:{line_number}: duplicate genome {genome}")
            assignments[genome] = phylogroup_label(fields[column])
    missing = allowed - assignments.keys()
    if missing:
        raise ValueError(f"{path}: missing phylogroups for {sorted(missing)}")
    return pd.Series(assignments, name="ClermonTyping phylogroup").reindex(genomes)


def load_ani_matrix(path: Path, genomes: pd.Index, rng: random.Random) -> pd.DataFrame:
    allowed = set(genomes)
    selected = {}
    tie_counts = {}
    retained_rows = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: missing genome IDs")
            query, reference = fields[:2]
            if query not in allowed or reference not in allowed:
                continue
            if len(fields) != 5:
                raise ValueError(f"{path}:{line_number}: expected five columns")
            try:
                ani = float(fields[2])
                matched, total = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid ANI/fragments") from error
            if not np.isfinite(ani) or not 0 <= ani <= 100:
                raise ValueError(f"{path}:{line_number}: ANI must be within [0, 100]")
            if not 0 <= matched <= total or total <= 0:
                raise ValueError(f"{path}:{line_number}: invalid fragment counts")
            retained_rows += 1
            pair = tuple(sorted((query, reference)))
            if pair not in selected or total > selected[pair][1]:
                selected[pair] = (ani, total)
                tie_counts[pair] = 1
            elif total == selected[pair][1]:
                # Reservoir sampling also handles more than two tied rows fairly.
                tie_counts[pair] += 1
                if rng.randrange(tie_counts[pair]) == 0:
                    selected[pair] = (ani, total)

    matrix = pd.DataFrame(np.nan, index=genomes, columns=genomes)
    for (query, reference), (ani, _) in selected.items():
        matrix.at[query, reference] = ani
        matrix.at[reference, query] = ani
    missing = np.argwhere(np.triu(matrix.isna().to_numpy()))
    if len(missing):
        examples = [(genomes[i], genomes[j]) for i, j in missing[:10]]
        raise ValueError(
            f"{path}: {len(missing)} missing unordered pairs (including self matches); "
            f"first examples: {examples}"
        )
    LOGGER.info(
        "Retained %d ANI rows -> %d unordered pairs, including %d self matches; "
        "%d pairs had tied maximum fifth-column values",
        retained_rows, len(selected), len(genomes),
        sum(count > 1 for count in tie_counts.values()),
    )
    return matrix


def load_mash_evidence(path: Path) -> dict:
    result = {
        "mash_phylogroup": "Unknown", "mash_identity": np.nan,
        "mash_reference": "", "mash_runner_up_identity": np.nan,
        "mash_status": "missing",
    }
    if not path.exists():
        LOGGER.warning("Missing Mash evidence: %s", path)
        return result
    best_by_group = {}
    with path.open() as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields:
                continue
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_number}: expected Mash screen columns")
            try:
                identity = float(fields[0])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid Mash identity") from error
            if not np.isfinite(identity) or not 0 <= identity <= 1:
                raise ValueError(f"{path}:{line_number}: invalid Mash identity")
            parts = Path(fields[4]).name.rsplit("_", 2)
            if len(parts) != 3 or not parts[1]:
                raise ValueError(f"{path}:{line_number}: reference lacks phylogroup tag")
            group = parts[1]
            hit = (identity, fields[4])
            if group not in best_by_group or hit > best_by_group[group]:
                best_by_group[group] = hit
    if not best_by_group:
        LOGGER.warning("No Mash hits: %s", path)
        result["mash_status"] = "no_hits"
        return result
    ranked = sorted(best_by_group, key=lambda group: (-best_by_group[group][0], group))
    group = ranked[0]
    identity, reference = best_by_group[group]
    runner_up = best_by_group[ranked[1]][0] if len(ranked) > 1 else np.nan
    tied = len(ranked) > 1 and identity == runner_up
    result.update(
        mash_phylogroup="Unknown" if tied else group,
        mash_identity=identity, mash_reference=reference,
        mash_runner_up_identity=runner_up,
        mash_status="tied_groups" if tied else (
            "supported" if identity >= MASH_MIN_IDENTITY else "low_identity"
        ),
    )
    return result


def load_typing_evidence(path: Path, output_dir: Path, genomes: pd.Index) -> pd.DataFrame:
    records = {}
    with path.open() as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields or genome_id(fields[0]) not in genomes:
                continue
            genome = genome_id(fields[0])
            if len(fields) != 6 or genome in records:
                raise ValueError(f"{path}:{line_number}: invalid or duplicate typing record")
            parsed = []
            for value in fields[1:4]:
                try:
                    markers = ast.literal_eval(value)
                except (ValueError, SyntaxError) as error:
                    raise ValueError(f"{path}:{line_number}: invalid marker list") from error
                if not isinstance(markers, list) or not all(isinstance(x, str) for x in markers):
                    raise ValueError(f"{path}:{line_number}: expected a list of marker strings")
                parsed.append(markers)
            if len(parsed[1]) != 4 or not set(parsed[1]) <= {"+", "-"}:
                raise ValueError(f"{path}:{line_number}: invalid quadruplex pattern")
            if Path(fields[5]).name != fields[5] or genome_id(
                fields[5].removesuffix("_mash_screen.tab")
            ) != genome:
                raise ValueError(f"{path}:{line_number}: Mash filename does not match genome")
            record = {
                "pcr_products": ";".join(parsed[0]),
                "pcr_arpA_chuA_yjaA_TspE4C2": "".join(parsed[1]),
                "pcr_specific_markers": ";".join(parsed[2]),
                **load_mash_evidence(output_dir / fields[5]),
                "xml_status": "missing", "xml_primer_queries": np.nan,
                "xml_queries_without_hits": "",
            }
            xml_path = output_dir / (Path(fields[0]).name + ".xml")
            if not xml_path.exists():
                LOGGER.warning("Missing primer BLAST XML audit evidence: %s", xml_path)
            else:
                root = ET.parse(xml_path).getroot()
                if root.findtext("BlastOutput_program") != "blastn":
                    raise ValueError(f"{xml_path}: expected BLASTN XML")
                database = root.findtext("BlastOutput_db")
                if not database or genome_id(database) != genome:
                    raise ValueError(f"{xml_path}: BLAST database does not match genome")
                queries = root.findall("./BlastOutput_iterations/Iteration")
                if not queries:
                    raise ValueError(f"{xml_path}: no primer query records")
                no_hits = []
                for query in queries:
                    name = query.findtext("Iteration_query-def")
                    if not name:
                        raise ValueError(f"{xml_path}: missing primer query name")
                    if not query.findall(".//Hsp"):
                        no_hits.append(name)
                record.update(
                    xml_status="parsed", xml_primer_queries=len(queries),
                    xml_queries_without_hits=";".join(no_hits),
                )
            records[genome] = record
    missing = set(genomes) - records.keys()
    if missing:
        raise ValueError(f"{path}: missing typing evidence for {sorted(missing)}")
    return pd.DataFrame.from_dict(records, orient="index").reindex(genomes)


def load_tree_distances(path: Path, genomes: pd.Index) -> pd.DataFrame:
    tree = Phylo.read(path, "newick")
    tips = {}
    for tip in tree.get_terminals():
        if not tip.name:
            raise ValueError(f"{path}: unnamed tree tip")
        name = genome_id(tip.name)
        if name in tips:
            raise ValueError(f"{path}: duplicate normalized tree tip {name}")
        tips[name] = tip
    missing = set(genomes) - tips.keys()
    if missing:
        raise ValueError(f"{path}: missing genomes: {sorted(missing)}")
    for clade in tree.find_clades():
        length = clade.branch_length
        if length is None and clade is tree.root:
            continue
        if length is None or not np.isfinite(length) or length < 0:
            raise ValueError(f"{path}: non-root branches require finite nonnegative lengths")
    distances = np.zeros((len(genomes), len(genomes)))
    for i, query in enumerate(genomes):
        for j in range(i):
            distances[i, j] = distances[j, i] = tree.distance(tips[query], tips[genomes[j]])
    return pd.DataFrame(distances, index=genomes, columns=genomes)


def tree_support(
    query: str, references: list[str], candidate: str,
    distances: pd.DataFrame, labels: pd.Series,
) -> dict:
    ranked = sorted(references, key=lambda ref: (distances.at[query, ref], ref))
    nearest = ranked[0] if ranked else None
    same = next((r for r in ranked if labels[r] == candidate), None)
    other = next((r for r in ranked if labels[r] != candidate), None)
    same_distance = distances.at[query, same] if same else np.nan
    other_distance = distances.at[query, other] if other else np.nan
    # A zero-distance tie is deliberately not treated as strong separation.
    ratio = same_distance / other_distance if other_distance > 0 else np.nan
    return {
        "tree_group": labels[nearest] if nearest else "Unknown",
        "tree_reference": nearest or "",
        "tree_distance": distances.at[query, nearest] if nearest else np.nan,
        "tree_candidate_reference": same or "", "tree_candidate_distance": same_distance,
        "tree_competitor_group": labels[other] if other else "Unknown",
        "tree_competitor_reference": other or "", "tree_competitor_distance": other_distance,
        "tree_ratio": ratio,
    }


def summarize_references(
    query: str, references: list[str], ani: pd.DataFrame, distances: pd.DataFrame,
    labels: pd.Series, genomovars: pd.Series,
) -> dict:
    evidence = {
        "ani_group": "Unknown", "ani_reference": "", "ani": np.nan,
        "ani_competitor_group": "Unknown", "ani_competitor": np.nan,
        "ani_margin": np.nan,
        "support_genomovars": 0, "support_genomes": 0, "supporting_neighbors": "",
        **tree_support(query, [], "Unknown", distances, labels),
    }
    if not references:
        return evidence
    ranked_ani = sorted(references, key=lambda ref: (-ani.at[query, ref], ref))
    best_ani = ranked_ani[0]
    group = labels[best_ani]
    ani_other = next((r for r in ranked_ani if labels[r] != group), None)
    supporters = [r for r in ranked_ani if labels[r] == group and ani.at[query, r] > EXTERNAL_ANI]
    evidence.update(
        ani_group=group, ani_reference=best_ani, ani=ani.at[query, best_ani],
        ani_competitor_group=labels[ani_other] if ani_other else "Unknown",
        ani_competitor=ani.at[query, ani_other] if ani_other else np.nan,
        ani_margin=ani.at[query, best_ani] - ani.at[query, ani_other] if ani_other else np.nan,
        **tree_support(query, references, group, distances, labels),
        support_genomovars=genomovars.loc[supporters].nunique(),
        support_genomes=len(supporters),
        supporting_neighbors=";".join(
            f"{r}|{genomovars[r]}|ANI={ani.at[query, r]:.4f}|tree={distances.at[query, r]:.6g}"
            for r in supporters
        ),
    )
    return evidence


def predict_phylogroups(
    ani: pd.DataFrame, metadata: pd.DataFrame, labels: pd.Series,
    distances: pd.DataFrame, typing: pd.DataFrame,
) -> pd.DataFrame:
    genomes = metadata.index
    for name, frame in (("ANI", ani), ("tree", distances)):
        if not frame.index.equals(genomes) or not frame.columns.equals(genomes):
            raise ValueError(f"{name} matrix must match the metadata genome order")
    if not labels.index.equals(genomes) or not typing.index.equals(genomes):
        raise ValueError("Typing inputs must match the metadata genome order")
    genomovars = metadata["genomovar_complete"]
    group_conflicts = {
        group: labels.loc[part.index][labels.loc[part.index] != "Unknown"].nunique() > 1
        for group, part in metadata.groupby("genomovar_complete")
    }
    trusted = [
        genome for genome in genomes
        if labels[genome] != "Unknown" and not group_conflicts[genomovars[genome]]
        and typing.at[genome, "mash_status"] == "supported"
        and typing.at[genome, "mash_phylogroup"] == labels[genome]
    ]
    LOGGER.info("Using %d original PCR/Mash-concordant references; no predicted labels reused", len(trusted))
    records = []
    for genome in genomes:
        group = genomovars[genome]
        peers = [r for r in genomes if r != genome and genomovars[r] == group]
        external = [r for r in trusted if genomovars[r] != group]
        ext = summarize_references(genome, external, ani, distances, labels, genomovars)
        all_refs = [r for r in trusted if r != genome]
        local = [
            r for r in all_refs if genomovars[r] == group and ani.at[genome, r] >= WITHIN_ANI
        ] if not group_conflicts[group] else []
        ani_reasons = []
        if local:
            method = "within_genomovar"
            candidate = labels[local[0]]
            tree_references = all_refs
        else:
            method = "external_genomovars"
            candidate = ext["ani_group"]
            tree_references = external
            if not external:
                ani_reasons.append("no_external_references")
            if not ext["ani"] > EXTERNAL_ANI:
                ani_reasons.append("external_ANI_not_above_98.5")
            if not ext["ani_margin"] >= MIN_ANI_MARGIN - 1e-10:
                ani_reasons.append("insufficient_ANI_separation")
            if ext["support_genomovars"] < MIN_EXTERNAL_GENOMOVARS:
                ani_reasons.append("no_supporting_external_genomovar")
        ani_vote = candidate if not ani_reasons else "Unknown"
        nearest = min(
            tree_references, key=lambda ref: (distances.at[genome, ref], ref), default=None
        )
        tree_candidate = labels[nearest] if nearest is not None else "Unknown"
        tree_evidence = tree_support(genome, tree_references, tree_candidate, distances, labels)
        tree_vote = (
            tree_candidate if tree_evidence["tree_ratio"] <= MAX_TREE_RATIO else "Unknown"
        )
        mash_vote = (
            typing.at[genome, "mash_phylogroup"]
            if typing.at[genome, "mash_status"] == "supported" else "Unknown"
        )
        votes = {"ANI": ani_vote, "Tree": tree_vote, "Mash": mash_vote}
        counts = Counter(vote for vote in votes.values() if vote != "Unknown")
        consensus = next((group for group, count in counts.items() if count >= 2), "Unknown")
        singleton_fallback = (
            consensus == "Unknown" and not peers and labels[genome] != "Unknown"
        )
        predicted = labels[genome] if singleton_fallback else consensus
        if singleton_fallback:
            basis = status = reason = "singleton_original_fallback"
        elif consensus != "Unknown":
            basis = "two_of_three_consensus"
            reason = "two_of_three_agreement"
            status = "predicted_unassigned" if labels[genome] == "Unknown" else (
                "supported_original" if labels[genome] == consensus else "suggested_revision"
            )
        else:
            basis = status = "unresolved"
            reason = "no_two_of_three_agreement"
            if not peers:
                reason += ";singleton_original_unassigned"
        supporting_methods = ";".join(
            name for name, vote in votes.items() if consensus != "Unknown" and vote == consensus
        )
        records.append({
            "genome": genome, "genomovar": group, "genomovar_size": len(peers) + 1,
            "original_phylogroup": labels[genome], "predicted_phylogroup": predicted,
            "manuscript_phylogroup": labels[genome] if predicted == "Unknown" else predicted,
            "candidate_phylogroup": candidate, "status": status, "method": method,
            "reason": reason, "assignment_basis": basis,
            "consensus_phylogroup": consensus, "supporting_methods": supporting_methods,
            "singleton_fallback": singleton_fallback,
            "ani_vote": ani_vote,
            "ani_vote_reason": ";".join(ani_reasons) if ani_reasons else "qualified",
            "tree_vote": tree_vote,
            "tree_vote_reason": "qualified" if tree_vote != "Unknown" else (
                "no_tree_references" if nearest is None else "insufficient_tree_separation"
            ),
            "mash_vote": mash_vote, "mash_vote_reason": typing.at[genome, "mash_status"],
            "genomovar_conflict": group_conflicts[group], "trusted_reference": genome in trusted,
            "peer_original_calls": ";".join(f"{key}:{value}" for key, value in sorted(Counter(labels.loc[peers]).items())),
            "within_supporting_neighbors": ";".join(
                f"{r}|{labels[r]}|ANI={ani.at[genome, r]:.4f}|tree={distances.at[genome, r]:.6g}"
                for r in local
            ),
            **{f"decision_{key}": value for key, value in tree_evidence.items()},
            **{f"external_{key}": value for key, value in ext.items()},
        })
    return pd.DataFrame(records).set_index("genome").join(typing)


def contrasting_text(color: str) -> str:
    rgb = np.asarray(to_rgb(color))
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    luminance = float(linear @ np.array([0.2126, 0.7152, 0.0722]))
    black_contrast = (luminance + 0.05) / 0.05
    white_contrast = 1.05 / (luminance + 0.05)
    return "black" if black_contrast >= white_contrast else "white"


def annotate_blocks(ax, labels: list[str], colors: list[str], row: int) -> None:
    start = 0
    while start < len(labels):
        end = start + 1
        while end < len(labels) and labels[end] == labels[start]:
            end += 1
        ax.text(
            (start + end) / 2, row + 0.5, labels[start],
            ha="center", va="center", rotation=90 if end - start < 4 else 0,
            fontsize=4.5, color=contrasting_text(colors[start]), clip_on=True,
        )
        start = end


def plot_clustermap(
    matrix: pd.DataFrame, metadata: pd.DataFrame, phylogroups: pd.Series, output: Path,
    predicted: pd.Series | None = None,
) -> None:
    distances = 100.0 - matrix.to_numpy()
    non_100_self = np.count_nonzero(np.diag(distances))
    if non_100_self:
        LOGGER.info(
            "Preserving %d non-100 self-match ANI values in the heatmap; "
            "setting clustering self-distances to zero", non_100_self,
        )
    np.fill_diagonal(distances, 0)
    tree = linkage(squareform(distances, checks=True), method="complete", optimal_ordering=True)
    groups = sorted((set(phylogroups) | (set(predicted) if predicted is not None else set())) - {"Unknown"})
    palette = dict(zip(groups, map(to_hex, sns.color_palette("colorblind", len(groups)))))
    palette["Unknown"] = "#999999"
    color_bars = pd.DataFrame({
        "Genomovar": metadata["genomovar_complete_colors"],
        "ClermonTyping phylogroup": phylogroups.map(palette),
    }).reindex(matrix.columns)
    annotation_labels = [metadata["genomovar_complete"], phylogroups]
    if predicted is not None:
        color_bars["Predicted phylogroup"] = predicted.map(palette)
        annotation_labels.append(predicted)

    # Allow approximately 20 pixels per genome at 300 DPI.
    side = max(12, len(matrix) * 0.08)
    with plt.rc_context({"figure.dpi": DPI, "savefig.dpi": DPI}):
        grid = sns.clustermap(
            matrix, row_linkage=tree, col_linkage=tree, col_colors=color_bars,
            cmap="viridis", vmin=float(matrix.to_numpy().min()), vmax=100,
            figsize=(side, side), dendrogram_ratio=0.10, colors_ratio=0.025,
            cbar_pos=(0.025, 0.90, 0.018, 0.065),
            cbar_kws={"label": "ANI (%)"},
            xticklabels=True, yticklabels=True, rasterized=True,
            tree_kws={"linewidths": 0.5},
        )
        order = grid.dendrogram_col.reordered_ind
        ordered_genomes = matrix.columns[order]
        for row, labels in enumerate(annotation_labels):
            annotate_blocks(
                grid.ax_col_colors,
                labels.loc[ordered_genomes].tolist(),
                color_bars.iloc[order, row].tolist(),
                row,
            )
        grid.ax_col_colors.set_xticks([])
        grid.ax_col_colors.tick_params(axis="y", labelsize=8, length=0, pad=3)
        grid.ax_heatmap.tick_params(axis="both", labelsize=4, length=0, pad=2)
        plt.setp(grid.ax_heatmap.get_xticklabels(), rotation=90)
        plt.setp(grid.ax_heatmap.get_yticklabels(), rotation=0)
        grid.ax_heatmap.set_xlabel("Genome", fontsize=10, labelpad=5)
        grid.ax_heatmap.set_ylabel("Genome", fontsize=10, labelpad=5)
        grid.ax_cbar.tick_params(labelsize=6)
        grid.ax_cbar.set_ylabel("ANI (%)", fontsize=8, labelpad=4)
        grid.fig.suptitle(
            f"EcoZUR ANI and ClermonTyping ({len(matrix)} genomes)",
            fontsize=12, y=1.025,
        )
        grid.savefig(output, dpi=DPI, bbox_inches="tight")
        plt.close(grid.fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, help="Random seed for ties (default: generate and record a random seed)")
    parser.add_argument(
        "--output", type=Path, default=RESULTS_DIR / "ClermonTyping_ANI_clustermap.png",
        help="Output PNG (default: ../results, relative to this script)",
    )
    parser.add_argument(
        "--predictions-output", type=Path,
        default=RESULTS_DIR / "ClermonTyping_phylogroup_predictions.tsv",
        help="Per-genome prediction and supporting evidence TSV",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if args.output.suffix.lower() != ".png":
        parser.error("--output must have a .png extension")
    if args.predictions_output.suffix.lower() != ".tsv":
        parser.error("--predictions-output must have a .tsv extension")
    protected_inputs = {
        (INPUT_DIR / "ClermonTyping_metadata.tsv").resolve(),
        (INPUT_DIR / "ClermonTyping_results").resolve(),
        (INPUT_DIR / "EcoZUR_genomovar_metadadta.tsv").resolve(),
        (INPUT_DIR / "EcoZUR.ani").resolve(),
        TREE_PATH.resolve(),
    }
    if args.predictions_output.resolve() in protected_inputs or args.output.resolve() in protected_inputs:
        parser.error("Output paths must not overwrite an input file")
    try:
        metadata = load_metadata(INPUT_DIR / "EcoZUR_genomovar_metadadta.tsv")
        phylogroups = load_phylogroups(INPUT_DIR / "ClermonTyping_results", metadata.index)
        raw_groups = load_phylogroups(
            INPUT_DIR / "ClermonTyping_metadata.tsv", metadata.index, raw_metadata=True
        )
        mismatches = phylogroups.index[phylogroups.ne(raw_groups)].tolist()
        if mismatches:
            raise ValueError(f"ClermonTyping input files disagree for {mismatches}")
        LOGGER.info("Using %d metadata-listed genomes; ClermonTyping inputs agree", len(metadata))
        seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**32)
        LOGGER.info("ANI tie-breaking seed: %d", seed)
        matrix = load_ani_matrix(
            INPUT_DIR / "EcoZUR.ani",
            metadata.index, random.Random(seed),
        )
        typing = load_typing_evidence(
            INPUT_DIR / "ClermonTyping_metadata.tsv", INPUT_DIR / "ClermonTyping_output",
            metadata.index,
        )
        tree_distances = load_tree_distances(TREE_PATH, metadata.index)
        predictions = predict_phylogroups(matrix, metadata, phylogroups, tree_distances, typing)
        predictions["ani_tie_seed"] = seed
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.predictions_output.parent.mkdir(parents=True, exist_ok=True)
        plot_clustermap(
            matrix, metadata, phylogroups, args.output, predictions["predicted_phylogroup"]
        )
        predictions.to_csv(args.predictions_output, sep="\t", na_rep="NA", float_format="%.8g")
    except (OSError, ValueError, ET.ParseError) as error:
        parser.exit(1, f"ERROR: {error}\n")
    LOGGER.info("Saved %s at %d DPI", args.output.resolve(), DPI)
    LOGGER.info("Prediction outcomes: %s", predictions["status"].value_counts().to_dict())
    LOGGER.info("Saved prediction evidence: %s", args.predictions_output.resolve())


if __name__ == "__main__":
    main()
