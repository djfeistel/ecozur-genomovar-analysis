# EcoZUR phylogroup assessment using ANI, ClermonTyping, and a core-gene tree

## Purpose and scope

This analysis evaluates ClermonTyping phylogroup membership for the [247 *E. coli* genomes that passed QA](../../genome-quality-assessment/results/contamination_vs_completeness.png), including both genomes with and without an assigned ClermonTyping phylogroup assignments.
It combines precomputed fastANI comparisons, genomovar membership grouped using 99.5% ANI,
ClermonTyping PCR calls, Mash reference comparisons, and patristic distances from our core-gene phylogeny.

The phylgorup consensus predictions are **hypotheses**, not confirmed
taxonomic corrections or calibrated probabilities. Original assignments are
preserved alongside predictions. Any two qualified methods agreeing determine
consensus. Without consensus, a singleton retains its assigned original call
as an explicitly flagged fallback; other cases remain `Unknown`.

This package reproduces the downstream assessment and figure. It does not rerun
fastANI, BLAST, Mash, ClermonTyping, genome assembly, or phylogenetic inference.
It requires no original assembly files or network access at analysis runtime.

## Package layout

```text
phylogroup-assessment/
  input/
    EcoZUR.ani
    EcoZUR_genomovar_metadadta.tsv
    ClermonTyping_results
    ClermonTyping_metadata.tsv
    EcoZUR_fasttree.tre
    ClermonTyping_output/
      <assembly>.fa_mash_screen.tab     (247 files)
      <assembly>.fa.xml                 (247 files)
    SHA256SUMS
  scripts/
    plot_clermontyping_ani.py
    run_ani_phylgroup_prediction.sh
    requirements.txt
    readme.md
  results/
    ClermonTyping_ANI_clustermap.png
    ClermonTyping_phylogroup_predictions.tsv
```

The Python script resolves its input and default output locations relative to
its own location. The whole package can be moved without editing absolute
paths. Keep the `input`, `scripts`, and `results` layout together.

## Execution

```bash
bash run_ani_phylgroup_prediction.sh
```

Outputs are written to [../results/](../results/). Rerunning the launcher
**overwrites these two generated outputs**, but never modifies the inputs.
Failures return a nonzero exit status; messages are written to standard error.

Equivalent direct invocation, from this directory:

## Input data and provenance

The inputs were copied from the original EcoZUR analysis; the source files were
left unchanged.

From this directory, verify the input copies with:

```bash
(cd ../input && sha256sum --check SHA256SUMS)
```

### Genome membership and colors

[EcoZUR_genomovar_metadadta.tsv](../input/EcoZUR_genomovar_metadadta.tsv)

1. `genome`
2. `genomovar_complete`
3. `genomovar_complete_colors`

This table defines the complete analysis population. Genomes absent from it
are not evaluated or used as references, even if they occur in another input.
Existing genomovar labels are used as supplied; genomovars are not recalculated.

### fastANI input

[EcoZUR.ani](../input/EcoZUR.ani)

Columns represent: query genome, reference genome, ANI percent, bidirectional mapped fragment count, and total fragment count of the query. 
The source file covers 248 genomes, but only the 247 metadata-listed genomes enter this analysis as B37-6 shows <95% ANi with the other genomes

### ClermonTyping inputs

- [ClermonTyping_results](../input/ClermonTyping_results): a headerless two-column
  table of genome IDs and original phylogroups.
- [ClermonTyping_metadata.tsv](../input/ClermonTyping_metadata.tsv): a headerless
  six-column table containing the assembly filename, detected PCR products,
  four-marker pattern, additional specific markers, original phylogroup, and
  associated Mash table filename.
- [ClermonTyping_output/](../input/ClermonTyping_output/): only the 247 Mash
  tables and 247 BLAST XML files used for included genomes were copied.


ClermonTyping identifies plausible PCR products using primer alignments,
orientation, spacing, and mismatch criteria, then uses marker patterns to assign
phylogroups. The quadruplex order is **arpA, chuA, yjaA, TspE4.C2**.

The XML files are BLAST results against primers, not separate phylogroup calls.
This script records query counts and primers with no hits; it does not rerun
the PCR interpretation. A BLAST hit alone does not establish a valid amplicon,
and absence of a reported primer hit does not prove biological gene absence.
PCR and XML details are audit evidence, not independent votes.

The Mash tables contain identity, shared hashes, median multiplicity, p-value,
reference ID, and an optional comment. Reference filenames encode their
phylogroup. The script selects the highest-identity hit within each group and
then the best group overall. A unique best group with Mash identity **>=0.95**
is considered supported. Equal best identities across groups are unresolved.
Mash screen identity is **not fastANI, a classification probability, or an
independent biological confirmation**. Winner-take-all Mash screen processing
can also affect the apparent separation between hits.

The paths embedded in XML database fields and Mash reference IDs are historical
identifiers. The script does not open those original database/reference paths.

**PCR calls and Mash calls are separate results.** For example, E5-5 has a D
PCR call in [ClermonTyping_metadata.tsv](../input/ClermonTyping_metadata.tsv),
but its highest Mash hit is the E-labeled `O157H7_SAKAI_E_.fasta` reference
(identity 0.984247). The metadata's phylogroup column is the PCR call; it does
not replace the Mash result referenced by the last column. This discrepancy
is retained rather than changing either input.

### Core-gene tree

[EcoZUR_fasttree.tre](../input/EcoZUR_fasttree.tre)

Only included genomes are compared.

## Methods

### 1. Construct one symmetric ANI matrix

For each pairwise genome comparison:

1. Retain the row with the largest **fifth-column** value (total query fragments corresponding to the larger genome in teh comparison)
2. If multiple rows share that maximum, choose one uniformly at random using
   the seeded random generator.
3. Place its ANI in both corresponding matrix cells. Do not average reciprocal
   ANI values.

Self-match values are retained as reported. The included dataset yields a complete 247-by-247 matrix: 30,628
pairs, including 247 self pairs.

### 2. Establish original-label reference candidates

A genome can supply a reference label only if:

- Its original phylogroup is assigned.
- Its supported Mash group agrees with that original call.
- All assigned original calls in its genomovar agree with one another.

Unknown calls do not create a genomovar conflict, but cannot supply labels.
Mixed-label genomovars are excluded from the reference set in their entirety.
The supplied inputs yield 228 reference candidates. These are operationally
trusted candidates, not an independently validated gold standard.

The reference set is built once. Predictions are never recycled as reference
labels, and no genome can use itself as a reference. Every included genome is
assessed, whether its original call was assigned or unknown.

### 3. Obtain an ANI vote

If a genome has eligible same-genomovar reference peers at **ANI >=99.5%**, and
its genomovar has no conflicting assigned original calls, their common
phylogroup becomes the **ANI vote**. At least one eligible peer is sufficient,
including in a two-member genomovar. This is not a majority vote over
conflicting original calls.

Otherwise, use eligible references **outside the entire focal genomovar**.
This external route applies to singletons, mixed-label genomovars, and cases
without eligible high-ANI peers. Examine **all** eligible external references,
not one preselected genomovar. The group of the highest-ANI reference supplies
the ANI vote only if:

- Best ANI is **strictly greater than 98.5%**.
- At least **one external genomovar** contains a candidate-group reference
  above 98.5% ANI.
- Best candidate ANI exceeds the best competing-group ANI by **at least 0.1
  percentage point**. A tiny numerical tolerance handles floating-point
  representation at this boundary.

If these criteria fail, or there is no represented competing group from which
to calculate the margin, **ANI abstains**. An abstention is recorded as
`ani_vote = Unknown`; it does not veto agreement between Tree and Mash.

One genomovar is the minimum support, not a maximum search size. The full
support count is still reported, without letting multiple genomes inflate the
distinct-genomovar count. External evidence is computed for every genome, even
when the within-genomovar route supplies the ANI vote. It is not a fourth vote
and no longer acts as a separate veto.

### 4. Obtain a Tree vote

For each genome pair, patristic distance is the sum of branch lengths along the
connecting path in the supplied Newick tree. Visual adjacency of tips in a
drawing is not used. Branch lengths must be finite and nonnegative; all included
genomes must have unique, normalized tree-tip matches.

For the within-genomovar route, tree comparisons use all eligible reference
genomes except the focal genome. For the external route, they also exclude
the entire focal genomovar.

The Tree method proposes its **own closest reference group**, independently
of the ANI candidate's identity or ANI qualification. Let `d_same` be the
distance to this group's nearest reference, and `d_other` the distance to the
nearest reference in any competing group. Tree casts a vote only if:

```text
d_same / d_other <= 0.8
```

Thus the Tree candidate reference must be at least 20% closer. Missing competing
groups, zero-distance cross-group ties, and insufficient separation do not
pass and produce `tree_vote = Unknown`. This abstention does not veto agreement
between ANI and Mash. No universal absolute tree-distance cutoff is used. Because the tree was
inferred after removing invariable sites, its branch-length scale should not
be equated with an unfiltered whole-genome substitution rate.

The method does not use bootstrap support, infer ancestral phylogroup labels,
or require monophyly. It compares distances to labeled reference tips.

### 5. Obtain a Mash vote

The focal genome's unique best Mash reference group casts a vote only at
**identity >=0.95**. Missing files, no hits, a lower identity, or an exact tie
between groups produce `mash_vote = Unknown` (an abstention).

Mash no longer vetoes a prediction. ANI and Tree can agree despite a
contradictory Mash vote. The original PCR call is not a fourth vote; it is
used for reference eligibility and the singleton fallback below.

### 6. Require two-of-three agreement, then apply the singleton fallback

Each method contributes at most one qualified vote. `Unknown` abstentions are
excluded from voting; two abstentions never count as agreement.

| Qualified agreement | Decision |
|---|---|
| ANI + Tree agree | Assign that group, regardless of Mash |
| ANI + Mash agree | Assign that group, regardless of Tree |
| Tree + Mash agree | Assign that group, regardless of ANI |
| All three agree | Assign that group |
| No group has at least two votes | Apply the fallback rules below |

No method has an unconditional veto or numerical priority. In particular,
ANI + Mash can outweigh Tree, while Tree + Mash can compensate
for an ANI abstention. Even a qualified ANI vote can lose to
agreeing Tree and Mash votes. Method quality thresholds are unchanged apart
from reducing the external ANI support minimum from two genomovars to one.

Without consensus:

1. If the focal genomovar has exactly one member and the original PCR
   phylogroup is assigned, use that original phylogroup as the final assignment.
   Record `assignment_basis = singleton_original_fallback`,
   `status = singleton_original_fallback`, and `singleton_fallback = True`.
2. Otherwise, record `predicted_phylogroup = Unknown` and `status = unresolved`.
   A singleton whose original call is unknown cannot obtain a label by fallback.

The fallback is a reporting rule, **not evidence of two-method agreement**.
For fallback rows, `consensus_phylogroup` remains `Unknown` and
`supporting_methods` is empty. Neither consensus nor fallback assignments are
fed back into the reference set.

### 7. Plot the ANI matrix

The Seaborn clustermap uses **complete linkage** on distance `100 - ANI`, with
optimal leaf ordering and the same linkage on both axes. Only clustering
self-distances are set to zero; the displayed ANI diagonal is preserved.

The displayed dendrogram is **ANI-derived**, not the supplied core-gene tree.
Three aligned annotation bars show:

1. Supplied genomovar labels and colors.
2. Original ClermonTyping phylogroups.
3. Final assessed phylogroups, including the explicitly recorded singleton
   original-call fallbacks.

Original and predicted phylogroups share a color palette; `Unknown` is gray.
Text is placed once per contiguous annotation block using black or white,
whichever provides greater contrast. The PNG is saved at **300 DPI**.
The third bar retains the label "Predicted phylogroup", and includes
fallback assignments as well as two-of-three consensus. The TSV's status/basis fields must be consulted to
distinguish a fallback from a supported consensus; they share phylogroup colors.

## Reading the TSV

[ClermonTyping_phylogroup_predictions.tsv](../results/ClermonTyping_phylogroup_predictions.tsv)
contains one row per included genome. Important fields include:

| Fields | Meaning |
|---|---|
| `genome`, `genomovar`, `genomovar_size` | Genome membership and group size |
| `original_phylogroup` | Preserved original assignment |
| `candidate_phylogroup` | ANI proposal before qualification; may differ from consensus or fallback |
| `ani_vote`, `tree_vote`, `mash_vote` | Qualified method votes; `Unknown` means abstention |
| `ani_vote_reason`, `tree_vote_reason`, `mash_vote_reason` | Qualification or abstention explanations |
| `consensus_phylogroup`, `supporting_methods` | Two-of-three result and agreeing methods; no consensus means `Unknown` and an empty method list |
| `assignment_basis`, `singleton_fallback` | Distinguish consensus, singleton original-call fallback, and unresolved cases |
| `predicted_phylogroup` | Final assignment: consensus, singleton fallback, or `Unknown` |
| `manuscript_phylogroup` | Final assignment when available; original ClermonTyping call otherwise |
| `status` | Outcome category listed below |
| `method`, `reason` | ANI/reference-pool route and final consensus/fallback explanation |
| `genomovar_conflict`, `trusted_reference` | Original-label conflict and reference eligibility |
| `peer_original_calls`, `within_supporting_neighbors` | Same-genomovar evidence |
| `decision_tree_*` | Tree evidence evaluated around the Tree method's own candidate |
| `external_*` | Best external ANI, competing group, support counts, neighbors, and tree evidence |
| `mash_*` | Best reference-group evidence and its status |
| `pcr_*`, `xml_*` | Original marker calls and primer-hit audit summaries |
| `ani_tie_seed` | Recorded random seed; 13 for launcher runs |

Outcome categories:

- `supported_original`: supported prediction agrees with the original call.
- `suggested_revision`: supported prediction differs from an assigned original.
- `predicted_unassigned`: an originally unknown genome receives a prediction.
- `singleton_original_fallback`: no two-method consensus; retain the singleton's
  assigned original call without claiming it was validated.
- `unresolved`: required evidence is absent, insufficient, or conflicting.

The positive reason `two_of_three_agreement` describes a rule-based decision,
not a confidence probability. `external_*` evidence is recorded even for genomes
using the within-genomovar route. `NA` indicates unavailable numeric evidence;
it must not be interpreted as zero. Neighbor lists contain genome IDs, group
information, ANI percentages, and tree distances. `external_tree_*` describes
support relative to the external ANI candidate, whereas `decision_tree_*`
describes the Tree method's own nearest-group candidate. These can differ.

## Singleton limitation and assignments used in the manuscript

Singletons have no within-genomovar peer and therefore depend on external
references for their ANI and Tree votes. Limited sampling or conflicting
signals can still prevent two-method consensus. The revised original-call
fallback produces a final assignment for such a singleton when an original
call exists, but does not remove this evidential limitation.

Of the 51 singletons in the current data, **50 have two-of-three consensus** and
**one (E5-5) uses the singleton original-call fallback**. Non-singletons and
originally unassigned singletons can still remain unresolved in other datasets.

For the results used in the manuscript, use the explicit
`manuscript_phylogroup` column in the evidence TSV:

Consensus revisions and newly assigned genomes use the agreed phylogroup.
Singleton fallbacks use the original call in both the PNG and manuscript column.
If a non-singleton remains unresolved, its PNG label stays `Unknown` while the
manuscript column retains its original call, preserving the earlier manuscript
reporting policy. If both calls are unknown, the manuscript value also remains
`Unknown`. Manuscript selection never changes votes or reference eligibility.

## Verified packaged run

| Outcome | Genomes |
|---|---:|
| `supported_original` | 234 |
| `suggested_revision` | 7 |
| `predicted_unassigned` | 5 |
| `singleton_original_fallback` | 1 |
| `unresolved` | 0 |

Key predictions:

- Q245-2 (g69): D to A.
- All eight g99 members: B1; the four original A calls are suggested revisions.
- Q91 and Q23-1: A to B1.
- B309-1, Q31-6, Q33-6, and Q34-6: Unknown to B1.
- Q38-1: Unknown to A.

There are 246 consensus assignments (244 with all three votes and two with
two votes) and one singleton fallback.

## Method and format references

- Beghain et al. (2018), *ClermonTyping: an easy-to-use and accurate in silico
  method for Escherichia genus strain phylotyping*.
  [doi:10.1099/mgen.0.000192](https://doi.org/10.1099/mgen.0.000192).
- [ClermonTyping documentation](https://github.com/A-BN/ClermonTyping).
- [ClermonTyping PCR interpretation](https://github.com/A-BN/ClermonTyping/blob/master/bin/clermont.py).
- [ClermonTyping Mash reporting](https://github.com/A-BN/ClermonTyping/blob/master/bin/add_mash_minimal.R).
- [Mash screen documentation](https://mash.readthedocs.io/en/latest/tutorials.html).

These sources explain input formats and tool behavior. The combined
ANI/genomovar/tree decision rules above are custom to this assessment, not an
official ClermonTyping classification method.
