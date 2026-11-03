"""Synthetic regression tests; run python3 -B -m unittest test_plot_clermontyping_ani."""

import io
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_clermontyping_ani as script


class PredictionTests(unittest.TestCase):
    def setUp(self):
        self.genomes = pd.Index(["Q", "A1", "A2", "B1"])
        self.metadata = pd.DataFrame({
            "genomovar_complete": ["gq", "g1", "g2", "g3"],
            "genomovar_complete_colors": ["#ffffff", "#111111", "#222222", "#333333"],
        }, index=self.genomes)
        self.labels = pd.Series(["Unknown", "A", "A", "B1"], index=self.genomes)
        self.typing = pd.DataFrame({
            "mash_status": ["supported"] * 4,
            "mash_phylogroup": ["A", "A", "A", "B1"],
        }, index=self.genomes)
        self.ani = pd.DataFrame(97.0, index=self.genomes, columns=self.genomes)
        self.distances = pd.DataFrame(0.5, index=self.genomes, columns=self.genomes)
        for genome in self.genomes:
            self.ani.at[genome, genome] = 100
            self.distances.at[genome, genome] = 0
        self.set_pair(self.ani, "Q", "A1", 99.0)
        self.set_pair(self.ani, "Q", "A2", 98.8)
        self.set_pair(self.distances, "Q", "A1", 0.1)
        self.set_pair(self.distances, "Q", "A2", 0.12)

    @staticmethod
    def set_pair(matrix, a, b, value):
        matrix.at[a, b] = matrix.at[b, a] = value

    def predict(self):
        return script.predict_phylogroups(
            self.ani, self.metadata, self.labels, self.distances, self.typing
        )

    def test_external_prediction(self):
        row = self.predict().loc["Q"]
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.status, "predicted_unassigned")
        self.assertEqual(row.manuscript_phylogroup, "A")
        self.assertEqual(row.external_support_genomovars, 2)
        self.assertEqual(row.assignment_basis, "two_of_three_consensus")
        self.assertEqual(row.supporting_methods, "ANI;Tree;Mash")

    def test_within_threshold_is_inclusive(self):
        self.metadata.at["Q", "genomovar_complete"] = "g1"
        self.set_pair(self.ani, "Q", "A1", 99.5)
        row = self.predict().loc["Q"]
        self.assertEqual(row.method, "within_genomovar")
        self.assertEqual(row.predicted_phylogroup, "A")
        self.set_pair(self.ani, "Q", "A1", 99.4999)
        row = self.predict().loc["Q"]
        self.assertEqual(row.method, "external_genomovars")
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.ani_vote, "A")
        self.assertNotEqual(row.external_ani_reference, "A1")

    def test_external_threshold_is_strict(self):
        for ref in ("A1", "A2"):
            self.set_pair(self.ani, "Q", ref, 98.5)
        row = self.predict().loc["Q"]
        self.assertEqual(row.ani_vote, "Unknown")
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.supporting_methods, "Tree;Mash")
        for ref in ("A1", "A2"):
            self.set_pair(self.ani, "Q", ref, 98.5001)
        self.assertEqual(self.predict().at["Q", "predicted_phylogroup"], "A")
        self.assertEqual(self.predict().at["Q", "ani_vote"], "A")

    def test_ani_margin_and_tie(self):
        self.set_pair(self.ani, "Q", "B1", 98.9)
        self.assertEqual(self.predict().at["Q", "ani_vote"], "A")
        self.set_pair(self.ani, "Q", "B1", 98.9001)
        self.assertEqual(self.predict().at["Q", "ani_vote"], "Unknown")
        self.set_pair(self.ani, "Q", "B1", 99.0)
        self.assertEqual(self.predict().at["Q", "ani_vote"], "Unknown")
        self.assertEqual(self.predict().at["Q", "predicted_phylogroup"], "A")

    def test_tree_separation_boundary_and_conflict(self):
        self.set_pair(self.distances, "Q", "A1", 0.4)
        self.set_pair(self.distances, "Q", "A2", 0.45)
        self.assertEqual(self.predict().at["Q", "tree_vote"], "A")
        self.set_pair(self.distances, "Q", "A1", 0.4001)
        row = self.predict().loc["Q"]
        self.assertEqual(row.tree_vote, "Unknown")
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.supporting_methods, "ANI;Mash")
        self.set_pair(self.distances, "Q", "B1", 0.05)
        row = self.predict().loc["Q"]
        self.assertEqual(row.tree_vote, "B1")
        self.assertEqual(row.predicted_phylogroup, "A")

    def test_zero_tree_distance_ties_abstain(self):
        self.set_pair(self.distances, "Q", "A1", 0.0)
        self.set_pair(self.distances, "Q", "B1", 0.0)
        row = self.predict().loc["Q"]
        self.assertEqual(row.tree_vote, "Unknown")
        self.assertEqual(row.predicted_phylogroup, "A")

    def test_ani_tree_agreement_overrides_mash_conflict(self):
        self.typing.at["Q", "mash_phylogroup"] = "B1"
        row = self.predict().loc["Q"]
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.mash_vote, "B1")
        self.assertEqual(row.supporting_methods, "ANI;Tree")

    def test_missing_query_mash_is_not_a_negative_hit(self):
        self.typing.at["Q", "mash_status"] = "missing"
        self.typing.at["Q", "mash_phylogroup"] = "Unknown"
        self.assertEqual(self.predict().at["Q", "predicted_phylogroup"], "A")
        self.assertEqual(self.predict().at["Q", "mash_vote"], "Unknown")

    def test_reference_requires_mash_agreement(self):
        self.typing.at["A2", "mash_phylogroup"] = "B1"
        row = self.predict().loc["Q"]
        self.assertEqual(row.external_support_genomovars, 1)
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertFalse(self.predict().at["A2", "trusted_reference"])

    def test_genomovar_replication_does_not_inflate_support(self):
        self.metadata.at["A2", "genomovar_complete"] = "g1"
        row = self.predict().loc["Q"]
        self.assertEqual(row.external_support_genomes, 2)
        self.assertEqual(row.external_support_genomovars, 1)
        self.assertEqual(row.ani_vote, "A")
        self.assertEqual(row.predicted_phylogroup, "A")

    def test_self_never_supplies_reference_evidence(self):
        self.labels["Q"] = "A"
        self.set_pair(self.ani, "Q", "A1", 98.3)
        row = self.predict().loc["Q"]
        self.assertEqual(row.external_ani, 98.8)
        self.assertNotEqual(row.external_ani_reference, "Q")
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertNotIn("Q|", row.external_supporting_neighbors)

    def test_known_calls_are_audited_not_copied(self):
        self.labels["Q"] = "B1"
        row = self.predict().loc["Q"]
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.status, "suggested_revision")
        self.assertEqual(row.manuscript_phylogroup, "A")
        self.assertEqual(self.labels["Q"], "B1")

    def test_singleton_fallback_is_distinct_from_consensus(self):
        self.labels["Q"] = "B1"
        self.typing.at["Q", "mash_phylogroup"] = "B1"
        for ref in ("A1", "A2", "B1"):
            self.set_pair(self.distances, "Q", ref, 0.5)
        row = self.predict().loc["Q"]
        self.assertEqual(row.candidate_phylogroup, "A")
        self.assertEqual(row.predicted_phylogroup, "B1")
        self.assertEqual(row.manuscript_phylogroup, "B1")
        self.assertEqual(row.status, "singleton_original_fallback")
        self.assertTrue(row.singleton_fallback)
        self.assertEqual(row.consensus_phylogroup, "Unknown")
        self.assertEqual(row.supporting_methods, "")

    def test_manuscript_remains_unknown_when_both_calls_unknown(self):
        self.typing.at["Q", "mash_phylogroup"] = "B1"
        for ref in ("A1", "A2", "B1"):
            self.set_pair(self.distances, "Q", ref, 0.5)
        row = self.predict().loc["Q"]
        self.assertEqual(row.candidate_phylogroup, "A")
        self.assertEqual(row.predicted_phylogroup, "Unknown")
        self.assertEqual(row.original_phylogroup, "Unknown")
        self.assertEqual(row.manuscript_phylogroup, "Unknown")
        self.assertFalse(row.singleton_fallback)

    def test_manuscript_keeps_supported_original(self):
        self.labels["Q"] = "A"
        row = self.predict().loc["Q"]
        self.assertEqual(row.status, "supported_original")
        self.assertEqual(row.predicted_phylogroup, "A")
        self.assertEqual(row.manuscript_phylogroup, "A")

    def test_conflicting_genomovar_is_excluded_from_references(self):
        self.metadata.at["B1", "genomovar_complete"] = "g1"
        row = self.predict().loc["Q"]
        self.assertEqual(row.external_support_genomovars, 1)
        self.assertEqual(row.predicted_phylogroup, "Unknown")

    def test_predictions_do_not_become_references(self):
        self.labels["A1"] = "Unknown"
        self.labels["A2"] = "Unknown"
        row = self.predict().loc["Q"]
        self.assertEqual(row.predicted_phylogroup, "Unknown")
        self.assertNotIn("A1|", row.external_supporting_neighbors)
        self.assertNotIn("A2|", row.external_supporting_neighbors)

    def test_tree_mash_can_overrule_a_different_qualified_ani_vote(self):
        self.typing.at["Q", "mash_phylogroup"] = "B1"
        self.set_pair(self.distances, "Q", "B1", 0.01)
        row = self.predict().loc["Q"]
        self.assertEqual(row.ani_vote, "A")
        self.assertEqual(row.tree_vote, "B1")
        self.assertEqual(row.predicted_phylogroup, "B1")
        self.assertEqual(row.supporting_methods, "Tree;Mash")

    def test_non_singleton_has_no_original_call_fallback(self):
        self.metadata.at["Q", "genomovar_complete"] = "g1"
        self.labels["Q"] = "B1"
        self.typing.at["Q", "mash_phylogroup"] = "D"
        self.set_pair(self.distances, "Q", "B1", 0.01)
        row = self.predict().loc["Q"]
        self.assertEqual((row.ani_vote, row.tree_vote, row.mash_vote), ("A", "B1", "D"))
        self.assertEqual(row.predicted_phylogroup, "Unknown")
        self.assertFalse(row.singleton_fallback)
        self.assertEqual(row.manuscript_phylogroup, "B1")

    def test_ineligible_mash_evidence_does_not_vote(self):
        for ref in ("A1", "A2", "B1"):
            self.set_pair(self.distances, "Q", ref, 0.5)
        for status in ("missing", "no_hits", "low_identity", "tied_groups"):
            with self.subTest(status=status):
                self.typing.at["Q", "mash_status"] = status
                row = self.predict().loc["Q"]
                self.assertEqual(row.mash_vote, "Unknown")
                self.assertEqual(row.predicted_phylogroup, "Unknown")

    def test_singleton_with_no_qualified_votes_retains_original(self):
        self.labels["Q"] = "B1"
        self.typing["mash_status"] = "missing"
        result = self.predict()
        row = result.loc["Q"]
        self.assertEqual((row.ani_vote, row.tree_vote, row.mash_vote), ("Unknown",) * 3)
        self.assertEqual(row.predicted_phylogroup, "B1")
        self.assertTrue(row.singleton_fallback)
        self.assertFalse(result.trusted_reference.any())

    def test_within_group_not_rejected_for_absent_external_class(self):
        self.metadata.at["Q", "genomovar_complete"] = "g1"
        self.labels["A1"] = "C"
        self.typing.loc[["Q", "A1"], "mash_phylogroup"] = "C"
        self.set_pair(self.ani, "Q", "A1", 99.7)
        self.set_pair(self.distances, "Q", "A1", 0.01)
        self.metadata.loc["A3"] = ["g4", "#444444"]
        self.labels["A3"] = "A"
        self.typing.loc["A3"] = ["supported", "A"]
        self.ani.loc["A3"] = 97.0
        self.ani["A3"] = 97.0
        self.ani.at["A3", "A3"] = 100
        self.distances.loc["A3"] = 0.5
        self.distances["A3"] = 0.5
        self.distances.at["A3", "A3"] = 0
        self.set_pair(self.ani, "Q", "A3", 98.7)
        self.set_pair(self.distances, "Q", "A3", 0.15)
        row = self.predict().loc["Q"]
        self.assertEqual(row.external_ani_group, "A")
        self.assertEqual(row.external_support_genomovars, 2)
        self.assertEqual(self.predict().at["Q", "predicted_phylogroup"], "C")


class InputAndPlotTests(unittest.TestCase):
    def test_ani_example_duplicate_self_and_ties(self):
        ids = pd.Index(["B001-5", "E58-3"])
        text = (
            "B001-5 E58-3 99.8577 1523 1665\n"
            "E58-3 B001-5 99.9312 1529 1565\n"
            "B001-5 B001-5 100 1659 1665\n"
            "B001-5 B001-5 100 1659 1665\n"
            "E58-3 E58-3 99.9999 1565 1565\n"
            "excluded B001-5 97 1 2000\n"
        )
        with patch.object(Path, "open", return_value=io.StringIO(text)):
            result = script.load_ani_matrix(Path("synthetic.ani"), ids, random.Random(42))
        self.assertEqual(result.at["B001-5", "E58-3"], 99.8577)
        np.testing.assert_array_equal(result.values, result.values.T)
        self.assertEqual(result.at["E58-3", "E58-3"], 99.9999)
        observed = set()
        for seed in range(20):
            tied = text.replace("99.9312 1529 1565", "99.9312 1529 1665")
            with patch.object(Path, "open", return_value=io.StringIO(tied)):
                observed.add(script.load_ani_matrix(Path("synthetic.ani"), ids, random.Random(seed)).iat[0, 1])
        self.assertEqual(observed, {99.8577, 99.9312})

    def test_missing_ani_pairs_fail(self):
        with patch.object(Path, "open", return_value=io.StringIO("a a 100 1 1\n")):
            with self.assertRaisesRegex(ValueError, "missing unordered pairs"):
                script.load_ani_matrix(Path("synthetic.ani"), pd.Index(["a", "b"]), random.Random(42))

    def test_tree_names_and_branch_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tree.tre"
            path.write_text("(B001_5:0.1,E58-3:0.2,excluded:0.4);")
            distances = script.load_tree_distances(path, pd.Index(["B001-5", "E58-3"]))
            self.assertAlmostEqual(distances.iat[0, 1], 0.3)
            self.assertEqual(distances.shape, (2, 2))
            for newick in (
                "(B001_5:0.1,B001-5:0.2,E58-3:0.3);",
                "(B001_5,E58-3:0.3);",
                "(B001_5:-0.1,E58-3:0.3);",
                "(B001_5:0.1,missing:0.3);",
            ):
                path.write_text(newick)
                with self.assertRaises(ValueError):
                    script.load_tree_distances(path, pd.Index(["B001-5", "E58-3"]))

    def test_mash_is_sorted_by_identity_and_ties_abstain(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "screen.tab"
            path.write_text(
                "0.80\t10/100\t1\t0\tmash_ref/ref_A_0.fasta\t\n"
                "0.98\t90/100\t1\t0\tmash_ref/ref_B1_.fasta\t\n"
            )
            hit = script.load_mash_evidence(path)
            self.assertEqual(hit["mash_phylogroup"], "B1")
            self.assertEqual(hit["mash_status"], "supported")
            path.write_text(
                "0.98\t90/100\t1\t0\tmash_ref/ref_A_0.fasta\t\n"
                "0.98\t90/100\t1\t0\tmash_ref/ref_B1_.fasta\t\n"
            )
            self.assertEqual(script.load_mash_evidence(path)["mash_status"], "tied_groups")
            with self.assertLogs(script.LOGGER, level="WARNING"):
                self.assertEqual(script.load_mash_evidence(Path(directory) / "absent")["mash_status"], "missing")

    def test_pcr_and_xml_are_audit_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            assembly = "Q_SPAdes_assembly.fa"
            metadata = base / "metadata.tsv"
            metadata.write_text(
                f"{assembly}\t['arpA']\t['+', '-', '-', '-']\t[]\tA\t"
                f"{assembly}_mash_screen.tab\n"
            )
            (base / f"{assembly}_mash_screen.tab").write_text(
                "0.98\t90/100\t1\t0\tmash_ref/ref_A_0.fasta\t\n"
            )
            xml = base / f"{assembly}.xml"
            xml.write_text(
                f"<BlastOutput><BlastOutput_program>blastn</BlastOutput_program>"
                f"<BlastOutput_db>/db/{assembly}</BlastOutput_db>"
                "<BlastOutput_iterations><Iteration>"
                "<Iteration_query-def>arpA_F</Iteration_query-def>"
                "<Iteration_hits><Hit><Hit_hsps><Hsp/></Hit_hsps></Hit></Iteration_hits>"
                "</Iteration><Iteration><Iteration_query-def>arpA_R</Iteration_query-def>"
                "<Iteration_hits/></Iteration></BlastOutput_iterations></BlastOutput>"
            )
            evidence = script.load_typing_evidence(metadata, base, pd.Index(["Q"]))
            self.assertEqual(evidence.at["Q", "xml_queries_without_hits"], "arpA_R")
            self.assertEqual(evidence.at["Q", "pcr_products"], "arpA")
            self.assertEqual(evidence.at["Q", "pcr_arpA_chuA_yjaA_TspE4C2"], "+---")
            xml.unlink()
            with self.assertLogs(script.LOGGER, level="WARNING"):
                evidence = script.load_typing_evidence(metadata, base, pd.Index(["Q"]))
            self.assertEqual(evidence.at["Q", "xml_status"], "missing")
            self.assertEqual(evidence.at["Q", "mash_status"], "supported")

    def test_annotation_rows_share_palette_and_preserve_order(self):
        ids = pd.Index(["q", "a", "b"])
        matrix = pd.DataFrame([[100, 99.8, 97], [99.8, 100, 97.1], [97, 97.1, 100]], index=ids, columns=ids)
        metadata = pd.DataFrame({
            "genomovar_complete": ["g1", "g1", "g2"],
            "genomovar_complete_colors": ["#212121", "#212121", "#e9e9e9"],
        }, index=ids)
        original = pd.Series(["Unknown", "A", "B1"], index=ids)
        predicted = pd.Series(["A", "A", "Unknown"], index=ids)
        captured = []
        real_clustermap = script.sns.clustermap

        def capture(*args, **kwargs):
            kwargs["figsize"] = (4, 4)
            grid = real_clustermap(*args, **kwargs)
            captured.append((grid, kwargs))
            return grid

        with patch.object(script.sns, "clustermap", side_effect=capture):
            with patch("seaborn.matrix.ClusterGrid.savefig"):
                script.plot_clustermap(matrix, metadata, original, Path("unused.png"), predicted)
        grid, kwargs = captured[0]
        bars = kwargs["col_colors"]
        self.assertEqual(list(bars.columns), ["Genomovar", "ClermonTyping phylogroup", "ANI-predicted phylogroup"])
        self.assertEqual(bars.at["a", "ClermonTyping phylogroup"], bars.at["q", "ANI-predicted phylogroup"])
        order = matrix.columns[grid.dendrogram_col.reordered_ind]
        expected = predicted.loc[order].tolist()
        actual = [text.get_text() for text in grid.ax_col_colors.texts if text.get_position()[1] == 2.5]
        self.assertEqual(actual, [value for i, value in enumerate(expected) if i == 0 or value != expected[i - 1]])
        self.assertEqual(script.contrasting_text("#212121"), "white")
        self.assertEqual(script.contrasting_text("#e9e9e9"), "black")
        plt.close("all")


if __name__ == "__main__":
    unittest.main()
