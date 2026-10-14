"""Unit test for the alignment statistic, the PAM gate and the CFD model."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import alignment as aln  # noqa: E402
import cfd  # noqa: E402
import score as score_module  # noqa: E402
from config import CONFIG  # noqa: E402
from pam import enumerate_candidate, matches, motif_mask  # noqa: E402


def decode(matrix):
    return ["".join(chr(code) for code in row) for row in matrix]


class TestReverseComplement(unittest.TestCase):
    def test_involution(self):
        row = ["ACGT-N", "TTTTAA", "MYRKWA"]
        matrix = aln.matrix_of(row)
        twice = aln.reverse_complement(aln.reverse_complement(matrix))
        self.assertEqual(decode(twice), row)

    def test_gap_maps_to_gap(self):
        matrix = aln.reverse_complement(aln.matrix_of(["ACGT-"]))
        self.assertEqual(decode(matrix), ["-ACGT"])

    def test_ambiguity_maps_to_its_complement(self):
        matrix = aln.reverse_complement(aln.matrix_of(["RYKMWS"]))
        self.assertEqual(decode(matrix), ["SWKMRY"])


class TestColumnStatistic(unittest.TestCase):
    def entropy_of(self, row):
        frequency, _ = aln.column_frequency(aln.matrix_of(row))
        return aln.shannon_entropy(frequency)

    def test_uniform_column_is_one(self):
        # Five symbol uniform column, log base 5, so exactly 1.
        self.assertAlmostEqual(self.entropy_of(["A", "C", "G", "T", "-"])[0], 1.0, 9)

    def test_pure_column_is_zero(self):
        self.assertAlmostEqual(self.entropy_of(["A", "A", "A", "A"])[0], 0.0, 9)

    def test_ambiguity_shrinks_the_denominator(self):
        frequency, gap = aln.column_frequency(aln.matrix_of(["A", "A", "N"]))
        self.assertAlmostEqual(frequency[0][0], 1.0, 9)
        self.assertAlmostEqual(gap[0], 0.0, 9)

    def test_consensus_ignores_gap(self):
        frequency, _ = aln.column_frequency(aln.matrix_of(["A", "A", "A", "-", "-"]))
        self.assertEqual(aln.consensus(frequency), "A")

    def test_entropy_is_strand_symmetric(self):
        row = ["ACGTACGTAC", "ACGTACGTAG", "A-GTACGTAC"]
        forward = self.entropy_of(row)
        reverse = aln.shannon_entropy(
            aln.column_frequency(aln.reverse_complement(aln.matrix_of(row)))[0])
        np.testing.assert_allclose(forward, reverse[::-1], atol=1e-12)

    def test_running_mean_of_a_constant_is_that_constant(self):
        value = np.full(50, 0.3)
        np.testing.assert_allclose(aln.running_mean(value, 20), value, atol=1e-12)


class TestMotif(unittest.TestCase):
    def test_iupac_match(self):
        self.assertTrue(matches("AGG", "NGG"))
        self.assertTrue(matches("TGG", "NGG"))
        self.assertFalse(matches("AGA", "NGG"))
        self.assertTrue(matches("TTTA", "TTTV"))
        self.assertFalse(matches("TTTT", "TTTV"))

    def test_mask_rejects_an_unknown_code(self):
        with self.assertRaises(ValueError):
            motif_mask("NGZ")


class TestGate(unittest.TestCase):
    """A 23 column toy alignment holds exactly one window."""

    def toy(self, gap_row=0, pam_row=10, depth=10, length=23):
        row = []
        for index in range(depth):
            base = ["A"] * length
            if index < gap_row:
                base[0] = "-"
            base[21] = base[22] = "G" if index < pam_row else "A"
            row.append("".join(base))
        return aln.matrix_of(row)

    def run_gate(self, matrix, max_gap_fraction=0.10, min_pam_fraction=0.80):
        frequency, gap = aln.column_frequency(matrix)
        return enumerate_candidate(
            matrix, gap, aln.consensus(frequency), "NGG", "+",
            spacer_length=20, pam_side="3prime",
            max_gap_fraction=max_gap_fraction, min_pam_fraction=min_pam_fraction)

    def test_gap_gate_reject(self):
        self.assertEqual(self.run_gate(self.toy(gap_row=5)), [])

    def test_gap_gate_pass(self):
        self.assertEqual(len(self.run_gate(self.toy())), 1)

    def test_pam_fraction_gate_reject(self):
        self.assertEqual(self.run_gate(self.toy(pam_row=5)), [])

    def test_pam_fraction_gate_pass(self):
        found = self.run_gate(self.toy(pam_row=9))
        self.assertEqual(len(found), 1)
        self.assertAlmostEqual(found[0].pam_fraction, 0.9, 9)

    def test_pam_conservation_counts_the_exact_consensus_pam(self):
        # Half the panel reads AGG, half TGG; the consensus PAM is one of
        # them, so exact conservation is 0.5 while motif presence is 1.0.
        row = []
        for index in range(10):
            base = ["A"] * 23
            base[20] = "A" if index < 5 else "T"
            base[21] = base[22] = "G"
            row.append("".join(base))
        found = self.run_gate(aln.matrix_of(row))
        self.assertEqual(len(found), 1)
        self.assertAlmostEqual(found[0].pam_fraction, 1.0, 9)
        self.assertAlmostEqual(found[0].pam_conservation, 0.5, 9)

    def test_minus_strand_column_map_back_to_the_plus_frame(self):
        matrix = self.toy()
        reverse = aln.reverse_complement(matrix)
        frequency, gap = aln.column_frequency(reverse)
        found = enumerate_candidate(
            reverse, gap, aln.consensus(frequency), "NGG", "-",
            spacer_length=20, pam_side="3prime",
            max_gap_fraction=1.0, min_pam_fraction=0.0)
        for entry in found:
            self.assertTrue(np.all(np.diff(entry.plus_spacer_column) == 1))
            self.assertTrue(entry.plus_spacer_column.min() >= 0)
            self.assertTrue(entry.plus_spacer_column.max() < matrix.shape[1])


class TestCFD(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mismatch, cls.pam = cfd.load_table(
            CONFIG["mismatch_table_path"], CONFIG["pam_table_path"])

    def test_perfect_match_equals_the_pam_score(self):
        for guide in ("ACGTACGTACGTACGTACGT", "GATTACAGATTACAGATTAC"):
            self.assertAlmostEqual(
                cfd.score(guide, guide, "GG", self.mismatch, self.pam),
                self.pam["GG"], 9)

    def test_gg_pam_score_is_near_one(self):
        self.assertAlmostEqual(self.pam["GG"], 1.0, 1)

    def test_index_direction_is_seed_heavy(self):
        """Position 20 is PAM proximal. If the index ran the other way the
        seed would look like the most tolerant end, not the harshest."""
        penalty = []
        for position in range(20):
            pair = [(x, y) for x in cfd.BASE for y in cfd.BASE if x != y]
            mean = sum(cfd.activity(position, x, y, self.mismatch, self.pam)
                       for x, y in pair) / len(pair)
            penalty.append(1.0 - mean)
        self.assertGreater(sum(penalty[12:]) / 8, sum(penalty[:8]) / 8)

    def test_activity_table_agrees_with_activity(self):
        table = cfd.activity_table(self.mismatch, self.pam)
        for position in (0, 7, 19):
            for index, guide in enumerate(cfd.BASE):
                for target in cfd.BASE:
                    self.assertAlmostEqual(
                        table[position, index, ord(target)],
                        cfd.activity(position, guide, target,
                                     self.mismatch, self.pam), 12)

    def test_gap_and_ambiguity_score_zero(self):
        table = cfd.activity_table(self.mismatch, self.pam)
        for code in ("-", "N", "R"):
            self.assertEqual(table[0, 0, ord(code)], 0.0)

    def test_activity_table_is_none_outside_the_published_length(self):
        self.assertIsNone(cfd.activity_table(self.mismatch, self.pam, 23))


class TestScore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mismatch, cls.pam = cfd.load_table(
            CONFIG["mismatch_table_path"], CONFIG["pam_table_path"])
        cls.table = cfd.activity_table(cls.mismatch, cls.pam)

    def candidate(self, row):
        matrix = aln.matrix_of(row)
        frequency, gap = aln.column_frequency(matrix)
        found = enumerate_candidate(
            matrix, gap, aln.consensus(frequency), "NGG", "+",
            spacer_length=20, pam_side="3prime",
            max_gap_fraction=1.0, min_pam_fraction=0.0)
        return matrix, found[0]

    def test_identical_panel_scores_the_pam_score(self):
        row = ["ACGTACGTACGTACGTACGTAGG"] * 8
        matrix, entry = self.candidate(row)
        self.assertAlmostEqual(
            score_module.cfd_mean(matrix, entry, "NGG", self.table, self.pam),
            self.pam["GG"], 9)

    def test_a_variant_with_a_broken_pam_contributes_zero(self):
        # The majority must keep a motif matching PAM, otherwise the consensus
        # PAM itself fails the gate and there is no candidate to score.
        row = ["ACGTACGTACGTACGTACGTAGG"] * 7 + ["ACGTACGTACGTACGTACGTAAA"] * 3
        matrix, entry = self.candidate(row)
        self.assertEqual(self.pam["AA"], 0.0)
        self.assertAlmostEqual(
            score_module.cfd_mean(matrix, entry, "NGG", self.table, self.pam),
            self.pam["GG"] * 0.7, 9)

    def test_a_gap_in_the_protospacer_contributes_zero(self):
        row = ["ACGTACGTACGTACGTACGTAGG"] * 7 + ["ACGTACGTACG-ACGTACGTAGG"] * 3
        matrix, entry = self.candidate(row)
        self.assertAlmostEqual(
            score_module.cfd_mean(matrix, entry, "NGG", self.table, self.pam),
            self.pam["GG"] * 0.7, 9)

    def test_entropy_mean_is_a_plain_mean(self):
        entropy = np.arange(30, dtype=float)
        column = np.arange(5, 25)
        self.assertAlmostEqual(
            score_module.entropy_mean(entropy, column), float(np.mean(column)), 12)


if __name__ == "__main__":
    unittest.main()
