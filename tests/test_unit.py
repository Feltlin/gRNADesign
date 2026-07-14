import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from io_fasta import reverse_complement
from columns import column_frequencies, entropy, consensus
from windows import enumerate_windows


class TestRevcomp(unittest.TestCase):
    def test_involution(self):
        rows = ["ACGT-N", "TTTTAA", "MYRKWA"]
        rc_twice = reverse_complement(reverse_complement(rows))
        self.assertEqual(rows, rc_twice)

    def test_basic_complement(self):
        rows = ["ACGT-"]
        rc = reverse_complement(rows)
        self.assertEqual(rc, ["-ACGT"])


class TestEntropy(unittest.TestCase):
    def test_uniform_column_is_one(self):
        # 5-symbol uniform column -> entropy 1.0 (log base 5, already normalized)
        rows = ["A", "C", "G", "T", "-"]
        freqs, _ = column_frequencies(rows)
        H = entropy(freqs)
        self.assertAlmostEqual(H[0], 1.0, places=9)

    def test_pure_column_is_zero(self):
        rows = ["A", "A", "A", "A"]
        freqs, _ = column_frequencies(rows)
        H = entropy(freqs)
        self.assertAlmostEqual(H[0], 0.0, places=9)

    def test_ambiguity_shrinks_denominator_not_counted(self):
        # 'N' should be excluded entirely, not counted as any symbol
        rows = ["A", "A", "N"]
        freqs, gapfrac = column_frequencies(rows)
        self.assertAlmostEqual(freqs[0]["A"], 1.0, places=9)
        self.assertAlmostEqual(gapfrac[0], 0.0, places=9)


class TestConsensus(unittest.TestCase):
    def test_majority_acgt_ignores_gap(self):
        rows = ["A", "A", "A", "-", "-"]
        freqs, _ = column_frequencies(rows)
        self.assertEqual(consensus(freqs), "A")


class TestWindowGates(unittest.TestCase):
    def _toy_alignment(self, spacer_col_gapfrac=0.0, gg_fraction=1.0, n_rows=10):
        """Build a minimal 23-column alignment (L=23, exactly one window)."""
        L = 23
        n_gap_rows = round(spacer_col_gapfrac * n_rows)
        n_gg_rows = round(gg_fraction * n_rows)
        rows = []
        for r in range(n_rows):
            row = ["A"] * L
            if r < n_gap_rows:
                row[0] = "-"  # inject gaps into column 0 only
            row[21] = "G" if r < n_gg_rows else "A"
            row[22] = "G" if r < n_gg_rows else "A"
            rows.append("".join(row))
        return rows

    def test_gap_gate_rejects_high_gap_column(self):
        rows = self._toy_alignment(spacer_col_gapfrac=0.5, gg_fraction=1.0)
        freqs, gapfrac = column_frequencies(rows)
        cons = consensus(freqs)
        wins = enumerate_windows(cons, freqs, gapfrac, rows, "+", g_max=0.10, p_pam_min=0.80)
        self.assertEqual(wins, [])

    def test_gap_gate_passes_low_gap_column(self):
        rows = self._toy_alignment(spacer_col_gapfrac=0.0, gg_fraction=1.0)
        freqs, gapfrac = column_frequencies(rows)
        cons = consensus(freqs)
        wins = enumerate_windows(cons, freqs, gapfrac, rows, "+", g_max=0.10, p_pam_min=0.80)
        self.assertEqual(len(wins), 1)

    def test_gg_gate_rejects_low_gg_fraction(self):
        rows = self._toy_alignment(spacer_col_gapfrac=0.0, gg_fraction=0.5)
        freqs, gapfrac = column_frequencies(rows)
        cons = consensus(freqs)
        wins = enumerate_windows(cons, freqs, gapfrac, rows, "+", g_max=0.10, p_pam_min=0.80)
        self.assertEqual(wins, [])

    def test_gg_gate_passes_high_gg_fraction(self):
        rows = self._toy_alignment(spacer_col_gapfrac=0.0, gg_fraction=0.9)
        freqs, gapfrac = column_frequencies(rows)
        cons = consensus(freqs)
        wins = enumerate_windows(cons, freqs, gapfrac, rows, "+", g_max=0.10, p_pam_min=0.80)
        self.assertEqual(len(wins), 1)
        self.assertAlmostEqual(wins[0].gg_fraction, 0.9, places=9)


if __name__ == "__main__":
    unittest.main()
