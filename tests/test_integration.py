"""Integration test on the real FASTA + CFD tables."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import CONFIG
from gRNADesign import run_pipeline

_ACGT = set("ACGT")


class TestIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = run_pipeline(
            CONFIG["fasta_path"], CONFIG["mm_pickle_path"], CONFIG["pam_pickle_path"],
            CONFIG["g_max"], CONFIG["p_pam_min"],
        )

    def test_at_least_one_candidate(self):
        self.assertGreaterEqual(len(self.rows), 1)

    def test_every_candidate_passes_gg_gate(self):
        for row in self.rows:
            self.assertGreaterEqual(row["PAM_Conservation"], CONFIG["p_pam_min"])

    def test_every_candidate_is_clean_acgt_20mer(self):
        for row in self.rows:
            g = row["_spacer"]
            self.assertEqual(len(g), 20)
            self.assertTrue(set(g).issubset(_ACGT), f"guide {g!r} has non-ACGT bases")

    def test_grna_5to3_displays_spacer_then_pam_in_parens(self):
        for row in self.rows:
            spacer, display = row["_spacer"], row["gRNA_5to3"]
            self.assertTrue(display.startswith(spacer + "("))
            self.assertTrue(display.endswith(")"))
            pam = display[len(spacer) + 1:-1]
            self.assertEqual(len(pam), 3)
            self.assertTrue(set(pam).issubset(_ACGT), f"PAM {pam!r} has non-ACGT bases")

    def test_deterministic_across_runs(self):
        rows_again = run_pipeline(
            CONFIG["fasta_path"], CONFIG["mm_pickle_path"], CONFIG["pam_pickle_path"],
            CONFIG["g_max"], CONFIG["p_pam_min"],
        )
        self.assertEqual(self.rows, rows_again)


if __name__ == "__main__":
    unittest.main()
