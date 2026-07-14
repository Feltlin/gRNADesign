"""CFD validation -- must pass before trusting any score downstream.
See README.md sections 4.3 (equations) and 9 (why these tests matter).
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import CONFIG
from cfd import load_cfd, calc_cfd, cfd_position_weights

_ACGT = "ACGT"


class TestCFDValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mm_scores, cls.pam_scores = load_cfd(
            CONFIG["mm_pickle_path"], CONFIG["pam_pickle_path"]
        )

    def test_perfect_match_equals_pam_score(self):
        for g in ("ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT", "GATTACAGATTACAGATTAC"):
            score = calc_cfd(g, g, "GG", self.mm_scores, self.pam_scores)
            self.assertAlmostEqual(score, self.pam_scores["GG"], places=9)

    def test_gg_pam_score_near_one(self):
        print(f"pam_scores['GG'] = {self.pam_scores['GG']}")
        self.assertAlmostEqual(self.pam_scores["GG"], 1.0, places=1)

    def test_position_weights_heavier_toward_pam_proximal_seed(self):
        w = cfd_position_weights(self.mm_scores, self.pam_scores)
        print(f"cfd_position_weights() = {w}")
        # Position 20 (index 19) is PAM-proximal (the seed); position 1
        # (index 0) is PAM-distal. If indices were reversed, the seed
        # region would show the *lightest* weights instead.
        distal_mean = sum(w[:8]) / 8
        seed_mean = sum(w[12:]) / 8
        self.assertGreater(
            seed_mean, distal_mean,
            "cfd_position_weights() is not seed-heavy at the PAM-proximal end "
            "-- CFD index direction is likely reversed (see README.md section 9)",
        )


if __name__ == "__main__":
    unittest.main()
