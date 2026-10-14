"""End to end test on a synthetic panel the test builds itself.

No network is touched: the MFE stage is exercised through its parser and its
cache, never by submitting a job.
"""

import os
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mfe  # noqa: E402
import report  # noqa: E402
from config import CONFIG  # noqa: E402
from gRNADesign import design  # noqa: E402
from test import synthetic  # noqa: E402

ACGT = set("ACGT")


def run(directory, **override):
    reference = synthetic.reference_sequence()
    panel_path, reference_path = synthetic.write(
        directory, reference, synthetic.panel(reference))
    (directory / "ann.gff3").write_text(synthetic.annotation_gff3(len(reference)))
    argument = dict(
        alignment_path=panel_path,
        reference_path=reference_path,
        annotation_path=directory / "ann.gff3",
        mismatch_table_path=CONFIG["mismatch_table_path"],
        pam_table_path=CONFIG["pam_table_path"],
        pam="NGG", pam_side="3prime", spacer_length=20,
        max_gap_fraction=0.10, min_pam_fraction=0.80,
        reference_region=None, core_gap_fraction=0.50,
        log=lambda *_: None,
    )
    argument.update(override)
    alignment_path = argument.pop("alignment_path")
    reference_path = argument.pop("reference_path")
    return design(alignment_path, reference_path, **argument)


class TestPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.row, cls.context = run(Path(cls.directory.name))

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def test_at_least_one_candidate(self):
        self.assertGreaterEqual(len(self.row), 1)

    def test_every_guide_is_a_clean_acgt_spacer(self):
        for entry in self.row:
            self.assertEqual(len(entry["_spacer"]), 20)
            self.assertTrue(set(entry["_spacer"]) <= ACGT, entry["_spacer"])

    def test_the_display_string_is_spacer_then_bracketed_pam(self):
        """gRNA5to3 must read SPACER(PAM): the bracket mark the PAM as the
        motif that has to be in the target, not part of what you order."""
        for entry in self.row:
            shown = entry["gRNA5to3"]
            spacer, _, rest = shown.partition("(")
            self.assertEqual(spacer, entry["_spacer"])
            self.assertTrue(rest.endswith(")"), shown)
            pam = rest[:-1]
            self.assertEqual(len(pam), 3)
            self.assertTrue(pam.endswith("GG"), pam)

    def test_every_score_is_in_range(self):
        for entry in self.row:
            for column in ("PAM_Conservation", "CFD_Score", "Entropy_Mean"):
                self.assertGreaterEqual(entry[column], 0.0, column)
                self.assertLessEqual(entry[column], 1.0, column)

    def test_the_reference_span_is_at_least_the_spacer(self):
        for entry in self.row:
            if entry["Ref_Start"] is None:
                continue
            self.assertGreaterEqual(
                entry["Ref_End"] - entry["Ref_Start"] + 1, 20)

    def test_annotation_is_attached(self):
        hit = [entry for entry in self.row if entry["Gene"]]
        self.assertTrue(hit, "no candidate picked up the synthetic gene")

    def test_alignment_and_reference_span_describe_the_same_thing(self):
        """Align_Start/End and Ref_Start/End must both cover the protospacer,
        not the whole spacer plus PAM window, or the two coordinate system
        quietly disagree on a minus strand guide."""
        for entry in self.row:
            self.assertEqual(entry["Align_End"] - entry["Align_Start"] + 1, 20)
            if entry["Ref_Start"] is not None:
                self.assertGreaterEqual(
                    entry["Ref_End"] - entry["Ref_Start"] + 1, 20)

    def test_coordinate_ascend_on_both_strand(self):
        """Genomic convention: start <= end whatever the strand. Orientation
        is carried by Derived_Strand, never by reversing the coordinate."""
        for entry in self.row:
            self.assertLess(entry["Align_Start"], entry["Align_End"])
            if entry["Ref_Start"] is not None:
                self.assertLess(entry["Ref_Start"], entry["Ref_End"])

    def test_strand_label_is_plus_or_minus(self):
        self.assertTrue({e["Derived_Strand"] for e in self.row} <= {"+", "-"})

    def test_the_table_round_trips(self):
        frame = report.to_frame(self.row)
        self.assertEqual(list(frame.columns), list(report.COLUMN))
        self.assertEqual(len(frame), len(self.row))

    def test_deterministic_across_run(self):
        with tempfile.TemporaryDirectory() as directory:
            again, _ = run(Path(directory))
        self.assertEqual(self.row, again)


class TestOrientation(unittest.TestCase):
    def test_a_reverse_complemented_panel_gives_the_same_answer(self):
        """A FASTA records nothing about which strand it holds, so a panel can
        arrive flipped. It must come back to the same guide, not to garbage."""
        from Bio.Seq import Seq

        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            forward, _ = run(directory)

            reference = synthetic.reference_sequence()
            flipped = [str(Seq(entry).reverse_complement())
                       for entry in synthetic.panel(reference)]
            synthetic.write(directory, reference, flipped, name="rc.fasta")
            reverse, _ = run(directory, alignment_path=directory / "rc.fasta")

        self.assertEqual([e["gRNA5to3"] for e in forward],
                         [e["gRNA5to3"] for e in reverse])
        self.assertEqual([e["Ref_Start"] for e in forward],
                         [e["Ref_Start"] for e in reverse])


class TestWithoutReference(unittest.TestCase):
    def test_no_reference_still_score_every_candidate(self):
        """The reference only supplies coordinate and annotation. Every score
        comes from the panel, so they must survive its absence."""
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            with_reference, _ = run(directory)
            without, _ = run(directory, reference_path=None,
                             annotation_path=None)

        self.assertEqual([e["gRNA5to3"] for e in with_reference],
                         [e["gRNA5to3"] for e in without])
        self.assertEqual([e["CFD_Score"] for e in with_reference],
                         [e["CFD_Score"] for e in without])
        for entry in without:
            self.assertIsNone(entry["Ref_Start"])
            self.assertEqual(entry["Gene"], "")


class TestUnalignedInput(unittest.TestCase):
    """Equal row length is necessary for an alignment but not sufficient."""

    def test_varying_length_is_refused(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "ragged.fasta"
            path.write_text(">a\nACGTACGT\n>b\nACGT\n")
            with self.assertRaises(ValueError) as caught:
                run(Path(name), alignment_path=path)
            self.assertIn("not an alignment", str(caught.exception))

    def test_equal_length_but_unaligned_is_warned_about(self):
        """Sequence merely trimmed to a common length pass the length check
        and then score nonsense, because the column do not correspond. That
        has to be surfaced rather than silently scored."""
        rng = random.Random(5)
        row = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(20)]
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            path = directory / "unaligned.fasta"
            path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(row)))
            said = []
            try:
                run(directory, alignment_path=path, log=said.append)
            except Exception:  # noqa: BLE001 - the warning is the point
                pass
        self.assertTrue(any("not actually aligned" in line for line in said),
                        f"no warning raised; log was {said}")

    def test_a_properly_aligned_gapless_panel_is_not_warned_about(self):
        """A conserved panel with no indel is legitimate and common, so the
        absence of gap must not on its own trigger the warning."""
        said = []
        with tempfile.TemporaryDirectory() as name:
            run(Path(name), log=said.append)
        self.assertFalse(any("not actually aligned" in line for line in said))


class TestFoldParser(unittest.TestCase):
    """The RNAfold result format, parsed without touching the network."""

    def parse(self, text):
        return mfe._RESULT_PATTERN.search(text)

    def test_reads_the_structure_and_the_energy(self):
        match = self.parse(
            ">guide\nCUGUGGGUCUACCACACACA\n.(((((.....))))).... ( -6.70)\n")
        self.assertEqual(match.group(1), ".(((((.....)))))....")
        self.assertAlmostEqual(float(match.group(2)), -6.70, 6)

    def test_reads_an_unstructured_guide(self):
        match = self.parse(">guide\nAAAA\n.... (  0.00)\n")
        self.assertEqual(match.group(1), "....")
        self.assertAlmostEqual(float(match.group(2)), 0.0, 6)

    def test_reads_a_two_digit_energy(self):
        match = self.parse("((((....)))) (-12.90)\n")
        self.assertAlmostEqual(float(match.group(2)), -12.90, 6)


class TestFoldCache(unittest.TestCase):
    def test_a_cached_guide_never_reaches_a_backend(self):
        """The cache key is sequence|scaffold|temperature, so a hit must be
        served without the server being contacted at all."""
        import json

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory) / "cache.json"
            cache.write_text(json.dumps({"ACGT||37": ["....", -1.5]}))
            result = mfe.fold(
                ["ACGT"], backend="server", url="http://invalid.invalid",
                temperature=37.0, cache_path=cache, log=lambda *_: None)
        self.assertEqual(result["ACGT"], ("....", -1.5))

    def test_a_different_temperature_is_a_different_key(self):
        import json

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory) / "cache.json"
            cache.write_text(json.dumps({"ACGT||37": ["....", -1.5]}))
            result = mfe.fold(
                ["ACGT"], backend="server", url="http://invalid.invalid",
                temperature=25.0, cache_path=cache, log=lambda *_: None)
        self.assertEqual(result["ACGT"], mfe.FAILED)


if __name__ == "__main__":
    unittest.main()
