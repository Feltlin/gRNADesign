"""Test for the reference reader, the column map and the gene lookup."""

import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import alignment as aln  # noqa: E402
import annotation  # noqa: E402

GFF = """##gff-version 3
##sequence-region chr1 1 100
chr1\tlocal\tgene\t11\t20\t.\t+\t.\tName=alpha
chr1\tlocal\tCDS\t31\t40\t.\t-\t.\tName=beta;note=ignored
chr1\tlocal\tmRNA\t1\t100\t.\t+\t.\tName=transcript
"""

BED = "chr1\t10\t20\talpha\t0\t+\nchr1\t30\t40\tbeta\t0\t-\n"

FASTA = ">chr1\n" + ("ACGT" * 25) + "\n"

# Three nested non coding feature, to check that the narrowest is reported
# first and the broadest becomes the context.
NESTED = """##gff-version 3
chr1\tlocal\tregion\t1\t100\t.\t+\t.\tName=whole
chr1\tlocal\tregion\t40\t60\t.\t+\t.\tName=middle
chr1\tlocal\tTATA_box\t45\t49\t.\t+\t.\tName=box
"""


class TestReader(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name)
        (self.path / "ref.fasta").write_text(FASTA)
        (self.path / "ann.gff3").write_text(GFF)
        (self.path / "ann.bed").write_text(BED)
        (self.path / "ann2.gff3").write_text(NESTED)

    def tearDown(self):
        self.directory.cleanup()

    def load(self, name, **kwargs):
        return annotation.load_reference(
            self.path / "ref.fasta", annotation_path=self.path / name, **kwargs)

    def test_gff_is_converted_to_zero_based_half_open(self):
        reference = self.load("ann.gff3")
        alpha = next(f for f in reference.feature if f.name == "alpha")
        self.assertEqual((alpha.start, alpha.end), (10, 20))

    def test_gff_strand_and_kind(self):
        reference = self.load("ann.gff3")
        beta = next(f for f in reference.feature if f.name == "beta")
        self.assertEqual(beta.strand, -1)
        self.assertEqual(beta.kind, "CDS")

    def test_bed_is_already_zero_based(self):
        reference = self.load("ann.bed")
        alpha = next(f for f in reference.feature if f.name == "alpha")
        self.assertEqual((alpha.start, alpha.end), (10, 20))

    def test_alias_renames_without_moving_a_coordinate(self):
        reference = self.load("ann.gff3", alias={"alpha": "renamed"})
        renamed = next(f for f in reference.feature if f.name == "renamed")
        self.assertEqual((renamed.start, renamed.end), (10, 20))
        self.assertNotIn("alpha", [f.name for f in reference.feature])

    def test_unknown_annotation_format_is_refused(self):
        (self.path / "ann.txt").write_text("whatever")
        with self.assertRaises(ValueError):
            self.load("ann.txt")

    def test_a_coding_hit_report_the_gene(self):
        reference = self.load("ann.gff3")
        gene, element = annotation.annotate(reference, 31, 40)
        self.assertEqual(gene, "beta")

    def test_a_non_coding_hit_falls_back_to_the_enclosing_region(self):
        reference = self.load("ann2.gff3")
        gene, element = annotation.annotate(reference, 46, 50)
        # The broadest overlapping feature names where we are ...
        self.assertEqual(gene, "whole")
        # ... and the narrower one name what is actually hit, narrowest first.
        self.assertEqual(element, "box; middle")

    def test_element_is_ordered_narrowest_first(self):
        reference = self.load("ann2.gff3")
        _, element = annotation.annotate(reference, 46, 47)
        self.assertEqual(element.split("; "), ["box", "middle"])

    def test_annotate_of_an_unanchored_span_is_empty(self):
        reference = self.load("ann.gff3")
        self.assertEqual(annotation.annotate(reference, None, None), ("", ""))


class TestColumnMap(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "ref.fasta"
        self.path.write_text(">chr1\nTTTTTACGTACGTACGTACGTACGTTTTTT\n")

    def tearDown(self):
        self.directory.cleanup()

    def test_a_gapped_panel_maps_onto_the_reference(self):
        row = ["ACGTACGT--ACGTACGTACGT", "ACGTACGT--ACGTACGTACGT"]
        matrix = aln.matrix_of(row)
        frequency, gap = aln.column_frequency(matrix)
        consensus = aln.consensus(frequency)
        reference = annotation.load_reference(self.path)
        mapping, note = annotation.map_column(consensus, gap, reference)

        # The two fully gapped column anchor nothing.
        self.assertTrue(np.all(mapping[8:10] == -1))
        # Everything else is strictly increasing along the reference.
        anchored = mapping[mapping >= 0]
        self.assertTrue(np.all(np.diff(anchored) > 0))
        self.assertIn("chr1", note)

    def test_orientation_is_detected_both_way(self):
        """An alignment records nothing about which strand it was deposited
        on, so a panel can arrive reverse complemented. If that went
        undetected it would anchor nowhere and every coordinate would be
        wrong, so it must be caught before anything is scored."""
        # The sequence must not be its own reverse complement, or the two
        # orientation score the same and there is nothing to detect.
        asymmetric = Path(self.directory.name) / "asym.fasta"
        asymmetric.write_text(">chr1\nTTTTTGGGAACCCATTGCAGGCATTTTTT\n")
        row = ["GGGAACCCATTGCAGGCA"] * 2
        matrix = aln.matrix_of(row)
        frequency, gap = aln.column_frequency(matrix)
        reference = annotation.load_reference(asymmetric)

        orientation, forward, reverse = annotation.detect_orientation(
            aln.consensus(frequency), gap, reference)
        self.assertEqual(orientation, "forward")
        self.assertGreater(forward, reverse)

        flipped = aln.reverse_complement(matrix)
        frequency, gap = aln.column_frequency(flipped)
        orientation, forward, reverse = annotation.detect_orientation(
            aln.consensus(frequency), gap, reference)
        self.assertEqual(orientation, "reverse")
        self.assertGreater(reverse, forward)

    def test_span_is_one_based_inclusive(self):
        mapping = np.array([-1, 4, 5, 6, -1])
        self.assertEqual(annotation.span(mapping, [1, 2, 3]), (5, 7))

    def test_span_of_an_unanchored_column_set(self):
        mapping = np.array([-1, -1, -1])
        self.assertEqual(annotation.span(mapping, [0, 1]), (None, None))


class TestRegion(unittest.TestCase):
    """A repeated element is the case --reference-region exists for.

    A local alignment has no way to prefer one copy of a duplicated region
    over another, so a panel of that region anchors to whichever scores a
    fraction higher. The region flag removes the choice.
    """

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "ref.fasta"
        # One element present twice, 100 nt apart, with filler either side.
        element = "GGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTA"
        filler = "ATATATATATCCCCAAAATTTTGGGGTTTTAAAACCCCGGGG"
        self.element, self.filler = element, filler
        self.path.write_text(f">dup\n{element}{filler}{element}{filler}\n")

    def tearDown(self):
        self.directory.cleanup()

    def test_both_copy_are_reachable_and_the_frame_stay_absolute(self):
        size = len(self.element)
        stride = size + len(self.filler)

        first = annotation.load_reference(self.path, region=(1, size))
        self.assertEqual(first.offset, 0)
        self.assertEqual(first.sequence, self.element)

        second = annotation.load_reference(
            self.path, region=(stride + 1, stride + size))
        self.assertEqual(second.offset, stride)
        self.assertEqual(second.sequence, self.element)

        # Identical content, different offset: the region is the only thing
        # that decides which copy a coordinate is reported against.
        self.assertEqual(first.sequence, second.sequence)
        self.assertNotEqual(first.offset, second.offset)
        self.assertEqual(first.full_length, second.full_length)

    def test_region_include_both_end(self):
        """START and END are both inside the region, so the length is
        END - START + 1. An off by one here silently shift every coordinate."""
        whole = annotation.load_reference(self.path).sequence
        for start, end in ((1, 10), (5, 7), (3, 3)):
            region = annotation.load_reference(self.path, region=(start, end))
            self.assertEqual(len(region.sequence), end - start + 1)
            self.assertEqual(region.sequence, whole[start - 1:end])
            self.assertEqual(region.offset, start - 1)

    def test_a_coordinate_is_reported_in_full_reference_space(self):
        size = len(self.element)
        stride = size + len(self.filler)
        second = annotation.load_reference(
            self.path, region=(stride + 1, stride + size))
        matrix = aln.matrix_of([self.element] * 2)
        frequency, gap = aln.column_frequency(matrix)
        mapping, _ = annotation.map_column(
            aln.consensus(frequency), gap, second)
        # Column 0 is base 1 of the second copy, not base 1 of the reference.
        self.assertEqual(int(mapping[0]), stride)


if __name__ == "__main__":
    unittest.main()
