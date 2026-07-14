"""FASTA alignment loading & reverse-complementing."""

from Bio import SeqIO
from Bio.Seq import Seq


def read_alignment(path: str) -> list[str]:
    """
    Read a FASTA alignment. Returns list[str] rows, uppercased, all equal length.

    Sequences must be aligned (equal length).
    """
    with open(path) as fin:
        rows = [str(sequence.seq).upper() for sequence in SeqIO.parse(fin, "fasta")]

    lengths = {len(r) for r in rows}
    if len(lengths) != 1:
        raise ValueError(f"Alignment rows are not equal length: {lengths}")

    return rows


def reverse_complement(rows: list[str]) -> list[str]:
    """
    Reverse-complement every row. Never read a row backwards elsewhere.

    Bio.Seq.reverse_complement() correctly maps gaps to gaps and IUPAC
    ambiguity codes to their complementary code (M<->K, R<->Y, W<->W, N<->N).
    """
    return [str(Seq(sequence).reverse_complement()) for sequence in rows]
