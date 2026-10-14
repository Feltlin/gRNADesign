"""Alignment input and per-column statistic.

The alignment is held as a (depth x length) uint8 matrix of ASCII code, so
every column statistic below is a vectorised numpy operation rather than a
Python loop over row.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO

#: Symbol counted in a column. Anything else (an IUPAC ambiguity code) is
#: skipped: it shrinks the denominator instead of being counted as a base.
SYMBOL = ("A", "C", "G", "T", "-")
GAP = SYMBOL.index("-")
_LOG_BASE = np.log(len(SYMBOL))

_PAIR = "AT CG GC TA UA RY YR SS WW KM MK BV VB DH HD NN -- .. nn"


def _complement_table():
    """256 entry lookup mapping an ASCII base to its complement, IUPAC aware."""
    table = np.arange(256, dtype=np.uint8)
    for pair in _PAIR.split():
        table[ord(pair[0])] = ord(pair[1])
        table[ord(pair[0].lower())] = ord(pair[1].lower())
    return table


_COMPLEMENT = _complement_table()


@dataclass(frozen=True)
class Alignment:
    """An equal length FASTA alignment."""

    identifier: tuple
    matrix: np.ndarray  # (depth, length) uint8 ASCII

    @property
    def depth(self):
        return self.matrix.shape[0]

    @property
    def length(self):
        return self.matrix.shape[1]


def matrix_of(row, source="the input"):
    """Pack a list of equal length string into the (depth, length) matrix."""
    width = {len(entry) for entry in row}
    if not row:
        raise ValueError(f"no sequence in {source}")
    if len(width) != 1:
        raise ValueError(
            f"{source} is not an alignment: row length vary ({sorted(width)[:5]})")
    packed = np.frombuffer("".join(row).encode("ascii", "replace"), dtype=np.uint8)
    return packed.reshape(len(row), width.pop()).copy()


def read_alignment(path):
    """Read an aligned FASTA file, uppercased. Every row must be equal length."""
    identifier, row = [], []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            identifier.append(record.id)
            row.append(str(record.seq).upper())
    return Alignment(tuple(identifier), matrix_of(row, path))


def reverse_complement(matrix):
    """Reverse complement every row. Gap map to gap, ambiguity code map to the
    complementary code, so a minus strand search reads the same alphabet."""
    return _COMPLEMENT[matrix][:, ::-1]


def column_frequency(matrix):
    """Per column frequency over SYMBOL.

    Returns (frequency, gap_fraction) where frequency has shape
    (length, 5) and each row sums to 1 over the counted symbol.
    """
    count = np.stack(
        [(matrix == ord(symbol)).sum(axis=0) for symbol in SYMBOL]
    ).astype(float)
    total = count.sum(axis=0)
    frequency = (count / np.where(total == 0, 1.0, total)).T
    return frequency, frequency[:, GAP].copy()


def consensus(frequency):
    """Per column majority base over A, C, G, T; gap and ambiguity ignored.

    A tie break toward A, C, G, T order, which keeps the result deterministic.
    """
    index = np.argmax(frequency[:, :4], axis=1)
    return "".join(SYMBOL[position] for position in index)


def shannon_entropy(frequency):
    r"""Per column Shannon entropy over the 5 symbol alphabet, log base 5.

    H_j = -\sum_s p_{j,s} log_5 p_{j,s}, with 0 log 0 := 0.

    Taking the log in base 5 puts H in [0, 1] already, so it is never
    normalised again downstream.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(frequency > 0, frequency * np.log(frequency) / _LOG_BASE, 0.0)
    return -term.sum(axis=1)


def running_mean(value, window):
    """Centred running mean, with the window trimmed at either end so the
    result is defined over the whole length."""
    kernel = np.ones(int(window))
    total = np.convolve(value, kernel, mode="same")
    count = np.convolve(np.ones_like(value, dtype=float), kernel, mode="same")
    return total / count


def write_fasta(sequence, path, header="consensus", width=70):
    """Write one sequence as a single record FASTA, wrapped at `width`."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    line = [sequence[i:i + width] for i in range(0, len(sequence), width)]
    path.write_text(f">{header}\n" + "\n".join(line) + "\n")
    return path
