"""IUPAC PAM motif, and enumeration of the candidate window that survive the gate.

Every gate is evaluated for all start position at once as a numpy operation
over the (depth x length) alignment matrix.
"""

from dataclasses import dataclass

import numpy as np

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def motif_mask(motif):
    """(len(motif), 256) boolean lookup: mask[i, code] is True when the ASCII
    base `code` satisfies position i of the motif."""
    mask = np.zeros((len(motif), 256), dtype=bool)
    for position, symbol in enumerate(motif.upper()):
        if symbol not in IUPAC:
            raise ValueError(f"{symbol!r} in PAM {motif!r} is not an IUPAC code")
        for base in IUPAC[symbol]:
            mask[position, ord(base)] = True
    return mask


def matches(sequence, motif):
    """True when `sequence` satisfies the IUPAC `motif` at every position."""
    if len(sequence) != len(motif):
        return False
    return all(
        base.upper() in IUPAC[symbol]
        for base, symbol in zip(sequence, motif.upper())
    )


@dataclass(frozen=True)
class Candidate:
    """One surviving window.

    `spacer_column` and `pam_column` index the alignment the window was found
    in (the plus matrix, or the reverse complemented matrix for a minus strand
    hit), so they are what scoring read. `plus_spacer_column` and
    `plus_pam_column` index the plus alignment in ascending order, so they are
    what reporting and plotting read.
    """

    spacer: str
    pam: str
    strand: str
    spacer_column: np.ndarray
    pam_column: np.ndarray
    plus_spacer_column: np.ndarray
    plus_pam_column: np.ndarray
    pam_conservation: float
    pam_fraction: float

    @property
    def plus_start(self):
        """Leftmost plus strand alignment column of the protospacer.

        The protospacer, not the whole window: it is the same span that
        Ref_Start and Ref_End report, so the alignment and reference
        coordinate always describe the same thing.
        """
        return int(self.plus_spacer_column[0])

    @property
    def plus_end(self):
        """Rightmost plus strand alignment column of the protospacer."""
        return int(self.plus_spacer_column[-1])


def _offset(pam_side, spacer_length, pam_length):
    """Offset of the spacer and of the PAM from the start of the window."""
    if pam_side == "3prime":
        return 0, spacer_length
    if pam_side == "5prime":
        return pam_length, 0
    raise ValueError(f"pam_side must be '3prime' or '5prime', got {pam_side!r}")


def enumerate_candidate(
    matrix,
    gap_fraction,
    consensus_sequence,
    motif,
    strand,
    *,
    spacer_length,
    pam_side,
    max_gap_fraction,
    min_pam_fraction,
):
    """Slide a (spacer_length + len(motif)) window and keep the window passing
    every gate.

    Gate, in order:

    1. Gap. Every column of the window has a gap fraction at or below
       `max_gap_fraction`.
    2. Consensus PAM. The consensus base at the PAM column satisfy the motif.
    3. PAM presence. At least `min_pam_fraction` of the aligned sequence carry
       a motif matching PAM at that locus, so the site is cuttable across the
       panel rather than only in the consensus.

    `pam_conservation` is reported separately from gate 3: it is the fraction
    of the aligned sequence whose PAM equals the consensus PAM exactly.
    """
    length = matrix.shape[1]
    pam_length = len(motif)
    window = spacer_length + pam_length
    span = length - window + 1
    if span <= 0:
        return []

    spacer_offset, pam_offset = _offset(pam_side, spacer_length, pam_length)
    mask = motif_mask(motif)
    code = np.frombuffer(consensus_sequence.encode("ascii"), dtype=np.uint8)

    # Gate 1: sliding maximum of the gap fraction over the window.
    gap_window = np.lib.stride_tricks.sliding_window_view(gap_fraction, window)
    gap_ok = gap_window.max(axis=1) <= max_gap_fraction

    # Column `start + pam_offset + k` for every start at once is the slice
    # [pam_offset + k : pam_offset + k + span].
    consensus_ok = np.ones(span, dtype=bool)
    functional = np.ones((matrix.shape[0], span), dtype=bool)
    identical = np.ones((matrix.shape[0], span), dtype=bool)
    for k in range(pam_length):
        column = slice(pam_offset + k, pam_offset + k + span)
        base = matrix[:, column]
        consensus_ok &= mask[k][code[column]]
        functional &= mask[k][base]
        identical &= base == code[column]

    pam_fraction = functional.mean(axis=0)
    pam_conservation = identical.mean(axis=0)

    keep = np.flatnonzero(gap_ok & consensus_ok & (pam_fraction >= min_pam_fraction))

    candidate = []
    for start in keep:
        spacer_column = np.arange(start + spacer_offset,
                                  start + spacer_offset + spacer_length)
        pam_column = np.arange(start + pam_offset, start + pam_offset + pam_length)
        if strand == "+":
            plus_spacer, plus_pam = spacer_column, pam_column
        else:
            # Reverse complemented column j is plus column length - 1 - j.
            plus_spacer = (length - 1 - spacer_column)[::-1]
            plus_pam = (length - 1 - pam_column)[::-1]
        candidate.append(Candidate(
            spacer=consensus_sequence[spacer_column[0]:spacer_column[-1] + 1],
            pam=consensus_sequence[pam_column[0]:pam_column[-1] + 1],
            strand=strand,
            spacer_column=spacer_column,
            pam_column=pam_column,
            plus_spacer_column=plus_spacer,
            plus_pam_column=plus_pam,
            pam_conservation=float(pam_conservation[start]),
            pam_fraction=float(pam_fraction[start]),
        ))
    return candidate
