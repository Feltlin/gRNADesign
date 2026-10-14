"""CFD (Cutting Frequency Determination), the Doench 2016 / CRISPOR mismatch model.

Anything needing a mismatch tolerance number goes through this module. Never
hand roll a mismatch table key elsewhere.

The published model is SpCas9 specific and indexed over a 20 nt spacer, so
`activity_table` returns None for any other spacer length and the CFD column
is then left empty rather than guessed.
"""

import pickle

import numpy as np

BASE = ("A", "C", "G", "T")
SPACER_LENGTH = 20

# Used only to build a mismatch table key. The 'd' side is always DNA, and
# `sg` has already been through T -> U, so U must map back to A here.
_KEY_COMPLEMENT = {"A": "T", "T": "A", "U": "A", "C": "G", "G": "C"}


def load_table(mismatch_path, pam_path):
    """Load the mismatch and PAM table, identifying each by key shape rather
    than by filename, so a swapped argument order corrects itself."""
    with open(mismatch_path, "rb") as handle:
        first = pickle.load(handle)
    with open(pam_path, "rb") as handle:
        second = pickle.load(handle)

    def is_mismatch(table):
        key = next(iter(table))
        return ":" in key and "," in key

    def is_pam(table):
        key = next(iter(table))
        return len(key) == 2 and all(character in BASE for character in key)

    if is_mismatch(first) and is_pam(second):
        return first, second
    if is_mismatch(second) and is_pam(first):
        return second, first
    raise ValueError(
        "could not tell the mismatch table from the PAM table by key shape: "
        f"first key {next(iter(first))!r}, second key {next(iter(second))!r}"
    )


def score(intended, compared, pam, mismatch_score, pam_score):
    """Canonical CFD score of `compared` against the intended protospacer.

    intended: the on-target 20 nt protospacer, 5' to 3', DNA letter.
    compared: the 20 nt being scored, same frame and length.
    pam: the 2 nt key into `pam_score` (the GG of an NGG).
    """
    intended = intended.replace("T", "U")
    compared = compared.replace("T", "U")
    result = 1.0
    for position in range(len(intended)):
        if intended[position] == compared[position]:
            continue
        key = (f"r{intended[position]}:"
               f"d{_KEY_COMPLEMENT[compared[position]]},{position + 1}")
        result *= mismatch_score[key]
    return result * pam_score[pam]


def activity(position, guide_base, target_base, mismatch_score, pam_score):
    """Activity of a single position, 0 based along the spacer, in [0, 1].

    A gap means an indel at that column, which shifts the guide out of frame,
    so activity is 0 rather than a table lookup; the CFD table has no notion
    of an indel.
    """
    if target_base == guide_base:
        return 1.0
    if target_base == "-" or guide_base == "-":
        return 0.0

    filler = "A"
    intended = [filler] * SPACER_LENGTH
    compared = [filler] * SPACER_LENGTH
    intended[position] = guide_base
    compared[position] = target_base
    raw = score("".join(intended), "".join(compared), "GG", mismatch_score, pam_score)
    return raw / pam_score["GG"]


def activity_table(mismatch_score, pam_score, spacer_length=SPACER_LENGTH):
    """(spacer_length, 4, 256) lookup of per position activity.

    table[position, guide base index, target ASCII code] is the activity. Any
    target code that is not A, C, G or T stays 0, which is exactly the
    "an indel or an ambiguity kills the guide" rule of `activity`.

    Returns None when `spacer_length` is not 20, because the published CFD
    table is only indexed over position 1 to 20.
    """
    if spacer_length != SPACER_LENGTH:
        return None
    table = np.zeros((spacer_length, len(BASE), 256))
    for position in range(spacer_length):
        for index, guide_base in enumerate(BASE):
            for target_base in BASE:
                table[position, index, ord(target_base)] = activity(
                    position, guide_base, target_base, mismatch_score, pam_score
                )
    return table


def pam_lookup(pam_score):
    """(256 * 256) lookup of the CFD PAM factor, keyed by the two PAM base.

    A key the table does not list stays 0, so a variant without a usable PAM
    contributes no cutting.
    """
    table = np.zeros(256 * 256)
    for key, value in pam_score.items():
        if len(key) == 2:
            table[ord(key[0]) * 256 + ord(key[1])] = value
    return table
