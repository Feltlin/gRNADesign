"""Single report-frame ruler: plus-strand alignment column.

Scoring frame (guide 5'->3', CFD-indexed) and report frame (plus-strand
column) are kept deliberately separate — see README.md section 6.
"""


def plus_coord(start_col, strand, L):
    """Map a window's start column to the plus-strand report frame.

    start_col: leftmost column of the 23-col window, in the coordinate
        system it was found in (plus alignment if strand == '+', the
        reverse-complemented alignment if strand == '-').
    strand: '+' or '-' (which alignment the window was found in).
    L: alignment length (866).

    Returns (start_col_plus, derived_strand).

    derived_strand semantics: '+' = spacer read off the plus strand, PAM is
    the NGG on the plus strand, guide RNA binds the minus strand. This is a
    "derived-from" label, never a "binds-to" label.
    """
    if strand == "+":
        return start_col, "+"
    if strand == "-":
        return L - (start_col + 23), "-"
    raise ValueError(f"strand must be '+' or '-', got {strand!r}")
