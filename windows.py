"""Enumerate 23-col candidate windows and apply the gap + GG gates.

Gate order matters: these are the only two places a window can be dropped
before scoring. Everything that survives both gates gets scored.
"""

from collections import namedtuple

from coords import plus_coord

Window = namedtuple(
    "Window",
    ["spacer", "pam", "spacer_cols", "pam_cols", "strand", "start_col", "gg_fraction"],
)


def enumerate_windows(consensus, freqs, gapfrac, rows, strand, g_max, p_pam_min, L=None):
    """For each start column i in 0..L-23: window = cols i..i+22
    (spacer i..i+19; PAM N=i+20, G=i+21, G=i+22).

    1. Gap gate: keep only if every one of the 23 columns has gapfrac <= g_max.
    2. GG gate: keep only if the joint fraction of variants with G at i+21
       AND G at i+22 is >= p_pam_min.

    `start_col` on the returned Window is already mapped to the plus-strand
    report frame via coords.plus_coord.
    """
    if L is None:
        L = len(consensus)
    n_rows = len(rows)
    windows = []

    for i in range(0, L - 23 + 1):
        window_cols = range(i, i + 23)

        if any(gapfrac[c] > g_max for c in window_cols):
            continue

        g_col, g2_col = i + 21, i + 22
        gg_hits = sum(1 for row in rows if row[g_col] == "G" and row[g2_col] == "G")
        gg_fraction = gg_hits / n_rows if n_rows else 0.0
        if gg_fraction < p_pam_min:
            continue

        spacer_cols = list(range(i, i + 20))
        pam_cols = [i + 20, i + 21, i + 22]
        spacer = consensus[i:i + 20]
        pam = consensus[i + 20:i + 23]
        start_col, derived_strand = plus_coord(i, strand, L)

        windows.append(Window(
            spacer=spacer,
            pam=pam,
            spacer_cols=spacer_cols,
            pam_cols=pam_cols,
            strand=derived_strand,
            start_col=start_col,
            gg_fraction=gg_fraction,
        ))

    return windows
