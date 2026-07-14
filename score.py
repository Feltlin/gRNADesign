"""Scoring: weighted conservation (entropy) and expected CFD.

Kept as two separate columns downstream — they answer different questions
(conservation vs. mismatch tolerance) and must not be fused. See
README.md section 2.
"""

from cfd import mm_activity, calc_cfd

_ACGT = ("A", "C", "G", "T")


def weighted_entropy(window, H, w):
    """S = sum(w_i * H[spacer_cols[i]]) / sum(w_i). Weighted mean entropy in
    [0,1]; lower = more conserved where it matters (seed-weighted)."""
    num = sum(w[i] * H[window.spacer_cols[i]] for i in range(20))
    den = sum(w)
    return num / den


def expected_cfd(window, rows, freqs, mm_scores, pam_scores):
    """Returns (e_indep, e_exact), both in [0,1].

    e_indep: per-position activities combined assuming column independence.
    e_exact: mean CFD score against each variant's actual spacer+PAM, no
        independence assumption. Divergence between the two is a
        subtype-clustering signal.
    """
    g = window.spacer
    spacer_cols = window.spacer_cols

    # --- e_indep ---
    e_indep = 1.0
    for i in range(20):
        col = spacer_cols[i]
        f = freqs[col]
        a_i = sum(
            f[b] * mm_activity(i, g[i], b, mm_scores, pam_scores)
            for b in _ACGT
        )
        # gap contributes 0 via mm_activity(i, g[i], '-') = 0; ambiguity
        # codes were already excluded from freqs' denominator.
        e_indep *= a_i

    # --- e_exact ---
    pam_col_1, pam_col_2 = window.pam_cols[1], window.pam_cols[2]
    total = 0.0
    for row in rows:
        variant_spacer = "".join(row[c] for c in spacer_cols)
        if any(b not in _ACGT for b in variant_spacer):
            continue  # activity 0: gap/ambiguity in spacer window
        variant_pam = row[pam_col_1] + row[pam_col_2]
        if variant_pam not in pam_scores:
            continue  # activity 0: not a valid PAM key
        total += calc_cfd(g, variant_spacer, variant_pam, mm_scores, pam_scores)
    e_exact = total / len(rows) if rows else 0.0

    return e_indep, e_exact
