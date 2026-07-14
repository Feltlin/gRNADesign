"""Finalist diagnostics -- run only on a shortlist, never in the main table."""

from cfd import mm_activity

_ACGT = ("A", "C", "G", "T")


def weighted_entropy_profile(window, H, w):
    """Length-20 profile of w_i * H_i, one value per spacer position -- lets
    you see *where* the conservation signal (or lack of it) sits within a
    single guide, rather than only its aggregate S."""
    return [w[i] * H[window.spacer_cols[i]] for i in range(20)]


def per_variant_distance(window, rows, mm_scores, pam_scores):
    """CFD-weighted Hamming distance from each of the panel's variants to
    this window's consensus spacer. Distance at a position is
    1 - mm_activity(...) (0 if the variant matches consensus there, up to 1
    for a fully CFD-intolerant mismatch); a gap/ambiguity at that position
    is treated as maximal distance (1.0), same convention as expected_cfd's
    activity-0 rule for indels.

    Returns a list of length len(rows); summarize as median/90th percentile,
    or cluster, to see B sub-lineages.
    """
    g = window.spacer
    spacer_cols = window.spacer_cols
    distances = []
    for row in rows:
        d = 0.0
        for i in range(20):
            b = row[spacer_cols[i]]
            if b not in _ACGT:
                d += 1.0
                continue
            d += 1.0 - mm_activity(i, g[i], b, mm_scores, pam_scores)
        distances.append(d)
    return distances


def summarize_distances(distances):
    """median + 90th percentile of a per_variant_distance list."""
    if not distances:
        return {"median": None, "p90": None}
    s = sorted(distances)
    n = len(s)

    def percentile(p):
        idx = min(n - 1, int(round(p * (n - 1))))
        return s[idx]

    return {"median": percentile(0.5), "p90": percentile(0.9)}
