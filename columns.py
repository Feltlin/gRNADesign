"""Per-column stats over an alignment: frequencies, consensus, entropy.

Vectorized with numpy; the public interface (freqs as list[dict], gapfrac
and entropy as list[float], consensus as str) is unchanged so windows.py,
score.py, and diagnostics.py don't need to know numpy is involved.
"""

import numpy as np

_ACGT = ("A", "C", "G", "T")
_SYMBOLS = ("A", "C", "G", "T", "-")
_LOG5 = np.log(5)


def column_frequencies(rows):
    """Per column, count {A,C,G,T,-}; ambiguity chars are skipped (they shrink
    the denominator rather than being counted as anything).

    Returns (freqs, gapfrac):
      freqs[i] -> dict {'A':p,'C':p,'G':p,'T':p,'-':p} normalized over counted symbols
      gapfrac[i] -> freqs[i]['-']
    """
    if not rows:
        return [], []

    arr = np.array([list(row) for row in rows])  # shape (n_rows, L)

    # a symbol not in _SYMBOLS (an ambiguity code) matches none of these
    # masks, so it never enters `total` -- exactly the "shrink the
    # denominator" rule.
    counts = np.stack([(arr == s).sum(axis=0) for s in _SYMBOLS])  # (5, L)
    total = counts.sum(axis=0)  # (L,)

    with np.errstate(divide="ignore", invalid="ignore"):
        freq_arr = np.where(total > 0, counts / np.where(total == 0, 1, total), 0.0)

    freqs = [
        {s: float(freq_arr[si, i]) for si, s in enumerate(_SYMBOLS)}
        for i in range(arr.shape[1])
    ]
    gapfrac = freq_arr[_SYMBOLS.index("-")].tolist()
    return freqs, gapfrac


def consensus(freqs):
    """Per column, argmax over {A,C,G,T} (gap/ambiguity ignored). Ties broken
    by A,C,G,T order for determinism."""
    out = []
    for f in freqs:
        best = max(_ACGT, key=lambda b: (f[b], -_ACGT.index(b)))
        out.append(best)
    return "".join(out)


def entropy(freqs):
    """Per-column Shannon entropy over all 5 symbols {A,C,G,T,-}, log base 5
    (so it's already in [0,1] -- never normalize again). 0*log(0) := 0."""
    arr = np.array([[f[s] for s in _SYMBOLS] for f in freqs])  # (L, 5)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(arr > 0, np.log(arr) / _LOG5, 0.0)
    H = -(arr * log_p).sum(axis=1)
    return H.tolist()
