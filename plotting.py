"""Diagnostic plots: conservation landscape + candidate overlay.

A separate, optional stage -- run after the main table exists (or standalone
via `python plotting.py`). Produces four figures:

  1. entropy_by_position.png       -- raw per-column Shannon entropy, both
                                       strands, over the whole alignment.
  2. entropy_window20_unweighted.png -- plain mean entropy over every
                                       20-col spacer window (a continuous,
                                       ungated version of what windows.py
                                       only keeps a filtered slice of).
  3. entropy_window20_weighted.png -- CFD-seed-weighted mean entropy over
                                       every 20-col window -- this is exactly
                                       score.weighted_entropy's S, computed
                                       continuously instead of only at the
                                       gated candidate windows.
  4. gap_fraction_window20.png     -- mean gap ('-') fraction over every
                                       20-col spacer window -- i.e. how much
                                       of the panel has a real deletion
                                       there (gap count / total sequences),
                                       the same quantity windows.py's gap
                                       gate filters on, computed continuously
                                       instead of only checked pass/fail.

Both strands are mapped onto one shared plus-strand "Start_Pos" ruler (the
same coordinate convention as the main table -- see coords.py and
README.md section 6) so the curves and the candidate overlay are directly
comparable. Candidate points (from table.py's rows) are overlaid on all
three windowed plots (2, 3, 4): marker shape encodes Derived_Strand, marker
color encodes CFD_Score_Exact.
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd

POSITION_GRID_STEP = 100

from config import CONFIG
from io_fasta import read_alignment, reverse_complement
from columns import column_frequencies, entropy
from cfd import load_cfd, cfd_position_weights

SPACER_LEN = 20
PAM_LEN = 3
WINDOW_LEN = SPACER_LEN + PAM_LEN  # 23, matches windows.enumerate_windows


def compute_profiles(fasta_path):
    """Per-column Shannon entropy and gap fraction (gap count / total
    sequences) for the plus alignment and its reverse-complement, each
    length L."""
    rows_plus = read_alignment(fasta_path)
    rows_minus = reverse_complement(rows_plus)
    L = len(rows_plus[0])

    freqs_plus, gapfrac_plus = column_frequencies(rows_plus)
    freqs_minus, gapfrac_minus = column_frequencies(rows_minus)
    H_plus = np.array(entropy(freqs_plus))
    H_minus = np.array(entropy(freqs_minus))
    return H_plus, H_minus, np.array(gapfrac_plus), np.array(gapfrac_minus), L


def rolling_window_stat(H, weights, n_starts):
    """stat[i] = dot(H[i:i+len(weights)], weights), for i in 0..n_starts-1."""
    K = len(weights)
    return np.array([np.dot(H[i:i + K], weights) for i in range(n_starts)])


def windowed_profiles(H_plus, H_minus, L, w):
    """Unweighted and CFD-weighted 20-col rolling means, both strands mapped
    onto the shared plus-strand Start_Pos ruler.

    Domain is 0..L-23 inclusive (L-22 points) -- identical to the leftmost-
    column range enumerate_windows() scans for real 23-col candidate
    windows, so x here lines up exactly with the table's Start_Pos.

    The minus-strand rolling stat is computed in its own native 5'->3' col
    order (so the weight vector sees the correct orientation), then the
    *resulting* 1-D array is reversed once to land on the plus ruler -- the
    same single-reversal argument used for coords.plus_coord.
    """
    n_starts = L - WINDOW_LEN + 1
    uniform = np.full(SPACER_LEN, 1.0 / SPACER_LEN)
    w_norm = np.array(w) / sum(w)

    plus_unweighted = rolling_window_stat(H_plus, uniform, n_starts)
    plus_weighted = rolling_window_stat(H_plus, w_norm, n_starts)

    minus_unweighted = rolling_window_stat(H_minus, uniform, n_starts)[::-1]
    minus_weighted = rolling_window_stat(H_minus, w_norm, n_starts)[::-1]

    x = np.arange(n_starts)
    return x, plus_unweighted, minus_unweighted, plus_weighted, minus_weighted


def gap_windowed_profile(gapfrac_plus, gapfrac_minus, L):
    """Mean gap ('-') fraction over every 20-col spacer window, both strands
    mapped onto the same shared Start_Pos ruler as windowed_profiles (same
    reversal argument -- see that function's docstring)."""
    n_starts = L - WINDOW_LEN + 1
    uniform = np.full(SPACER_LEN, 1.0 / SPACER_LEN)

    plus_gap = rolling_window_stat(gapfrac_plus, uniform, n_starts)
    minus_gap = rolling_window_stat(gapfrac_minus, uniform, n_starts)[::-1]
    return plus_gap, minus_gap


def _scatter_candidates(ax, candidates, y_values, cmap="viridis"):
    """y_values: dict mapping id(candidate row) -> y coordinate."""
    plus = [c for c in candidates if c["Derived_Strand"] == "+"]
    minus = [c for c in candidates if c["Derived_Strand"] == "-"]
    scores = [c["CFD_Score_Exact"] for c in candidates]
    vmin, vmax = (min(scores), max(scores)) if scores else (0.0, 1.0)

    sc = None
    if plus:
        sc = ax.scatter(
            [c["Start_Pos"] for c in plus], [y_values[id(c)] for c in plus],
            c=[c["CFD_Score_Exact"] for c in plus], cmap=cmap, vmin=vmin, vmax=vmax,
            marker="o", edgecolor="black", linewidth=0.5, s=70, label="+ strand candidate",
            zorder=5,
        )
    if minus:
        sc2 = ax.scatter(
            [c["Start_Pos"] for c in minus], [y_values[id(c)] for c in minus],
            c=[c["CFD_Score_Exact"] for c in minus], cmap=cmap, vmin=vmin, vmax=vmax,
            marker="^", edgecolor="black", linewidth=0.5, s=70, label="- strand candidate",
            zorder=5,
        )
        sc = sc or sc2
    if sc is not None:
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("CFD_Score_Exact")
    ax.legend(loc="best")


def _apply_position_grid(ax, step=POSITION_GRID_STEP):
    """Major x-axis gridline every `step` alignment/report columns, so the
    conservation landscape can be read off against a fixed ruler."""
    ax.xaxis.set_major_locator(MultipleLocator(step))
    ax.grid(True, which="major", axis="x", linestyle="--", linewidth=0.5, alpha=0.5)


def plot_entropy_by_position(H_plus, H_minus, L, out_path):
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(np.arange(L), H_plus, color="tab:blue", linewidth=0.8, label="+ strand")
    ax.plot(np.arange(L), H_minus[::-1], color="tab:orange", linewidth=0.8,
            alpha=0.7, label="- strand (mapped to plus ruler)")
    ax.set_xlabel("Alignment column (plus-strand ruler)")
    ax.set_ylabel("Shannon entropy (log base 5, range [0,1])")
    ax.set_title("Per-column Shannon entropy across the HIV-1 LTR panel")
    ax.legend(loc="upper right")
    _apply_position_grid(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_window20_unweighted(x, plus_curve, minus_curve, candidates, out_path):
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, plus_curve, color="tab:blue", linewidth=1.0, label="+ strand")
    ax.plot(x, minus_curve, color="tab:orange", linewidth=1.0, alpha=0.7, label="- strand")
    if candidates:
        curve_by_strand = {"+": plus_curve, "-": minus_curve}
        y_values = {
            id(c): curve_by_strand[c["Derived_Strand"]][int(c["Start_Pos"])]
            for c in candidates
        }
        _scatter_candidates(ax, candidates, y_values)
    ax.set_xlabel("Start_Pos (plus-strand ruler)")
    ax.set_ylabel("Mean entropy over 20-col window")
    ax.set_title("Unweighted mean Shannon entropy, every 20-col spacer window")
    _apply_position_grid(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_window20_weighted(x, plus_curve, minus_curve, candidates, out_path):
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, plus_curve, color="tab:blue", linewidth=1.0, label="+ strand")
    ax.plot(x, minus_curve, color="tab:orange", linewidth=1.0, alpha=0.7, label="- strand")
    if candidates:
        # Use the table's own Weighted_Entropy value directly (rather than
        # re-indexing the curve) so the overlay can never drift from what's
        # actually reported for that candidate.
        y_values = {id(c): c["Weighted_Entropy"] for c in candidates}
        _scatter_candidates(ax, candidates, y_values)
    ax.set_xlabel("Start_Pos (plus-strand ruler)")
    ax.set_ylabel("Weighted_Entropy (CFD seed-weighted mean entropy)")
    ax.set_title("CFD-seed-weighted Shannon entropy, every 20-col spacer window")
    _apply_position_grid(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_gap_fraction_window20(x, plus_curve, minus_curve, candidates, out_path):
    """Mean gap ('-') fraction (gap count / total sequences) over every
    20-col spacer window -- the continuous version of the gap gate
    (windows.py only keeps a pass/fail slice of this at g_max)."""
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, plus_curve, color="tab:blue", linewidth=1.0, label="+ strand")
    ax.plot(x, minus_curve, color="tab:orange", linewidth=1.0, alpha=0.7, label="- strand")
    if candidates:
        curve_by_strand = {"+": plus_curve, "-": minus_curve}
        y_values = {
            id(c): curve_by_strand[c["Derived_Strand"]][int(c["Start_Pos"])]
            for c in candidates
        }
        _scatter_candidates(ax, candidates, y_values)
    ax.set_xlabel("Start_Pos (plus-strand ruler)")
    ax.set_ylabel("Mean gap fraction over 20-col window")
    ax.set_title("Gap ('-') fraction, every 20-col spacer window")
    _apply_position_grid(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def generate_all_plots(out_dir, fasta_path, mm_pickle_path, pam_pickle_path, candidates):
    """candidates: list of dicts (or DataFrame records) with at least
    Start_Pos, Derived_Strand, Weighted_Entropy, CFD_Score_Exact. Pass []
    to plot background curves only."""
    os.makedirs(out_dir, exist_ok=True)
    print(f"{len(candidates)} candidate(s) will be overlaid (from the table)")

    print("Loading CFD tables and computing position weights ...")
    mm_scores, pam_scores = load_cfd(mm_pickle_path, pam_pickle_path)
    w = cfd_position_weights(mm_scores, pam_scores)

    print(f"Computing per-column entropy and gap-fraction profiles from {fasta_path} ...")
    H_plus, H_minus, gapfrac_plus, gapfrac_minus, L = compute_profiles(fasta_path)
    print("Computing 20-col rolling window statistics (unweighted + CFD-weighted + gap) ...")
    x, plus_u, minus_u, plus_w, minus_w = windowed_profiles(H_plus, H_minus, L, w)
    plus_gap, minus_gap = gap_windowed_profile(gapfrac_plus, gapfrac_minus, L)

    p1 = os.path.join(out_dir, "entropy_by_position.png")
    p2 = os.path.join(out_dir, "entropy_window20_unweighted.png")
    p3 = os.path.join(out_dir, "entropy_window20_weighted.png")
    p4 = os.path.join(out_dir, "gap_fraction_window20.png")

    print(f"[1/4] Drawing per-column entropy plot -> {p1}")
    plot_entropy_by_position(H_plus, H_minus, L, p1)
    print(f"[2/4] Drawing unweighted 20-col window plot -> {p2}")
    plot_window20_unweighted(x, plus_u, minus_u, candidates, p2)
    print(f"[3/4] Drawing CFD-weighted 20-col window plot -> {p3}")
    plot_window20_weighted(x, plus_w, minus_w, candidates, p3)
    print(f"[4/4] Drawing gap-fraction 20-col window plot -> {p4}")
    plot_gap_fraction_window20(x, plus_gap, minus_gap, candidates, p4)

    print(f"Wrote plots:\n  {p1}\n  {p2}\n  {p3}\n  {p4}")
    return p1, p2, p3, p4


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate entropy/candidate diagnostic plots")
    parser.add_argument("--fasta", default=CONFIG["fasta_path"])
    parser.add_argument("--mm-pickle", default=CONFIG["mm_pickle_path"])
    parser.add_argument("--pam-pickle", default=CONFIG["pam_pickle_path"])
    parser.add_argument("--table", default=CONFIG["table_path"],
                         help="Existing table.tsv to overlay as candidates")
    parser.add_argument("--out-dir", default=CONFIG["plots_dir"])
    args = parser.parse_args(argv)

    try:
        df = pd.read_csv(args.table, sep="\t")
        candidates = df.to_dict("records")
    except FileNotFoundError:
        print(f"No table found at {args.table}; plotting background curves only.")
        candidates = []

    generate_all_plots(args.out_dir, args.fasta, args.mm_pickle, args.pam_pickle, candidates)


if __name__ == "__main__":
    main()
