"""Assemble and write the main output table.

Column names are the user-facing report format -- see README.md section 4
for what each one means. Internal code (Window, score.py) keeps its
original short/technical field names; only this module translates them
into the reader-facing table.
"""

import os

import pandas as pd

FIELDS = [
    "gRNA_5to3", "Start_Pos", "Derived_Strand", "PAM_Conservation",
    "Weighted_Entropy", "CFD_Score_Independent", "CFD_Score_Exact",
    "Offtarget_Count", "Offtarget_Repeat_Flag",
]


def assemble_table(windows, H, w, rows, freqs, mm_scores, pam_scores):
    """One row per surviving window. Offtarget_Count / Offtarget_Repeat_Flag
    are left as <NA> until the optional Cas-OFFinder step
    (offtarget.score_offinder + table.merge_offtarget) has run.

    gRNA_5to3 is a *display* string ("SPACER(PAM)") -- the PAM is shown in
    parentheses so it's visually obvious, immediately 3' of the spacer in
    the guide's own 5'->3' reading frame (true regardless of strand -- see
    README.md section 6's coordinate-system walkthrough). The row also carries a
    "_spacer" key with the bare 20-nt spacer used internally for off-target
    matching; it's dropped automatically when the table is written (it's
    not in table.FIELDS).
    """
    from score import weighted_entropy, expected_cfd

    out_rows = []
    for window in windows:
        S = weighted_entropy(window, H, w)
        e_indep, e_exact = expected_cfd(window, rows, freqs, mm_scores, pam_scores)
        out_rows.append({
            "gRNA_5to3": f"{window.spacer}({window.pam})",
            "_spacer": window.spacer,
            "Start_Pos": window.start_col,
            "Derived_Strand": window.strand,
            "PAM_Conservation": window.gg_fraction,
            "Weighted_Entropy": S,
            "CFD_Score_Independent": e_indep,
            "CFD_Score_Exact": e_exact,
            "Offtarget_Count": pd.NA,
            "Offtarget_Repeat_Flag": pd.NA,
        })
    return out_rows


def merge_offtarget(out_rows, offtarget_summary):
    """Fill in Offtarget_Count / Offtarget_Repeat_Flag from
    offtarget.score_offinder()'s per-guide summary, keyed by the bare
    spacer (row["_spacer"]), matching what run_offtarget() passed in as
    guides."""
    for row in out_rows:
        summary = offtarget_summary.get(row["_spacer"])
        if summary is None:
            continue
        row["Offtarget_Count"] = summary["offtarget_count"]
        row["Offtarget_Repeat_Flag"] = summary["n_repeat_hits"]
    return out_rows


def to_dataframe(out_rows):
    return pd.DataFrame(out_rows, columns=FIELDS)


def write_table(out_rows, path):
    """Write the assembled rows to a TSV file. Sorting/filtering is left to
    the caller."""
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    df = to_dataframe(out_rows)
    df.to_csv(path, sep="\t", index=False, na_rep="")
    return path
