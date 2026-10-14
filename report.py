"""The output table.

This module owns the reader facing column name. Internal code keeps its own
short field name; only this translation layer knows both.
"""

from pathlib import Path

import pandas as pd

COLUMN = (
    "gRNA5to3",
    "Derived_Strand",
    "Entropy_Mean",
    "PAM_Conservation",
    "CFD_Score",
    "MFE",
    "Fold_Structure",
    "Gene",
    "Element",
    "Ref_Start",
    "Ref_End",
    "Align_Start",
    "Align_End",
)

ROUND = {
    "Entropy_Mean": 4,
    "PAM_Conservation": 4,
    "CFD_Score": 4,
    "MFE": 2,
}


def build_row(candidate, entropy_mean, cfd_score, reference_span, gene, element):
    """One table row.

    gRNA5to3 is a display string, "SPACER(PAM)", read 5' to 3' in the guide's
    own frame. The PAM is bracketed because it is **not** part of the guide
    you order: it is the motif that has to be present in the target, next to
    the protospacer. Order the 20 nt outside the bracket.

    The row also carries a bare "_spacer" used internally for folding and
    look up. It is not in COLUMN, so it never reaches the file.

    Derived_Strand is a "read off" label, never a "binds to" label: '+' means
    the protospacer and its PAM are read off the plus strand of the reference,
    so the guide RNA itself pairs with the minus strand.
    """
    start, end = reference_span
    return {
        "gRNA5to3": f"{candidate.spacer}({candidate.pam})",
        "_spacer": candidate.spacer,
        "Derived_Strand": candidate.strand,
        "Entropy_Mean": entropy_mean,
        "PAM_Conservation": candidate.pam_conservation,
        "CFD_Score": cfd_score,
        "MFE": None,
        "Fold_Structure": None,
        "Gene": gene,
        "Element": element,
        "Ref_Start": start,
        "Ref_End": end,
        "Align_Start": candidate.plus_start,
        "Align_End": candidate.plus_end,
    }


def merge_fold(row, folded):
    """Fill MFE and Fold_Structure from {spacer: (structure, mfe)}."""
    for entry in row:
        structure, energy = folded.get(entry["_spacer"], (None, None))
        entry["Fold_Structure"] = structure
        entry["MFE"] = energy
    return row


def to_frame(row):
    frame = pd.DataFrame(row, columns=COLUMN)
    for column, digit in ROUND.items():
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce").round(digit)
    for column in ("Ref_Start", "Ref_End", "Align_Start", "Align_End"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("Int64")
    return frame.sort_values(
        ["Ref_Start", "Align_Start", "Derived_Strand"], na_position="last"
    ).reset_index(drop=True)


def write_table(row, path):
    """Write the assembled row as a TSV, sorted along the reference."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    to_frame(row).to_csv(path, sep="\t", index=False, na_rep="")
    return path
