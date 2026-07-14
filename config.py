"""CONFIG block: paths, enzyme, gates, thresholds. Edit here, not in callers."""

import os

_HERE = os.path.dirname(os.path.abspath(__file__))

CONFIG = {
    # --- inputs ---
    "fasta_path": os.path.join(_HERE, "input", "hiv-db.fasta"),
    "mm_pickle_path": os.path.join(_HERE, "data", "mismatch_score.pkl"),
    "pam_pickle_path": os.path.join(_HERE, "data", "pam_scores.pkl"),

    # --- nuclease (SpCas9) ---
    "pam": "NGG",
    "spacer_len": 20,
    "pam_len": 3,

    # --- window gates ---
    "g_max": 0.10,       # max gap fraction allowed in any of the 23 window columns
    "p_pam_min": 0.80,   # min fraction of variants with GG at the PAM positions

    # --- off-target (Cas-OFFinder) ---
    "genome_path": None,     # path to unmasked GRCh38, set before running emit_offinder
    "offinder_mm": 4,        # allowed mismatches for cas-offinder search
    "offtarget_cfd_threshold": 0.2,
    "repeat_hit_threshold": 50,

    # --- output ---
    "table_path": os.path.join(_HERE, "output", "gRNA_table.tsv"),
    "plots_dir": os.path.join(_HERE, "output", "plots"),
}
