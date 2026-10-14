"""Every tunable setting in one place. Edit here, not in a caller.

Each key is also exposed as a command line flag, so this file only sets the
default.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent

CONFIG = {
    # --- Input ---
    # No default. The alignment is a required argument and the reference is
    # given with --reference, because nothing about this tool is tied to one
    # organism. README.md, "Targeting The HIV-1 Provirus", is a worked example
    # with the command that fetch a real data set.
    "alignment_path": None,
    "reference_path": None,

    # The only data this repository ships: the published CFD table.
    "mismatch_table_path": ROOT / "data" / "mismatch_score.pkl",
    "pam_table_path": ROOT / "data" / "pam_scores.pkl",

    # --- Nuclease ---
    # PAM is an IUPAC motif. "3prime" places it immediately 3' of the spacer
    # (SpCas9 NGG); "5prime" places it immediately 5' of the spacer (Cas12a
    # TTTV). CFD scoring is SpCas9 derived and needs a 20 nt spacer.
    "pam": "NGG",
    "pam_side": "3prime",
    "spacer_length": 20,

    # --- Gate ---
    # A window survives only if every column is below max_gap_fraction, the
    # consensus PAM matches the motif, and at least min_pam_fraction of the
    # aligned sequence carries a motif matching PAM at that locus.
    "max_gap_fraction": 0.10,
    "min_pam_fraction": 0.80,

    # --- Reference anchoring ---
    # Restrict the reference search to a 1 based inclusive (start, end), or
    # None to search the whole reference. Needed when the panel matches a
    # repeated element, because a local alignment cannot choose between two
    # identical copy on its own. See README.md, "Anchoring The Panel".
    "reference_region": None,
    # A column is used to anchor the panel to the reference only when its gap
    # fraction is at or below this value.
    "core_gap_fraction": 0.50,

    # A separate GFF3, GTF or BED annotation, used when the reference sequence
    # is a plain FASTA or when its own annotation should be replaced. None
    # keeps whatever the reference file carries.
    "annotation_path": None,

    # Display name substituted for a recorded feature name, for a record whose
    # feature name are archaic or terse. The coordinate is never touched.
    # Example: {"sor 23K protein": "vif"} for the 1985 HIV-1 GenBank record.
    "feature_alias": {},

    # --- Minimum free energy ---
    # "auto" uses a local ViennaRNA install when importable and the RNAfold
    # web server otherwise. The server takes one sequence per job, so a run
    # costs one job per distinct guide; the cache makes a repeat run free.
    "mfe_backend": "auto",
    "mfe_url": "http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi",
    "mfe_temperature": 37.0,
    "mfe_worker": 4,
    "mfe_cache_path": ROOT / "data" / "mfe_cache.json",
    # Appended to the spacer before folding. Empty folds the spacer alone;
    # set the constant region of the sgRNA to fold the whole guide.
    "scaffold": "",

    # --- Plot ---
    "entropy_window": 20,
    # An optional horizontal guide on the CFD panel. None draws none: the y
    # axis already carry the value, so a threshold line is decoration unless
    # a particular cut off is being argued for.
    "cfd_reference_line": None,

    # --- Output ---
    "table_path": ROOT / "output" / "gRNA_table.tsv",
    "plot_path": ROOT / "output" / "gRNA_map.png",
    "consensus_path": ROOT / "output" / "consensus.fasta",
}
