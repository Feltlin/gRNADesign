#!/usr/bin/env python
"""HIV-1 LTR gRNA Designer -- CLI entry point.

Pipeline (see README.md section 2 for the numbered workflow and section 4
for the output column glossary):
    load alignment
      -> column frequencies (+ gap fraction)        [plus and RC alignment separately]
      -> consensus (per column, majority ACGT)
      -> per-column Shannon entropy                 [conservation track]
      -> enumerate 23-col windows, both strands
          gate: gap content low AND GG present      [drop failing windows here]
      -> per surviving window: S, E[CFD] indep/exact, GG fraction, plus-frame coordinate
      -> assemble table
      -> (optional, only if --genome or --offinder-results given) off-target
      -> (optional, only if --plots given) diagnostic plots

Off-target scoring is entirely optional: without a reference genome there is
nothing to search against, so Offtarget_Count / Offtarget_Repeat_Flag are
just left blank.

Usage:
    python gRNADesign.py in.fasta out.tsv
    python gRNADesign.py in.fasta out.tsv --genome /path/to/GRCh38_dir
    python gRNADesign.py in.fasta out.tsv --offinder-results cas_offinder_out.txt
    python gRNADesign.py in.fasta out.tsv --plots plots/
"""

import argparse
import shutil
import subprocess
import sys
import time

from config import CONFIG
from io_fasta import read_alignment, reverse_complement
from columns import column_frequencies, consensus, entropy
from cfd import load_cfd, cfd_position_weights
from windows import enumerate_windows
from table import assemble_table, write_table, merge_offtarget
from offtarget import emit_offinder, score_offinder



def run_pipeline(fasta_path, mm_path, pam_path, g_max, p_pam_min):
    t0 = time.time()

    print("Step 1/6: Loading alignment & building reverse complement")
    print(f"Reading {fasta_path}")
    rows_plus = read_alignment(fasta_path)
    L = len(rows_plus[0])
    print(f"  loaded {len(rows_plus)} sequences x {L} columns")

    print("Building reverse-complement alignment (for minus-strand search)")
    rows_minus = reverse_complement(rows_plus)
    print("  Done")



    print("Step 2/6: Loading CFD tables & deriving position weights")
    print(f"Reading {mm_path} & {pam_path} ...")
    mm_scores, pam_scores = load_cfd(mm_path, pam_path)
    print(f"  Mismatch table: {len(mm_scores)} keys; PAM table: {len(pam_scores)} keys")

    print("Computing per-position CFD weights (seed-heavy weighting vector) ...")
    w = cfd_position_weights(mm_scores, pam_scores)
    print("  Done")



    print("Step 3/6: Column statistics + window enumeration, both strands")
    all_windows = []
    for strand, rows in (("+", rows_plus), ("-", rows_minus)):
        print(f"[{strand} strand] computing base frequencies, consensus, entropy ...")
        freqs, gapfrac = column_frequencies(rows)
        cons = consensus(freqs)
        H = entropy(freqs)

        print(f"[{strand} strand] sliding the 23-col window and applying gap/PAM gates ...")
        wins = enumerate_windows(cons, freqs, gapfrac, rows, strand, g_max, p_pam_min, L=L)
        print(f"[{strand} strand] {len(wins)} windows survived both gates")

        all_windows.extend((win, freqs, rows, H) for win in wins)



    print("Step 4/6: Scoring surviving windows")
    print(f"Scoring {len(all_windows)} surviving window(s) "
          f"(Weighted_Entropy + CFD_Score_Independent/Exact for each) ...")
    out_rows = []
    for n, (win, freqs, rows, H) in enumerate(all_windows, start=1):
        out_rows.extend(assemble_table([win], H, w, rows, freqs, mm_scores, pam_scores))
        if n % 10 == 0 or n == len(all_windows):
            print(f"  scored {n}/{len(all_windows)}")

    print(f"\nTotal candidates: {len(out_rows)}  (pipeline took {time.time() - t0:.1f}s)")

    return out_rows


def run_offtarget(out_rows, args):
    """Optional stage: only runs if a genome to search,
    or an already-generated Cas-OFFinder result file to score provided.
    Mutates out_rows in place (fills Offtarget_Count / Offtarget_Repeat_Flag).
    """

    print("Step 5/6 (optional): Off-target scoring")
    guides = [r["_spacer"] for r in out_rows]
    print(f"{len(guides)} guide(s) to check; loading CFD tables ...")
    mm_scores, pam_scores = load_cfd(args.mm_pickle, args.pam_pickle)

    offinder_out_file = args.offinder_results

    if offinder_out_file is None and args.genome:
        query_path = args.offinder_query or (args.out + ".offinder_in.txt")
        print(f"Writing Cas-OFFinder query file (genome: {args.genome}, mm={args.offinder_mm}) ...")
        emit_offinder(guides, args.genome, query_path, mm=args.offinder_mm)
        print(f"  wrote {query_path}")

        cas_offinder_bin = shutil.which("cas-offinder")
        if cas_offinder_bin is None:
            print(
                "cas-offinder not found on PATH -- run it manually, then rerun "
                f"this tool with --offinder-results <output>:\n"
                f"  cas-offinder {query_path} C <cas_offinder_out>"
            )
            return
        offinder_out_file = query_path + ".out"
        print(f"Running Cas-OFFinder (this can take a while on a full genome): "
              f"{cas_offinder_bin} {query_path} C {offinder_out_file}")
        subprocess.run([cas_offinder_bin, query_path, "C", offinder_out_file], check=True)
        print("  Cas-OFFinder finished")

    if offinder_out_file is None:
        print("No genome and no existing results given -- nothing to score.")
        return

    print(f"Parsing and CFD-scoring hits from {offinder_out_file} ...")
    summary = score_offinder(
        offinder_out_file, guides, mm_scores, pam_scores,
        cfd_threshold=CONFIG["offtarget_cfd_threshold"],
        repeat_hit_threshold=CONFIG["repeat_hit_threshold"],
    )
    merge_offtarget(out_rows, summary)
    n_flagged = sum(1 for g in guides if summary[g]["n_repeat_hits"])
    print(f"  done: {sum(s['hits'] for s in summary.values())} total genomic hits parsed, "
          f"{n_flagged} guide(s) flagged as high-repeat")


def main(argv=None):
    parser = argparse.ArgumentParser(description="HIV-1 LTR gRNA Designer")
    parser.add_argument("fasta", nargs="?", default=CONFIG["fasta_path"],
                         help="Aligned FASTA panel (default: %(default)s)")
    parser.add_argument("out", nargs="?", default=CONFIG["table_path"],
                         help="Output TSV path (default: %(default)s)")
    parser.add_argument("--mm-pickle", default=CONFIG["mm_pickle_path"])
    parser.add_argument("--pam-pickle", default=CONFIG["pam_pickle_path"])
    parser.add_argument("--g-max", type=float, default=CONFIG["g_max"])
    parser.add_argument("--p-pam-min", type=float, default=CONFIG["p_pam_min"])

    offtarget_group = parser.add_argument_group(
        "off-target (optional -- skipped entirely if neither is given)")
    offtarget_group.add_argument(
        "--genome", default=CONFIG["genome_path"],
        help="Path to an unmasked reference genome (e.g. GRCh38) to search against")
    offtarget_group.add_argument(
        "--offinder-query", default=None,
        help="Where to write the Cas-OFFinder query file (default: <out>.offinder_in.txt)")
    offtarget_group.add_argument(
        "--offinder-results", default=None,
        help="An already-generated Cas-OFFinder output file to score, skipping the search")
    offtarget_group.add_argument("--offinder-mm", type=int, default=CONFIG["offinder_mm"])

    parser.add_argument("--plots", metavar="OUT_DIR", nargs="?",
                         const=CONFIG["plots_dir"], default=None,
                         help="Also write the 3 diagnostic plots to OUT_DIR "
                              "(defaults to %(const)s if OUT_DIR is omitted)")
    args = parser.parse_args(argv)

    out_rows = run_pipeline(args.fasta, args.mm_pickle, args.pam_pickle,
                             args.g_max, args.p_pam_min)

    if args.genome or args.offinder_results:
        run_offtarget(out_rows, args)
    else:
        print("Step 5/6 (optional): Off-target scoring")
        print("No --genome or --offinder-results given; skipping.")

    print("Step 6/6: Writing output table")
    write_table(out_rows, args.out)
    print(f"Wrote {len(out_rows)} rows to {args.out}")

    if args.plots:
        print("Bonus step (optional): Generating diagnostic plots")
        from plotting import generate_all_plots
        generate_all_plots(args.plots, args.fasta, args.mm_pickle, args.pam_pickle, out_rows)

    print("Done")

    return 0


if __name__ == "__main__":
    sys.exit(main())
