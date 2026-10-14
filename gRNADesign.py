#!/usr/bin/env python
"""gRNA Designer - command line entry point.

Design a guide RNA panel from an aligned sequence set, a PAM motif and an
annotated reference. See README.md, "Data Process", for the step below in
detail.

    read the alignment
      -> per column frequency, consensus and Shannon entropy
      -> anchor the panel to the annotated reference
      -> slide a (spacer + PAM) window over both strand
           gate: gap content low, consensus PAM matches, PAM present in the panel
      -> per surviving window: entropy mean, PAM conservation, mean CFD
      -> (optional) fold each guide for its MFE and structure
      -> write the table
      -> (optional) draw the candidate map
      -> (optional) write the consensus FASTA

Usage:
    python gRNADesign.py
    python gRNADesign.py input/panel.fasta output/gRNA_table.tsv
    python gRNADesign.py --pam NGG --reference input/K03455.gb --plot
    python gRNADesign.py --mfe --plot --consensus
"""

import argparse
import sys
import time

import numpy as np

import alignment as aln
import annotation
import cfd
import mfe as mfe_module
import plot as plot_module
import report
import score as score_module
from config import CONFIG
from pam import enumerate_candidate


def design(
    alignment_path,
    reference_path,
    *,
    mismatch_table_path,
    pam_table_path,
    pam,
    pam_side,
    spacer_length,
    max_gap_fraction,
    min_pam_fraction,
    reference_region,
    core_gap_fraction,
    annotation_path=None,
    feature_alias=None,
    log=print,
):
    """Run the pipeline and return (row, context) where context carries what
    the plot and the run log need."""
    started = time.time()

    log("Step 1/5  Read the alignment")
    panel = aln.read_alignment(alignment_path)
    log(f"  {alignment_path}")
    log(f"  {panel.depth} sequence x {panel.length} column")

    log("Step 2/5  Per column frequency, consensus and Shannon entropy")
    frequency, gap_fraction = aln.column_frequency(panel.matrix)
    consensus = aln.consensus(frequency)
    entropy = aln.shannon_entropy(frequency)
    log(f"  mean H {entropy.mean():.4f}, mean gap fraction {gap_fraction.mean():.4f}")

    # Equal row length is necessary for an alignment but not sufficient: a set
    # of unaligned sequence all trimmed to the same length passes that check
    # and then scores nonsense, because the column do not correspond. A real
    # alignment either carries gap or is highly conserved; carrying neither is
    # the signature of input that was never aligned.
    if entropy.mean() > 0.30 and gap_fraction.mean() < 0.001:
        log(f"  WARNING: mean column entropy is {entropy.mean():.2f} and there "
            f"is essentially no gap. That combination usually means the input "
            f"is not actually aligned, only equal in length. Align the panel "
            f"first (MAFFT, MUSCLE, or whatever your field uses) and rerun. "
            f"If the panel really is this diverse, ignore this.")

    log("Step 3/5  Anchor the panel to the annotated reference")
    reference = None
    mapping = np.full(panel.length, -1, dtype=np.int64)
    note = "not anchored, no reference given"

    if reference_path is None:
        log("  no --reference given: the position and annotation column are "
            "left empty. Every score is unaffected, since all of them are "
            "computed from the panel alone")
    else:
        reference = annotation.load_reference(
            reference_path, reference_region, annotation_path, feature_alias)
        log(f"  {reference_path}"
            + (f" + {annotation_path}" if annotation_path else ""))
        log(f"  {reference.identifier}  {reference.full_length} nt, "
            f"{len(reference.feature)} annotated span")

        orientation, forward_score, reverse_score = annotation.detect_orientation(
            consensus, gap_fraction, reference, core_gap_fraction)
        if orientation == "reverse":
            log(f"  the panel is the reverse complement of the reference "
                f"(probe score {reverse_score:.0f} against {forward_score:.0f}); "
                f"flipping it so every reported coordinate and strand is in "
                f"reference orientation")
            panel = aln.Alignment(panel.identifier,
                                  aln.reverse_complement(panel.matrix))
            frequency, gap_fraction = aln.column_frequency(panel.matrix)
            consensus = aln.consensus(frequency)
            entropy = aln.shannon_entropy(frequency)

        mapping, note = annotation.map_column(
            consensus, gap_fraction, reference, core_gap_fraction)
        log(f"  anchored to {note}")

        anchored = int((mapping >= 0).sum())
        if anchored < 0.5 * panel.length:
            log(f"  WARNING: only {anchored} of {panel.length} column anchored. "
                f"The panel and the reference may not be the same locus, or the "
                f"panel may cover only part of it. Check the span above before "
                f"trusting any coordinate.")

    log(f"Step 4/5  Slide the {spacer_length} + {len(pam)} nt window over both strand")
    mismatch_score, pam_score = cfd.load_table(mismatch_table_path, pam_table_path)
    table = cfd.activity_table(mismatch_score, pam_score, spacer_length)
    if table is None:
        log(f"  note: CFD is defined for a 20 nt spacer only, so the CFD column "
            f"is left empty for this {spacer_length} nt spacer")

    candidate = []
    for strand, matrix in (("+", panel.matrix),
                           ("-", aln.reverse_complement(panel.matrix))):
        if strand == "+":
            strand_frequency, strand_gap = frequency, gap_fraction
        else:
            strand_frequency, strand_gap = aln.column_frequency(matrix)
        found = enumerate_candidate(
            matrix, strand_gap, aln.consensus(strand_frequency), pam, strand,
            spacer_length=spacer_length,
            pam_side=pam_side,
            max_gap_fraction=max_gap_fraction,
            min_pam_fraction=min_pam_fraction,
        )
        log(f"  [{strand} strand] {len(found)} window survived every gate")
        candidate.extend((entry, matrix) for entry in found)

    log("Step 5/5  Score every surviving window")
    row = []
    for entry, matrix in candidate:
        reference_span = annotation.span(mapping, entry.plus_spacer_column)
        gene, element = (annotation.annotate(reference, *reference_span)
                         if reference is not None else ("", ""))
        row.append(report.build_row(
            entry,
            entropy_mean=score_module.entropy_mean(entropy, entry.plus_spacer_column),
            cfd_score=score_module.cfd_mean(matrix, entry, pam, table, pam_score),
            reference_span=reference_span,
            gene=gene,
            element=element,
        ))
    anchored = sum(1 for entry in row if entry["Ref_Start"] is not None)
    log(f"  {len(row)} candidate, {anchored} with a reference coordinate "
        f"({time.time() - started:.1f}s)")

    context = {
        "panel": panel,
        "consensus": consensus,
        "entropy": entropy,
        "gap_fraction": gap_fraction,
        "mapping": mapping,
        "reference": reference,
        "anchor_note": note,
    }
    return row, context


def add_fold(row, argument, log=print):
    """Fold every distinct guide and merge its MFE and structure in."""
    log("Bonus  Fold each guide for its minimum free energy")
    folded = mfe_module.fold(
        [entry["_spacer"] for entry in row],
        backend=argument.mfe_backend,
        url=argument.mfe_url,
        temperature=argument.mfe_temperature,
        worker=argument.mfe_worker,
        cache_path=argument.mfe_cache,
        scaffold=argument.scaffold,
        log=log,
    )
    return report.merge_fold(row, folded)


def draw(row, context, argument, log=print):
    """Draw the candidate map."""
    log("Bonus  Draw the candidate map")
    frame = report.to_frame(row)
    mapping = context["mapping"]
    reference = context["reference"]

    running = aln.running_mean(context["entropy"], argument.entropy_window)
    anchored = mapping >= 0
    position = mapping[anchored] + 1
    value = running[anchored]
    order = np.argsort(position)

    side = "3" if argument.pam_side == "3prime" else "5"
    best = frame["CFD_Score"].max()
    descriptor = reference.description.split(";")[0].strip()
    span = f"{reference.offset + 1}-{reference.offset + len(reference.sequence)}"
    subtitle = (
        f"{descriptor}  ·  {reference.identifier} {span}\n"
        f"{context['panel'].depth} aligned sequence  ·  PAM {argument.pam} "
        f"{side}' of a {argument.spacer_length} nt spacer  ·  "
        f"{len(frame)} candidate past the gate, best CFD {best:.2f}"
    )
    footnote = (
        "Entropy and CFD are independent: entropy is how conserved the locus is "
        "across the panel, CFD is how well the one consensus guide still cuts "
        "every variant at that locus.\n"
        "Arrow points 5' to 3' along the strand the protospacer is read from. "
        "Score come from the panel alignment only, not from a genome wide off "
        "target search; screen a shortlist against the host genome before "
        "ordering.\n"
        f"Anchored to {context['anchor_note']}."
    )
    path = plot_module.draw(
        frame, position[order], value[order], reference, argument.plot,
        window=argument.entropy_window,
        cfd_line=argument.cfd_line,
        title="gRNA candidate map",
        subtitle=subtitle,
        footnote=footnote,
    )
    log(f"  wrote {path}")
    return path


def parse(argv=None):
    parser = argparse.ArgumentParser(
        description="Design a guide RNA panel from an aligned sequence set, "
                    "a PAM motif and an annotated reference.")
    parser.add_argument("alignment",
                        help="aligned FASTA panel, every row the same length")
    parser.add_argument("out", nargs="?", default=str(CONFIG["table_path"]),
                        help="output TSV path (default: %(default)s)")
    parser.add_argument("--reference", default=CONFIG["reference_path"],
                        help="annotated reference, GenBank or FASTA. Without "
                             "it the position and annotation column are left "
                             "empty and every score is unchanged")
    parser.add_argument("--annotation", default=CONFIG["annotation_path"],
                        help="separate GFF3, GTF or BED annotation; use it when "
                             "the reference is a plain FASTA")
    parser.add_argument("--reference-region", nargs=2, type=int, metavar=("START", "END"),
                        default=CONFIG["reference_region"],
                        help="RESTRICT anchoring to this span of the reference, "
                             "1 based with both end included. Omit it to search "
                             "the whole reference, which is the usual case. Use "
                             "it only when the panel matches a repeated element "
                             "and you must say which copy")

    nuclease = parser.add_argument_group("nuclease")
    nuclease.add_argument("--pam", default=CONFIG["pam"],
                          help="IUPAC PAM motif (default: %(default)s)")
    nuclease.add_argument("--pam-side", choices=("3prime", "5prime"),
                          default=CONFIG["pam_side"],
                          help="side of the spacer the PAM sits on "
                               "(default: %(default)s)")
    nuclease.add_argument("--spacer-length", type=int, default=CONFIG["spacer_length"])

    gate = parser.add_argument_group("gate")
    gate.add_argument("--max-gap-fraction", type=float,
                      default=CONFIG["max_gap_fraction"])
    gate.add_argument("--min-pam-fraction", type=float,
                      default=CONFIG["min_pam_fraction"])
    gate.add_argument("--core-gap-fraction", type=float,
                      default=CONFIG["core_gap_fraction"])

    folding = parser.add_argument_group(
        "minimum free energy (off unless --mfe is given)")
    folding.add_argument("--mfe", action="store_true",
                         help="fold every guide; uses a local ViennaRNA install "
                              "when available, otherwise the public RNAfold server")
    folding.add_argument("--mfe-backend", choices=("auto", "local", "server"),
                         default=CONFIG["mfe_backend"])
    folding.add_argument("--mfe-url", default=CONFIG["mfe_url"])
    folding.add_argument("--mfe-temperature", type=float,
                         default=CONFIG["mfe_temperature"])
    folding.add_argument("--mfe-worker", type=int, default=CONFIG["mfe_worker"])
    folding.add_argument("--mfe-cache", default=str(CONFIG["mfe_cache_path"]))
    folding.add_argument("--scaffold", default=CONFIG["scaffold"],
                         help="appended to the spacer before folding; empty "
                              "folds the spacer alone")

    output = parser.add_argument_group("extra output")
    output.add_argument("--plot", metavar="PATH", nargs="?",
                        const=str(CONFIG["plot_path"]), default=None,
                        help="draw the candidate map (default: %(const)s)")
    output.add_argument("--entropy-window", type=int, default=CONFIG["entropy_window"],
                        help="running mean window for the entropy line "
                             "(default: %(default)s)")
    output.add_argument("--cfd-line", type=float, default=CONFIG["cfd_reference_line"],
                        help="CFD reference line drawn on the map "
                             "(default: %(default)s)")
    output.add_argument("--consensus", metavar="PATH", nargs="?",
                        const=str(CONFIG["consensus_path"]), default=None,
                        help="write the panel consensus FASTA (default: %(const)s)")
    return parser.parse_args(argv)


def main(argv=None):
    argument = parse(argv)

    row, context = design(
        argument.alignment,
        argument.reference,
        mismatch_table_path=CONFIG["mismatch_table_path"],
        pam_table_path=CONFIG["pam_table_path"],
        pam=argument.pam,
        pam_side=argument.pam_side,
        spacer_length=argument.spacer_length,
        max_gap_fraction=argument.max_gap_fraction,
        min_pam_fraction=argument.min_pam_fraction,
        reference_region=argument.reference_region,
        core_gap_fraction=argument.core_gap_fraction,
        annotation_path=argument.annotation,
        feature_alias=CONFIG["feature_alias"],
    )

    if not row:
        print("No window survived the gate; nothing to write.")
        return 1

    if argument.mfe:
        add_fold(row, argument)
    else:
        print("Bonus  Minimum free energy: --mfe not given, column left empty")

    path = report.write_table(row, argument.out)
    print(f"Wrote {len(row)} row to {path}")

    if argument.plot and context["reference"] is None:
        print("Bonus  Candidate map: skipped, because the map is drawn against "
              "reference coordinate and no --reference was given")
    elif argument.plot:
        draw(row, context, argument)

    if argument.consensus:
        written = aln.write_fasta(
            context["consensus"], argument.consensus,
            header=f"consensus {context['panel'].depth} sequence")
        print(f"Wrote the {len(context['consensus'])} nt consensus to {written}")

    print("Done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
