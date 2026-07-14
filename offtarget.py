"""Off-target I/O for Cas-OFFinder. The search itself runs externally
(network may be off here) -- this module only builds the query file and
parses/scores the result.

Run command (document, do not execute):
    cas-offinder <in> C <out>
(use the unmasked GRCh38 build; 'C' selects the CPU device; mm=4 by default
because repeats matter -- LTRs resemble HERVs.)
"""

import os

from cfd import calc_cfd


def emit_offinder(guides, genome_path, out_path, mm=4):
    """Write a Cas-OFFinder input file.

    guides: iterable of 20-nt spacer strings (5'->3', ACGT).
    genome_path: path to the (unmasked) reference genome directory, written
        as-is on line 1 -- Cas-OFFinder resolves it at run time.
    mm: mismatches allowed, written on each query line.

    File format:
        line 1: genome_path
        line 2: pattern, N*20 + "NGG"
        line 3+: "<spacer>NGG <mm>" per guide
    """
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(genome_path + "\n")
        f.write("N" * 20 + "NGG" + "\n")
        for g in guides:
            f.write(f"{g}NGG {mm}\n")
    return out_path


def score_offinder(out_file, guides, mm_scores, pam_scores, cfd_threshold=0.2,
                    repeat_hit_threshold=50):
    """Parse Cas-OFFinder output and summarize per guide.

    Cas-OFFinder output columns (tab or whitespace separated):
        query  chrom  position  target_seq  strand  mismatches

    target_seq is the 23-nt genomic hit (20-nt protospacer + 3-nt PAM).
    For each hit, calc_cfd(guide, hit_protospacer, hit_pam2nt) scores it;
    offtarget_count = number of hits scoring above cfd_threshold.
    n_repeat_hits flags guides with more than repeat_hit_threshold total
    genomic matches (regardless of score).
    """
    guides = list(guides)
    guide_of_query = {f"{g}NGG": g for g in guides}

    summary = {
        g: {"hits": 0, "offtarget_count": 0, "n_repeat_hits": False}
        for g in guides
    }

    n_lines = 0
    n_matched = 0
    with open(out_file) as f:
        for line in f:
            n_lines += 1
            if n_lines % 100000 == 0:
                print(f"  ... {n_lines} lines parsed, {n_matched} matched a query guide so far")
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            query, chrom, position, target_seq, strand, mismatches = fields[:6]

            guide = guide_of_query.get(query)
            if guide is None:
                continue
            n_matched += 1

            target_seq = target_seq.upper()
            hit_protospacer = target_seq[:20]
            hit_pam2nt = target_seq[21:23]

            summary[guide]["hits"] += 1

            if hit_pam2nt not in pam_scores or any(b not in "ACGT" for b in hit_protospacer):
                continue
            cfd = calc_cfd(guide, hit_protospacer, hit_pam2nt, mm_scores, pam_scores)
            if cfd > cfd_threshold:
                summary[guide]["offtarget_count"] += 1

    for g in guides:
        summary[g]["n_repeat_hits"] = summary[g]["hits"] > repeat_hit_threshold

    print(f"  parsed {n_lines} output line(s), {n_matched} matched one of our {len(guides)} guides")
    return summary
