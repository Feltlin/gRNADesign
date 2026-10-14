"""Build a small synthetic panel and reference, so the end to end test owns
its own data and runs on a fresh clone.

The repository ships no sequence: `input/` is untracked because a real panel
runs to hundreds of MB. A test that depended on a real file would simply be
skipped most of the time, which is the same as not having it.
"""

import random

BASE = "ACGT"


def reference_sequence(length=400, seed=11):
    """A reference with a usable NGG density, planted deterministically."""
    rng = random.Random(seed)
    base = [rng.choice(BASE) for _ in range(length)]
    # Plant a PAM every 40 nt so there is something for the gate to find.
    for position in range(30, length - 3, 40):
        base[position + 1] = "G"
        base[position + 2] = "G"
    return "".join(base)


def panel(reference, depth=24, divergence=0.02, gap_block=(120, 126), seed=7):
    """A panel of `depth` row derived from the reference.

    Every row is the reference with a few substitution, plus one block of
    column that is mostly gap in most row, which is what the gap gate and the
    core column rule exist to handle.
    """
    rng = random.Random(seed)
    row = []
    for index in range(depth):
        base = list(reference)
        for position in range(len(base)):
            if rng.random() < divergence:
                base[position] = rng.choice(BASE)
        if index % 4:  # most row carry the gap, a few do not
            for position in range(*gap_block):
                base[position] = "-"
        row.append("".join(base))
    return row


def write(directory, reference, row, name="panel.fasta", reference_name="ref.fasta"):
    """Write both as FASTA and return (panel path, reference path)."""
    panel_path = directory / name
    reference_path = directory / reference_name
    panel_path.write_text(
        "".join(f">seq{index}\n{entry}\n" for index, entry in enumerate(row)),
        encoding="utf-8")
    reference_path.write_text(f">testref\n{reference}\n", encoding="utf-8")
    return panel_path, reference_path


def annotation_gff3(length=400):
    """A GFF3 over the synthetic reference, with a gene and a nested element."""
    return (
        "##gff-version 3\n"
        f"##sequence-region testref 1 {length}\n"
        "testref\ttest\tgene\t40\t200\t.\t+\t.\tName=geneA\n"
        "testref\ttest\tregion\t1\t400\t.\t+\t.\tName=whole locus\n"
        "testref\ttest\tTATA_box\t60\t64\t.\t+\t.\tName=TATA box\n"
    )
