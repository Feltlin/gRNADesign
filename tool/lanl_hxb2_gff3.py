#!/usr/bin/env python
"""Convert the LANL HXB2 landmark spreadsheet into GFF3.

This is an optional helper for the HIV-1 worked example, not part of the
pipeline. It exists because the LANL sheet is the authoritative source for
HIV-1 LTR regulatory element, and it is distributed as a base-by-base
spreadsheet rather than as an annotation file.

    curl -L https://www.hiv.lanl.gov/content/sequence/HIV/MAP/HXB2.xlsx \\
        -o input/HXB2.xlsx
    python tool/lanl_hxb2_gff3.py input/HXB2.xlsx input/HXB2.gff3

The sheet has one row per base of HXB2 (GenBank K03455, 9719 nt) and three
annotation source, each read here:

* "RNA/DNA features"       a label repeated on every base of a feature
* "RNA/DNA features notes" an explicit "<name> start" / "<name> end" marker
* the three "protein" column, one per reading frame

Nothing is invented: every coordinate is a row number of the sheet.

Needs openpyxl (pip install openpyxl).
"""

import sys
from collections import OrderedDict

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import openpyxl

#: Sequence Ontology type for a landmark whose name matches, first hit win.
#: Anything unmatched becomes a plain region, which is still drawn and still
#: reported, just not specifically typed.
KIND = (
    ("TATA", "TATA_box"),
    ("NF-kappa-B", "protein_binding_site"),
    ("Sp1", "protein_binding_site"),
    ("TCF-1", "protein_binding_site"),
    ("binding site", "protein_binding_site"),
    ("binding region", "protein_binding_site"),
    ("primer binding", "primer_binding_site"),
    ("Poly-A signal", "polyA_signal_sequence"),
    ("stem-loop", "stem_loop"),
    ("TAR", "stem_loop"),
    ("secondary structure", "stem_loop"),
    ("Response Element", "regulatory_region"),
    ("Responsive Element", "regulatory_region"),
)

#: Shorter display name. The sheet is verbose; a plot label is not.
RENAME = {
    "NF-kappa-B-I site": "NF-κB I",
    "NF-kappa-B-II site": "NF-κB II",
    "Sp1-I site": "Sp1 I",
    "Sp1-II site": "Sp1 II",
    "Sp1-III site": "Sp1 III",
    "TCF-1alpha binding site": "TCF-1α",
    "TATA Box": "TATA box",
    "Poly-A signal sequence": "poly(A) signal",
    "Lys tRNA primer binding site": "PBS",
    "TAR element": "TAR",
    "5' LTR U3": "U3 (5')",
    "5' LTR R repeat": "R (5')",
    "5' LTR U5": "U5 (5')",
    "3' LTR U3": "U3 (3')",
    "3' LTR R repeat": "R (3')",
    "3' LTR U5": "U5 (3')",
    "Rev Responsive Element (RRE) region": "RRE",
    "+1 mRNA start site": "+1 TSS",
    "Extensive secondary structure": "U5 structure",
    "Interferon-Stimulated Response Element": "ISRE",
    "Packaging signal RNA secondary structure": "packaging signal",
    "Central polypurine tract": "cPPT",
}


def kind_of(name):
    for needle, kind in KIND:
        if needle.lower() in name.lower():
            return kind
    return "region"


def runs(column, gap=0):
    """Turn a per base label column into (start, end, label), 1 based.

    `gap` tolerates that many unlabelled base inside a run. The protein
    column label one base per codon, so without a tolerance of 2 a single
    gene would come back as hundreds of separate three base fragment.
    """
    out, current, begin, last = [], None, None, None
    for position, label in enumerate(column, start=1):
        if not label:
            continue
        if label == current and position - last <= gap + 1:
            last = position
            continue
        if current:
            out.append((begin, last, current))
        current, begin, last = label, position, position
    if current:
        out.append((begin, last, current))
    return out


def marker(note):
    """Pair every "<name> start" with its "<name> end" down the sheet.

    A note cell may hold several marker separated by "; ". A marker with
    neither suffix is a single position landmark.
    """
    open_at, out = OrderedDict(), []
    for position, text in enumerate(note, start=1):
        if not text:
            continue
        for part in str(text).split(";"):
            part = part.strip()
            if not part:
                continue
            lower = part.lower()
            if lower.endswith(" start"):
                open_at.setdefault(lower[:-6], []).append((position, part[:-6]))
            elif lower.endswith(" end"):
                name = lower[:-4]
                if open_at.get(name):
                    begin, label = open_at[name].pop(0)
                    out.append((begin, position, label))
                else:
                    print(f"  warning: '{part}' at {position} has no start")
            else:
                out.append((position, position, part))
    for name, pending in open_at.items():
        for position, _ in pending:
            print(f"  warning: '{name} start' at {position} is never closed")
    return out


def main(source, destination, sheet_name="HXB2"):
    book = openpyxl.load_workbook(source, read_only=True, data_only=True)
    row = list(book[sheet_name].iter_rows(values_only=True))
    header = ["" if value is None else str(value).strip() for value in row[1]]
    body = row[2:]

    def column(index):
        return [
            "" if index >= len(entry) or entry[index] is None
            else str(entry[index]).strip().replace("-", "").strip()
            if str(entry[index]).strip() == "-" else str(entry[index]).strip()
            for entry in body
        ]

    position = [entry[0] for entry in body if entry[0] is not None]
    length = max(int(value) for value in position)

    feature = []
    # 1. the broad RNA/DNA feature label
    index = header.index("RNA/DNA features")
    feature += [(s, e, n, kind_of(n), "RNA/DNA features column")
                for s, e, n in runs(column(index))]
    # 2. the explicit start/end landmark
    index = header.index("RNA/DNA features notes")
    feature += [(s, e, n, kind_of(n), "RNA/DNA features notes column")
                for s, e, n in marker(column(index))]
    # 3. one protein column per reading frame
    for index, name in enumerate(header):
        if name == "protein":
            frame = len([x for x in header[:index] if x == "protein"]) + 1
            feature += [(s, e, n, "CDS",
                         f"protein column, reading frame {frame}")
                        for s, e, n in runs(column(index), gap=2)]

    feature = [entry for entry in feature if entry[2]]
    feature.sort(key=lambda entry: (entry[0], entry[1]))

    out = [
        "##gff-version 3",
        f"##sequence-region K03455.1 1 {length}",
        "# HIV-1 HXB2 annotation, converted from the LANL landmark spreadsheet",
        "# https://www.hiv.lanl.gov/content/sequence/HIV/MAP/HXB2.xlsx",
        "# by tool/lanl_hxb2_gff3.py. Coordinate are the sheet's own row number,",
        "# which is HXB2 / GenBank K03455 numbering, 1 based inclusive.",
    ]
    for start, end, name, kind, origin in feature:
        label = RENAME.get(name, name)
        out.append("\t".join([
            "K03455.1", "LANL", kind, str(start), str(end), ".", "+", ".",
            f"Name={label};source_column={origin}",
        ]))
    with open(destination, "w", encoding="utf-8") as handle:
        handle.write("\n".join(out) + "\n")

    print(f"wrote {len(feature)} feature to {destination}")
    for start, end, name, kind, _ in feature:
        if end <= 700:
            print(f"  {start:>5}-{end:<5} {kind:22s} {RENAME.get(name, name)}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
