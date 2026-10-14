"""The annotated reference, and the map from alignment column to reference position.

The alignment is gapped and the reference is not, so a column carries no
coordinate by itself. The panel is anchored by locally aligning its consensus
to the reference once, and every candidate coordinate is then read off that
one map.
"""

import urllib.parse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

_GENBANK_SUFFIX = {".gb", ".gbk", ".genbank", ".gbff"}
_GFF_SUFFIX = {".gff", ".gff3", ".gtf"}
_BED_SUFFIX = {".bed"}

#: Attribute key searched, in order, for the display name of a GFF feature.
_GFF_NAME_KEY = ("Name", "gene", "gene_name", "product", "locus_tag", "ID", "note")

#: Feature type naming a gene product. Everything else is treated as context
#: or as a regulatory element.
_GENE_KIND = ("CDS", "gene", "mat_peptide")

#: Feature type that is structural bookkeeping rather than something a guide
#: can usefully be said to hit.
_SKIP_KIND = ("source", "mRNA", "prim_transcript", "exon", "intron")


@dataclass(frozen=True)
class Feature:
    """One annotated span, 0 based half open, in full reference coordinate."""

    start: int
    end: int
    strand: int
    name: str
    kind: str
    seqid: str = ""


def _match_seqid(feature, identifier):
    """Keep the feature belonging to the reference record.

    A single sequence annotation often name its sequence differently from the
    FASTA header, so a total mismatch is read as "this file is about this one
    sequence" rather than as an empty annotation.
    """
    bare = identifier.split(".")[0].lower()
    kept = [entry for entry in feature
            if entry.seqid.split(".")[0].lower() in (bare, "")]
    return kept or feature


@dataclass(frozen=True)
class Reference:
    """The annotated reference, optionally restricted to a search region."""

    identifier: str
    description: str
    sequence: str          # the searched region only
    offset: int            # 0 based start of that region in the full reference
    full_length: int
    feature: tuple         # over the full reference


def _feature_name(feature):
    for key in ("gene", "product", "note", "standard_name", "label"):
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return feature.type


def _from_genbank(record):
    feature = []
    for entry in record.features:
        if entry.type == "source":
            continue
        name = _feature_name(entry)
        for part in entry.location.parts:
            feature.append(Feature(
                start=int(part.start),
                end=int(part.end),
                strand=int(part.strand or 1),
                name=name,
                kind=entry.type,
            ))
    return feature


def _parse_attribute(field):
    """Read a GFF3 (key=value;) or GTF (key "value";) attribute column."""
    attribute = {}
    for chunk in field.strip().strip(";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" in chunk:
            key, _, value = chunk.partition("=")
        else:
            key, _, value = chunk.partition(" ")
        attribute[key.strip()] = urllib.parse.unquote(value.strip().strip('"'))
    return attribute


def _from_gff(path):
    """Read GFF3 or GTF. Coordinate are 1 based inclusive on the file side."""
    feature = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            if line.startswith("##FASTA"):
                break
            continue
        field = line.split("\t")
        if len(field) < 8:
            continue
        attribute = _parse_attribute(field[8]) if len(field) > 8 else {}
        name = next((attribute[key] for key in _GFF_NAME_KEY if key in attribute),
                    field[2])
        feature.append(Feature(
            start=int(field[3]) - 1,
            end=int(field[4]),
            strand=-1 if field[6] == "-" else 1,
            name=name,
            kind=field[2],
            seqid=field[0],
        ))
    return feature


def _from_bed(path):
    """Read BED. Coordinate are already 0 based half open."""
    feature = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        field = line.split("\t") if "\t" in line else line.split()
        if len(field) < 3:
            continue
        feature.append(Feature(
            start=int(field[1]),
            end=int(field[2]),
            strand=-1 if len(field) > 5 and field[5] == "-" else 1,
            name=field[3] if len(field) > 3 else "feature",
            kind="gene",
            seqid=field[0],
        ))
    return feature


def load_reference(path, region=None, annotation_path=None, alias=None):
    """Read a reference sequence and its annotation.

    path: the sequence, as GenBank (annotation included) or FASTA.
    annotation_path: an optional separate GFF3, GTF or BED, which replaces
        whatever annotation `path` carried. This is the generic route: any
        reference shipping a standard annotation file works without the
        sequence having to be GenBank.
    region: a 1 based inclusive (start, end) restricting the anchoring search,
        which is how a panel matching a repeated element such as an LTR is
        pinned to one copy. Reported coordinate stays in full reference space.
    alias: {recorded name: display name}, applied to every feature name. It
        lets an archaic or terse record read in current gene name without the
        coordinate being touched.
    """
    path = Path(path)
    fmt = "genbank" if path.suffix.lower() in _GENBANK_SUFFIX else "fasta"
    with open(path, encoding="utf-8", errors="replace") as handle:
        record = next(SeqIO.parse(handle, fmt))

    if annotation_path is None:
        feature = _from_genbank(record) if fmt == "genbank" else []
    else:
        suffix = Path(annotation_path).suffix.lower()
        if suffix in _GFF_SUFFIX:
            feature = _from_gff(annotation_path)
        elif suffix in _BED_SUFFIX:
            feature = _from_bed(annotation_path)
        else:
            raise ValueError(
                f"{annotation_path} is not a GFF3, GTF or BED file; give the "
                "annotation in one of those, or use a GenBank reference")
        feature = _match_seqid(feature, record.id)

    if alias:
        feature = [
            Feature(entry.start, entry.end, entry.strand,
                    alias.get(entry.name, entry.name), entry.kind, entry.seqid)
            for entry in feature
        ]

    sequence = str(record.seq).upper()
    offset = 0
    if region is not None:
        start, end = region
        offset = max(0, int(start) - 1)
        sequence = sequence[offset:int(end)]

    return Reference(
        identifier=record.id,
        description=record.description,
        sequence=sequence,
        offset=offset,
        full_length=len(record.seq),
        feature=tuple(feature),
    )


def _aligner():
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def _core(consensus_sequence, gap_fraction, core_gap_fraction):
    return np.flatnonzero(np.asarray(gap_fraction) <= core_gap_fraction)


def detect_orientation(consensus_sequence, gap_fraction, reference,
                       core_gap_fraction=0.50, probe=800):
    """Decide whether the panel reads the same way round as the reference.

    An alignment is just a pile of sequence; nothing in it records which
    strand it was deposited on. A panel that is the reverse complement of the
    reference anchors almost nowhere, and every coordinate downstream would be
    wrong, so the two orientation are scored against each other up front.

    Only a probe of the first `probe` core base is scored, because the answer
    is obvious long before the whole panel has been aligned.

    Returns (orientation, forward score, reverse score).
    """
    core = _core(consensus_sequence, gap_fraction, core_gap_fraction)
    if core.size == 0:
        return "forward", 0.0, 0.0
    query = "".join(consensus_sequence[column] for column in core[:probe])

    aligner = _aligner()
    forward = float(aligner.score(reference.sequence, query))
    reverse = float(aligner.score(
        reference.sequence, str(Seq(query).reverse_complement())))
    return ("reverse" if reverse > forward else "forward"), forward, reverse


def map_column(consensus_sequence, gap_fraction, reference, core_gap_fraction=0.50):
    """Map every alignment column to a 0 based reference position, or to -1.

    Only a column whose gap fraction is at or below `core_gap_fraction` can
    anchor, since a column that is mostly gap carries no reliable base. Those
    core column are concatenated into a query and locally aligned to the
    reference; local alignment is used, not global, because the panel may
    cover only part of the reference.

    Returns (mapping, note) where mapping has one entry per column and note
    describes the span matched, for the run log.
    """
    length = len(consensus_sequence)
    mapping = np.full(length, -1, dtype=np.int64)

    core = _core(consensus_sequence, gap_fraction, core_gap_fraction)
    if core.size == 0:
        return mapping, "no column below the core gap fraction; nothing anchored"
    query = "".join(consensus_sequence[column] for column in core)
    best = _aligner().align(reference.sequence, query)[0]

    target_block, query_block = best.aligned
    for (target_start, target_end), (query_start, _) in zip(target_block, query_block):
        for step in range(target_end - target_start):
            mapping[core[query_start + step]] = reference.offset + target_start + step

    matched = mapping[mapping >= 0]
    identity = best.score / max(1, len(query))
    note = (
        f"{reference.identifier} {matched.min() + 1}-{matched.max() + 1} "
        f"({matched.size}/{length} column anchored, score {best.score:.0f} "
        f"over {len(query)} core base, {identity:.2f} per base)"
    )
    return mapping, note


def span(mapping, column):
    """Reference span covered by a set of alignment column, 1 based inclusive.

    Returns (None, None) when no column in the set is anchored.
    """
    position = mapping[np.asarray(column)]
    position = position[position >= 0]
    if position.size == 0:
        return None, None
    return int(position.min()) + 1, int(position.max()) + 1


def structural():
    """Feature type that is bookkeeping rather than something a guide hits."""
    return _SKIP_KIND


def overlapping(reference, start, end):
    """Every annotated feature overlapping a 1 based inclusive span, ordered
    most specific first, meaning the narrowest span first."""
    if start is None or end is None:
        return []
    low, high = start - 1, end  # 0 based half open
    hit = [feature for feature in reference.feature
           if feature.kind not in _SKIP_KIND
           and feature.start < high and low < feature.end]
    return sorted(hit, key=lambda feature: feature.end - feature.start)


def annotate(reference, start, end):
    """Describe what a span hits, as (gene, element).

    `gene` is the gene product the span sits in. With no coding feature there
    it falls back to the broadest enclosing region, so a guide in a non coding
    element still reports where it is rather than nothing.

    `element` is every other overlapping feature, narrowest first, which is
    what names the regulatory part actually being targeted: a guide over the
    HIV-1 LTR TATA box reads as gene "5' LTR", element "TATA box; U3".
    """
    hit = overlapping(reference, start, end)
    if not hit:
        return "", ""

    gene = [feature for feature in hit if feature.kind in _GENE_KIND]
    if gene:
        named, rest = gene, [f for f in hit if f.kind not in _GENE_KIND]
    else:
        named, rest = [hit[-1]], hit[:-1]  # the broadest is the context

    def join(feature):
        seen = []
        for entry in feature:
            if entry.name not in seen:
                seen.append(entry.name)
        return "; ".join(seen)

    return join(named), join(rest)
