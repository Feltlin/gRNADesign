# gRNA Designer

Design a guide RNA against a **population** of a target rather than against a
single reference sequence.

Give the tool an aligned sequence set, the PAM motif of your nuclease, and an
annotated reference. It slides a protospacer window along both strand, keeps
every window that is cuttable across the whole panel, and scores each survivor
on four independent axis. It then reports where each guide lands on the
reference and what it hits, and draws the result as one map.

Nothing in the pipeline is organism specific. The nuclease is whatever PAM
motif you give it, the reference is whatever annotated sequence you point at,
and the panel is any aligned FASTA. It was written for an OSU iGEM project
targeting the integrated HIV-1 provirus, and that case is carried through
section 11 as a worked example, but no step assumes it and none of that data
is part of this repository.


## Contents

1. [Overview](#1-overview)
2. [Biological Background](#2-biological-background)
3. [Install](#3-install)
4. [Usage](#4-usage)
5. [Input](#5-input)
6. [Output](#6-output)
7. [Data Process](#7-data-process)
8. [Coordinate & Strand](#8-coordinate--strand)
9. [File Map](#9-file-map)
10. [Test](#10-test)
11. [Targeting the HIV-1 Provirus](#11-targeting-the-hiv-1-provirus)
12. [Limitation](#12-limitation)
13. [Reference](#13-reference)


## 1. Overview

**Input.** An aligned FASTA panel, a PAM motif, and optionally an annotated
reference.

**Output.** A TSV of candidate guide, each scored on four independent axis and
located on the reference, plus a one figure map of the whole locus.

**What it does, in one pass:**

```
read the alignment
  -> per column frequency, consensus and Shannon entropy
  -> anchor the panel to the reference, checking orientation first
  -> slide a (spacer + PAM) window along both strand
       gate: gap content low, consensus PAM usable, PAM present across the panel
  -> score every survivor: entropy mean, PAM conservation, mean CFD
  -> fold each guide for its minimum free energy          (optional)
  -> report position, gene and regulatory element         (needs a reference)
  -> write the table, and draw the map                    (optional)
```

Everything is driven from the command line, with every default in one file,
`config.py`. Nothing is organism specific: section 11 works the HIV-1 provirus
through as an example, but no step assumes it.


## 2. Biological Background

### Why One Reference Is Not Enough

A guide designed against a single published genome often fails in practice,
because the thing being cut is not one genome. It is a population: the
quasispecies of a virus, a bacterial species with strain variation, a gene
family with paralogue. HIV-1 is the extreme case, mutating fast enough that one
infected person carries a swarm of variant rather than a single sequence.

Cas9 is sensitive enough for that to matter. One mismatch in the wrong place,
or a PAM that has drifted, is enough to stop the cut. So a guide that looks
perfect against the reference can miss most of what is actually there.

Designing against an **alignment** of many real sequence turns the problem into
a question about column: which stretch of this locus is conserved enough, and
tolerant enough, that one guide still cuts nearly everything in the panel.

### What Each Metric Measures

| Metric | The question it answers | Biology behind it |
| --- | --- | --- |
| **Shannon entropy** | How variable is this locus across the panel? | A conserved stretch is usually one the organism cannot easily change, because it is functionally constrained: a structural RNA element, a promoter motif, a catalytic residue. Those are harder places for an escape mutation to appear |
| **PAM conservation** | Does the PAM itself survive across the panel? | No PAM, no cut, however good the spacer. Cas9 needs that motif present in the target DNA, so the PAM is a point of failure separate from the protospacer |
| **CFD score** | Given the variation that *is* there, how much cutting survives? | Not every mismatch matters equally. One far from the PAM is often tolerated; one in the PAM proximal seed region usually is not. CFD puts a published number on that asymmetry |
| **Minimum free energy** | Does the guide fold up on itself? | A spacer that forms a strong hairpin is less available to load into the nuclease and to pair with its target, so it underperforms however good the target site is |

### Why They Are Kept Apart

Entropy and CFD are deliberately **not** fused into a single rank. They can
disagree, and the disagreement is the information: a locus can be variable in a
way CFD tolerates, or conserved in a way one mismatch still ruins. Collapsing
them hides exactly the case worth looking at. PAM conservation and MFE describe
different failure mode again.

Ranking is therefore left to you, because the right trade off depends on the
experiment.


## 3. Install

```bash
pip install -r requirements.txt
```

Required: `numpy`, `pandas`, `matplotlib`, `biopython`.

Optional:

- `ViennaRNA` makes folding local and much faster. Without it, `--mfe` uses the
  public RNAfold web server through the standard library, so nothing extra is
  needed.
- `openpyxl` is needed only by the converter in `tool/`, never by the pipeline.


## 4. Usage

### Command

```
python gRNADesign.py ALIGNMENT [OUT] [option]
```

`ALIGNMENT` is required. `OUT` defaults to `output/gRNA_table.tsv`. Everything
else is optional, and every default lives in `config.py`.

```bash
# A GenBank reference carries its own annotation, so one file is enough
python gRNADesign.py panel.fasta --reference genome.gb

# A FASTA reference with a separate GFF3, GTF or BED
python gRNADesign.py panel.fasta --reference genome.fasta --annotation genome.gff3

# No reference at all: score the panel, leave position and annotation empty
python gRNADesign.py panel.fasta

# A different nuclease
python gRNADesign.py panel.fasta --pam TTTV --pam-side 5prime --spacer-length 23
```

### Every Flag

| Flag | Default | Meaning |
| --- | --- | --- |
| `ALIGNMENT` | required | Aligned FASTA panel |
| `OUT` | `output/gRNA_table.tsv` | Output TSV |
| `--reference` | none | Annotated reference, GenBank or FASTA |
| `--annotation` | none | Separate GFF3, GTF or BED |
| `--reference-region START END` | none | **Restrict** anchoring to one span, 1 based, both end included. Omit to search the whole reference |
| `--pam` | `NGG` | IUPAC PAM motif |
| `--pam-side` | `3prime` | `3prime` (Cas9 like) or `5prime` (Cas12a like) |
| `--spacer-length` | `20` | Protospacer length |
| `--max-gap-fraction` | `0.10` | Gate 1: max gap in any window column |
| `--min-pam-fraction` | `0.80` | Gate 3: min panel fraction with a usable PAM |
| `--core-gap-fraction` | `0.50` | Max gap for a column to be allowed to anchor |
| `--mfe` | off | Fold every guide |
| `--mfe-backend` | `auto` | `auto`, `local` or `server` |
| `--mfe-url` | ViennaRNA server | RNAfold endpoint |
| `--mfe-temperature` | `37` | Folding temperature, °C |
| `--mfe-worker` | `4` | Concurrent server job |
| `--mfe-cache` | `data/mfe_cache.json` | Fold cache |
| `--scaffold` | empty | Appended to the spacer before folding |
| `--plot [PATH]` | off | Draw the candidate map |
| `--entropy-window` | `20` | Running mean window for the entropy line |
| `--cfd-line` | none | Optional guide line on the CFD panel |
| `--consensus [PATH]` | off | Write the panel consensus FASTA |

### A Complete Example

Every flag at once, against the HIV-1 data from section 11. None of those file
are in this repository; section 11 has the command that fetch them.

```bash
python gRNADesign.py \
    input/hiv-db_gap_squeeze.fasta \
    output/gRNA_table.tsv \
    --reference input/K03455.gb \
    --annotation input/HXB2.gff3 \
    --reference-region 1 634 \
    --pam NGG \
    --pam-side 3prime \
    --spacer-length 20 \
    --max-gap-fraction 0.10 \
    --min-pam-fraction 0.80 \
    --core-gap-fraction 0.50 \
    --mfe \
    --mfe-backend auto \
    --mfe-temperature 37 \
    --mfe-worker 4 \
    --mfe-cache data/mfe_cache.json \
    --scaffold "" \
    --plot output/gRNA_map.png \
    --entropy-window 20 \
    --cfd-line 0.8 \
    --consensus output/consensus.fasta
```

### The Same Run, Leaning On The Default

Every flag above whose value is already the default can simply be dropped. For
the HIV-1 case that leaves four, because four thing have no sensible default:
the reference, its annotation, which of the two LTR copy to pin to, and the two
stage that are off unless asked for.

```bash
python gRNADesign.py input/hiv-db_gap_squeeze.fasta \
    --reference input/K03455.gb \
    --annotation input/HXB2.gff3 \
    --reference-region 1 634 \
    --mfe --plot
```

This writes the same `output/gRNA_table.tsv` and `output/gRNA_map.png`. The
PAM `NGG` 3' of a 20 nt spacer, the 0.10 and 0.80 gate, the 20 bp entropy
window and the cache at `data/mfe_cache.json` all come from `config.py`.

Two cosmetic thing differ from the longer command, both because their flag is
off by default: the CFD panel has no dashed guide line, and no
`consensus.fasta` is written. Add `--cfd-line 0.8` or `--consensus` to get
them back.

**On the first run** the fold cache is empty, so `--mfe` submits one job per
distinct guide to the public RNAfold server and the step is the slow part of
the run:

```
Bonus  Fold each guide for its minimum free energy
  31 distinct guide, 0 cached, 31 to fold via server
  folded 10/31
  folded 20/31
  folded 30/31
  folded 31/31
```

That is the whole cost of the run: 31 guide took about two minute against the
public server, while every other step together took under a second. Every later
run reads `data/mfe_cache.json` instead and the step is instant:

```
  31 distinct guide, 31 cached, 0 to fold via server
```

Without `--mfe` the whole run takes well under a second, and the `MFE` and
`Fold_Structure` column are left empty. With a local ViennaRNA install there is
no server and no wait at all.


## 5. Input

### The Alignment

**The panel must already be aligned.** This tool does not align for you. Run
MAFFT, MUSCLE, Clustal or whatever your field uses first, and feed in the
result.

That is not a formality. Every number here is computed **down a column**: the
entropy of column $j$, the gap fraction of column $j$, the consensus base at
column $j$. That only means anything if column $j$ holds the same position in
every sequence. Alignment is what makes a column a position rather than an
index.

What happens if you feed in something unaligned:

- **Row of different length** → refused outright, naming the length it saw.
  This catches an ordinary unaligned FASTA, since real sequence rarely comes
  out exactly the same length.
- **Row of the same length that are not actually aligned** (for example every
  sequence trimmed to 600 nt) → this passes the length check, so the run
  **warns** instead:

  ```
  WARNING: mean column entropy is 0.41 and there is essentially no gap.
  That combination usually means the input is not actually aligned, only
  equal in length.
  ```

  The signature is high per column entropy together with no gap at all. A real
  alignment of related sequence either carries gap from the indel it absorbed,
  or is conserved enough to have low entropy. Carrying neither is the mark of
  sequence that were never aligned. It is a warning and not an error, because
  a genuinely very diverse panel can look the same; check the panel rather than
  the message.

Gap is expected and handled. An IUPAC ambiguity code is skipped rather than
guessed at: it shrinks the sample at that column instead of being counted as a
base it might not be.

**How the alignment was trimmed does not matter.** Anchoring is by sequence
content, not by column index, so a gap stripped, gap squeezed or sub-region
alignment all work and all report the same coordinate. See
[section 7.3](#73-step-3-anchor-the-panel-to-the-reference).

### The Reference

Optional. It supplies **coordinate and annotation only**; every score is
computed from the panel. Without `--reference` the run still produces a full
table, with `Ref_Start`, `Ref_End`, `Gene` and `Element` left empty, and
`--plot` is skipped because the map is drawn against reference coordinate.

### Which Annotation Format Work

The tool needs a reference **sequence** and **feature** with coordinate on it.
They can come from one file or two.

| You have | Pass | Note |
| --- | --- | --- |
| GenBank (`.gb`, `.gbk`, `.gbff`, `.genbank`) | `--reference` | Sequence and annotation in one file. Simplest case |
| FASTA + GFF3 (`.gff`, `.gff3`) | `--reference` + `--annotation` | The usual pairing from Ensembl, NCBI or a browser |
| FASTA + GTF (`.gtf`) | `--reference` + `--annotation` | Read as GFF3 with GTF style attribute |
| FASTA + BED (`.bed`) | `--reference` + `--annotation` | Minimal: every feature becomes a gene |
| FASTA alone | `--reference` | Coordinate reported, `Gene` and `Element` empty |

`--annotation` **replaces** whatever annotation the reference carried rather
than adding to it, so that giving a GenBank file and its own GFF3 export does
not duplicate every feature.

Things worth knowing before you hit them:

- **Coordinate convention is handled.** GFF3 and GTF are 1 based inclusive,
  BED is 0 based half open; both land on the same internal span, and there is
  a test for exactly this.
- **The sequence name need not match.** If no feature matches the reference
  name, the file is read as being about that one sequence rather than silently
  coming back empty.
- **Only the first record of the reference is read.** Give a single sequence
  reference, or split it first.
- **BED carries no feature type**, so every BED feature becomes a gene and
  nothing reaches `Element`. Use GFF3 if you want the regulatory element split
  out.
- **A compound location is split into its part**, so a guide in an intron is
  not reported as hitting the gene.
- **Feature naming** comes from `Name`, then `gene`, `gene_name`, `product`,
  `locus_tag`, `ID`, `note`. First one present wins.
- **`feature_alias`** in `config.py` substitutes display name for a record
  whose feature name are archaic or terse. It never moves a coordinate.


## 6. Output

### The Table

One TSV row per surviving candidate.

| Column | Meaning |
| --- | --- |
| `gRNA5to3` | `SPACER(PAM)`, read 5' to 3' in the guide's own frame |
| `Derived_Strand` | `+` or `-`; the strand the protospacer is read off |
| `Entropy_Mean` | Mean Shannon entropy over the spacer column |
| `PAM_Conservation` | Fraction of the panel carrying that exact PAM |
| `CFD_Score` | Mean CFD of the guide against every sequence in the panel |
| `MFE` | Minimum free energy of the folded guide, kcal/mol |
| `Fold_Structure` | That fold in dot bracket notation |
| `Gene` | Gene product hit, or the enclosing region when there is none |
| `Element` | Narrower feature hit, narrowest first |
| `Ref_Start`, `Ref_End` | Protospacer span on the reference, 1 based inclusive |
| `Align_Start`, `Align_End` | The same span as alignment column, for tracing a row back to the panel |

#### What To Order

`gRNA5to3` reads `GACAAGATATCCTTGATCTG(TGG)`. **Order the 20 nt outside the
bracket.** The bracketed PAM is not part of the guide: it is the motif that has
to be present in the target DNA, immediately 3' of the protospacer. It is shown
because a guide is meaningless without knowing its PAM is there, and bracketed
so it can never be pasted into an order by mistake.

Both parts are written 5' to 3' in the guide's own frame, so a minus strand
guide reads just as directly as a plus strand one. The sequence shown is always
the sequence to synthesise, whatever `Derived_Strand` says.

#### How To Read The Coordinate

All four coordinate are **plus strand of the reference**, and **`Start` is
always the smaller number**, on both strand. `Ref_*` and `Align_*` describe the
same 20 nt protospacer, one in reference coordinate and one in alignment
column, so they can be compared directly.

Start is deliberately *not* the 5' end of the guide. Every genomic format
(GFF3, BED, GenBank, SAM, VCF) reports `start ≤ end` with orientation carried
by a separate strand field, and following that convention is what lets the
output go straight into bedtools, IGV or a genome browser without being
rewritten. Reversing the coordinate on minus strand row would break all of
them.

To get the 5' end of the guide itself, read `Derived_Strand`:

| `Derived_Strand` | 5' end of the guide | 3' end, next to the PAM |
| --- | --- | --- |
| `+` | `Ref_Start` | `Ref_End` |
| `-` | `Ref_End` | `Ref_Start` |

### Reading The Number

| Column | Range | Lowest means | Highest means | Want |
| --- | --- | --- | --- | --- |
| `Entropy_Mean` | 0 to 1 | **0**: every sequence identical across all 20 nt | **1**: all five symbol equally likely at every position | **low** |
| `PAM_Conservation` | 0 to 1 | **0**: not one sequence carries the consensus PAM exactly | **1**: every sequence carries it | **high** |
| `CFD_Score` | 0 to 1 | **0**: the guide cuts nothing in the panel | **1**: cuts every sequence at full efficiency | **high** |
| `MFE` | ≤ 0 | strongly negative, about −8 or below for a 20 mer: folds hard on itself | **0**: no base pair at all, the spacer is free | **near 0** |

Two caveats about the extremes.

`Entropy_Mean` is on a log base 5 scale over `{A, C, G, T, -}`, so 1 is the
entropy of a column where all five symbol are equally likely. A real alignment
rarely passes about 0.5. Compare candidate to each other, not to the
theoretical 1.

`CFD_Score` is a mean over the whole panel **including** the sequence the guide
fails on, so 0.90 does not mean "90% efficient". It means closer to "the
average sequence in this panel is cut at 90% of the efficiency of a perfect
match". One sequence with a broken PAM contributes a hard 0 and pulls the mean
down, which is intended.

### The Plot

`--plot` draws one figure, four panel sharing the reference coordinate:

1. **CFD score** per candidate as a stem. Height is the score, marker area is
   the MFE, colour is the strand.
2. **Shannon entropy** as a running window mean line: the conservation
   landscape. A dip is a conserved stretch.
3. **Reference annotation**. A broad feature is a labelled bar, broadest at the
   bottom; a short element is a tick with a leader up to its label.
4. **Guide footprint** as an arrow pointing 5' to 3' along the strand the
   protospacer is read from. Plus strand above the rule, minus below. Each
   strand is packed into as many sub lane as it takes for no two arrow to
   touch, so a cluster of guide a base apart stays countable.

The strand colour is shared between panel 1 and panel 4, so a stem can be read
straight down onto its own arrow.

Design follows the Swiss International Typographic manner: one grid, flush left
text, a single sans serif in two weight, hairline rule, ink and paper with two
muted accent. Every mark carries data. There is no threshold line on the CFD
panel by default, because the y axis already states the score; pass
`--cfd-line 0.8` if you want one drawn.


## 7. Data Process

Notation. The alignment has $N$ sequence and $L$ column. $x_{k,j}$ is the
symbol of sequence $k$ at column $j$. The alphabet counted is

$$\Sigma = \lbrace \mathtt{A}, \mathtt{C}, \mathtt{G}, \mathtt{T}, \mathtt{-} \rbrace .$$

### 7.1 Step 1: Read The Alignment

The panel is read into an $N \times L$ matrix of byte, which makes every step
below a vector operation rather than a loop over sequence. Unequal row length
is an error, not a warning: a ragged "alignment" silently destroys every column
statistic that follows.

The minus strand is handled by reverse complementing the whole matrix once and
running the identical search on it, so nothing downstream ever reads a row
backward. Complementing is IUPAC aware: gap maps to gap, and an ambiguity code
maps to its complementary code ($\mathtt{R} \leftrightarrow \mathtt{Y}$,
$\mathtt{K} \leftrightarrow \mathtt{M}$, $\mathtt{N} \leftrightarrow \mathtt{N}$).

### 7.2 Step 2: Column Statistic

For each column $j$, count every symbol of $\Sigma$ and normalise:

$$p_{j,s} = \frac{n_{j,s}}{\sum_{s' \in \Sigma} n_{j,s'}}, \qquad s \in \Sigma .$$

A symbol outside $\Sigma$ matches no counter, so an ambiguity code never enters
the denominator.

**Gap fraction** is the gap component, $g_j = p_{j,\mathtt{-}}$.

**Consensus** is the majority base, ignoring gap and ambiguity, with a tie
broken toward A, C, G, T order so a rerun gives the same answer:

$$c_j = \arg\max_{b \in \lbrace \mathtt{A}, \mathtt{C}, \mathtt{G}, \mathtt{T} \rbrace} p_{j,b} .$$

**Shannon entropy** over all five symbol, in base 5:

$$H_j = -\sum_{s \in \Sigma} p_{j,s} \log_5 p_{j,s}, \qquad 0 \log 0 := 0 .$$

Base 5 is not cosmetic. A column with all five symbol equally likely has
entropy exactly 1 and a pure column exactly 0, so $H_j \in [0,1]$ by
construction and is never renormalised later. Counting the gap as a symbol is
also deliberate: a column that is half gap is genuinely unreliable to target,
and that should show up as entropy rather than be hidden.

Entropy is strand symmetric. Complementing permutes the symbol and reversing
permutes the column, and neither changes a sum over the whole alphabet, so $H$
is computed once on the plus strand and reused for a minus strand hit.

### 7.3 Step 3: Anchor The Panel To The Reference

This is the step most likely to be silently wrong, so it is worth
understanding.

**The problem.** An alignment is a pile of sequence in its own coordinate
system. Column 1 of your panel is column 1 of your panel; it is not base 1 of
anything. Nothing inside a FASTA records which part of which genome it came
from, or even which strand it is on. So no correspondence is assumed. One is
derived, once, and printed.

**How.** Only a **core column**, one with $g_j \le$ `--core-gap-fraction`
(0.5), may anchor, since a column that is mostly gap carries no reliable base.
The consensus letter of the core column are concatenated into one string and
**locally** aligned to the reference, with an affine gap penalty:

$$s_{\text{match}} = +1, \quad s_{\text{mismatch}} = -1, \quad s_{\text{open}} = -10, \quad s_{\text{extend}} = -0.5 .$$

Local, not global, because the panel usually covers only part of the reference.
The alignment yields

$$m : \lbrace 0, \dots, L-1 \rbrace \to \lbrace 0, \dots, R-1 \rbrace \cup \lbrace -1 \rbrace,$$

mapping each column to a reference position, or to $-1$ when it anchors
nothing.

**Worked micro example.** Four row, eight column:

```
column          1  2  3  4  5  6  7  8
row 1           A  C  -  G  T  A  C  G
row 2           A  C  -  G  T  A  C  G
row 3           A  T  -  G  T  A  -  G
row 4           A  C  -  G  T  A  C  G
gap fraction   .0 .0 1.0 .0 .0 .0 .25 .0
consensus       A  C  (-) G  T  A  C  G
core column     *  *      *  *  *  *  *
```

Column 3 is all gap, so it is not a core column and cannot anchor. The query
built from the core column is `ACGTACG`. Suppose that aligns to reference
positions 101–107. Then:

```
column    1    2    3    4    5    6    7    8
maps to  101  102   -   103  104  105  106  107
```

Column 3 reports no coordinate; everything else does. A guide sitting on
columns 4–8 is reported as reference 103–107.

**Why the trimming of your alignment does not matter.** Anchoring is driven by
sequence content, never by column index. So:

- A **gap squeezed** alignment (LANL's Gap Strip/Squeeze, or any tool that
  drops mostly gap columns) is shorter and has *fewer* unanchored column, so it
  anchors more cleanly, not less. The coordinate it reports are the same.
- An alignment with a long **insertion** relative to the reference anchors
  everywhere except the inserted column, which report no coordinate.
- An alignment of a **sub-region** is found by the local alignment wherever it
  sits in the reference.
- A **differently gapped** alignment of the same sequence gives the same
  reference coordinate, because the gap column were never the thing being
  matched.

**What is printed, every run:**

```
anchored to K03455.1 1-634 (634/866 column anchored, score 590 over 634
core base, 0.93 per base)
```

Read that line. It names the reference span matched, how many column found a
home, and the per base match rate.

**The four way it can go wrong.**

1. **The panel is reverse complemented.** Nothing in a FASTA says which strand
   it holds. Both orientation are scored against a short probe before anything
   else, and the whole panel is flipped if the reverse wins, so a flipped panel
   gives byte for byte the same table as a forward one.
2. **The target is repeated.** A local alignment cannot choose between two
   identical copy, so the coordinate can land on the wrong one with a perfectly
   healthy score. `--reference-region START END` removes the choice. Any
   repeated element, gene family or segmental duplication has this problem.
   See "What `--reference-region` Is For" just below.
3. **The panel covers only part of the reference.** Normal, expected and
   handled. This is the whole reason the alignment is local rather than global.
   A 511 nt panel inside a 9719 nt genome anchors exactly, and the anchored
   span tells you where. You do **not** need any flag for this.
4. **The panel and the reference are not the same locus.** Then little anchors,
   and the run warns when under half the column do.

#### What `--reference-region` Is For

It **restricts** the search. It does not define it.

- **Omit it** and the whole reference is searched. That is the normal case,
  including when your panel covers only a small part of the reference. You do
  not need the flag to "select" a sub-region of the genome; the local alignment
  finds your panel wherever it sits.
- **Give it** only when the panel matches more than one place and you must say
  which. That is the repeat case.

`START` and `END` are 1 based with **both end included**, so `1 634` is 634
base, and the reported coordinate stays in **full reference space** regardless.
Pinning to `9086 9719` reports a hit there as 9100, not as 15.

For HIV-1 the flag is needed because the LTR is in HXB2 twice, at 1–634 (the
5' LTR, which is where `1 634` comes from) and at 9086–9719 (the 3' LTR), and
the two are near identical. Without the flag:

```
anchored to K03455.1 9086-9719   # the 3' copy, score 592
```

With `--reference-region 1 634`:

```
anchored to K03455.1 1-634       # the 5' copy, score 590
```

592 against 590. The margin is two points out of 634, so which copy wins is
effectively arbitrary, and the reported coordinate differ by about 9 kb. That
is exactly the situation the flag exists for.

By contrast, a panel of a **unique** locus needs nothing. A 511 nt panel taken
from *gag* anchors on its own:

```
anchored to K03455.1 790-1300 (511/511 column anchored, 1.00 per base)
```

### 7.4 Step 4: Window Enumeration And Gate

Let the spacer be $\ell$ nt and the PAM motif $M$ be $m$ nt, giving a window of
$\ell + m$. For a 3' PAM the window is spacer then PAM; for a 5' PAM it is PAM
then spacer. Every start position on both strand is tested against three gate.

**Gate 1, gap.** Every column of the window must be clean enough:

$$\max_{j \in \text{window}} g_j \le g_{\max} .$$

**Gate 2, consensus PAM.** The consensus must itself present a usable PAM,
written $c_{\text{PAM}} \models M$ for "satisfies the motif at every position".

**Gate 3, PAM presence.** The site must be cuttable across the panel, not only
in the consensus:

$$f_{\text{PAM}} = \frac{1}{N} \left| \lbrace k : x_{k,\text{PAM}} \models M \rbrace \right| \ge p_{\min} .$$

Gate 3 is a **motif** match, not an exact match: a variant reading `TGG` where
the consensus reads `AGG` is still cut by SpCas9, so it counts as present.
That is a different question from the one `PAM_Conservation` answers, which is
why both exist.

### 7.5 Step 5: Metric

#### Shannon Entropy Mean

A plain, unweighted mean over the spacer column:

$$\bar{H} = \frac{1}{\ell} \sum_{j \in \text{spacer}} H_j \in [0, 1] .$$

Unweighted on purpose: position dependent weighting is what CFD already does,
with a published and validated weight, so doing it twice would double count the
seed region and make the two column correlate by construction.

#### PAM Conservation

The fraction of the panel carrying the **exact** consensus PAM:

$$C_{\text{PAM}} = \frac{1}{N} \left| \lbrace k : x_{k,\text{PAM}} = c_{\text{PAM}} \rbrace \right| .$$

Since an exact match implies a motif match, $C_{\text{PAM}} \le f_{\text{PAM}}$
always. A large gap between the two says the PAM is present but its identity is
drifting, which matters wherever PAM identity affects efficiency.

#### CFD Score

The Cutting Frequency Determination model of Doench et al. 2016, as used by
CRISPOR. For the consensus guide $g$ against one variant $v^k$:

$$\mathrm{CFD}(g, v^k) = P(\pi^k) \prod_{i=1}^{\ell} a_i(g_i, v^k_i) ,$$

where $a_i$ is the activity retained by a mismatch at spacer position $i$,

$$a_i(x, y) = \begin{cases}
1 & x = y \\
0 & y \in \lbrace \mathtt{-} \rbrace \cup \text{ambiguity} \\
T\left[\, \mathtt{r}x : \mathtt{d}\bar{y},\, i \,\right] & \text{otherwise,}
\end{cases}$$

with $T$ the published mismatch table and $\bar{y}$ the complement of $y$. A
gap scores 0 rather than being looked up, because an indel shifts the guide out
of frame and the table has no notion of an indel.

$P$ is the PAM factor. For SpCas9 `NGG` it is the published PAM table on the
two 3' base; any other motif falls back to a presence test, because the CFD PAM
table is SpCas9 specific:

$$P(\pi) = \begin{cases}
T_{\text{PAM}}[\pi_{2:3}] & M = \mathtt{NGG} \\
\mathbb{1}[\pi \models M] & \text{otherwise.}
\end{cases}$$

The reported score is the mean over the whole panel:

$$\overline{\mathrm{CFD}} = \frac{1}{N} \sum_{k=1}^{N} \mathrm{CFD}(g, v^k) \in [0, 1] .$$

The denominator is the full $N$. A variant with a gap in the protospacer or a
broken PAM contributes 0 and still counts, so it pulls the mean down rather
than quietly disappearing from the average.

Because the published table is indexed over position 1 to 20 of an SpCas9
spacer, CFD is only defined for $\ell = 20$. At any other spacer length the
column is left empty and the run says so, rather than extrapolating.

### 7.6 Step 6: Fold

Each **distinct** guide is folded as RNA and its minimum free energy reported
with the structure in dot bracket notation. Folding uses a local ViennaRNA
install when one is importable and the public RNAfold web server otherwise, at
`--mfe-temperature` with the Turner 2004 parameter. By default the bare spacer
is folded; `--scaffold` appends the constant region so the whole sgRNA is
folded instead.

**One sequence per job.** The RNAfold server silently concatenates a multi
record FASTA into a single sequence and returns one quietly wrong answer, so
one job is submitted per distinct guide and nothing is ever batched.

**The cache.** Result are stored in `data/mfe_cache.json` (`--mfe-cache`), a
flat JSON object:

```json
{"AAGCACCATCCAAAGGTCAG||37": ["....(((.......)))...", -1.1]}
```

The key is `sequence|scaffold|temperature`, so changing the scaffold or the
temperature is a different key and refolds, while a rerun at the same setting
costs nothing and never touches the network. Two guide with the same spacer
fold once. Delete the file to force a full refold.

### 7.7 Step 7: Annotate

The protospacer column are mapped back through $m$ from step 3, and the span is
reported 1 based inclusive:

$$\text{Ref\_Start} = \min_{j \in \text{spacer},\, m(j) \ge 0} m(j) + 1, \qquad
\text{Ref\_End} = \max_{j \in \text{spacer},\, m(j) \ge 0} m(j) + 1 .$$

Every annotated feature overlapping that span is collected, sorted narrowest
first, and split in two:

- **`Gene`** is the overlapping gene product (`CDS`, `gene`, `mat_peptide`).
  With no coding feature there it falls back to the single broadest overlapping
  feature, so a guide in a non coding element still reports where it sits.
- **`Element`** is everything else that overlaps, narrowest first. This is the
  column that names the regulatory part being targeted.

Structural bookkeeping (`source`, `mRNA`, `prim_transcript`, `exon`, `intron`)
is skipped, because "this guide is inside the genomic mRNA" is true of almost
everything and says nothing useful.


## 8. Coordinate & Strand

Two frame are kept deliberately separate, and confusing them is the easiest way
to order the wrong oligo.

The **search frame** is whichever matrix the window was found in: the plus
matrix, or the reverse complemented matrix for a minus strand hit. Scoring
reads this frame, because the guide must be read 5' to 3' for the CFD position
index to mean anything.

The **report frame** is always the plus strand of the reference. Reverse
complemented column $j$ is plus column $L - 1 - j$, so a minus strand window
carries both index set and every coordinate in the table is plus frame.

`Derived_Strand` is a "read off" label, never a "binds to" label. `+` means the
protospacer and its PAM are read off the plus strand, so the guide RNA you
synthesise has the sequence printed in `gRNA5to3`, and the RNA itself pairs with
the minus strand. **The `gRNA5to3` column is always the sequence to order**, minus
the bracketed PAM, whatever the strand says.


## 9. File Map

| File | Job |
| --- | --- |
| `gRNADesign.py` | Command line entry; orchestrates the pipeline and owns the run log |
| `config.py` | Every default, in one place. The only file to edit to change behaviour |
| `alignment.py` | FASTA input, the byte matrix, reverse complement, column frequency, consensus, entropy, running mean |
| `pam.py` | IUPAC motif matching, window enumeration, the three gate, PAM conservation |
| `cfd.py` | The published CFD mismatch and PAM table, the score, and the vectorised activity lookup |
| `score.py` | The per candidate metric: entropy mean and mean CFD |
| `annotation.py` | Reference reader (GenBank, GFF3, GTF, BED), orientation detection, the column to coordinate map, the gene and element lookup |
| `mfe.py` | Folding: local ViennaRNA or the RNAfold server, with the on disk cache |
| `report.py` | Column name, rounding, sorting, and writing the TSV |
| `plot.py` | The candidate map |

| Folder | Holds |
| --- | --- |
| `data/` | The CFD mismatch and PAM table, and the fold cache. The only data shipped |
| `input/` | Your sequence, reference and annotation. Not tracked by git |
| `output/` | Regenerated by every run. Not tracked by git |
| `tool/` | Optional helper, not part of the pipeline |
| `test/` | The test suite |


## 10. Test

```bash
python -m unittest discover -s test -t .
```

58 test, well under a second, no network touched and no outside data needed:
the end to end test builds its own synthetic panel and reference, so it runs on
a fresh clone.

**`test/test_core.py`** — the maths.

- *Reverse complement* is its own inverse, maps gap to gap, and maps an IUPAC
  code to its complementary code.
- *Entropy* is 1 for a five symbol uniform column and 0 for a pure column,
  which pins the log base; an ambiguity code shrinks the denominator; entropy
  is strand symmetric, which is what licenses computing it once.
- *The gate* each reject and pass on a toy alignment built to sit either side
  of the threshold, and `PAM_Conservation` counts exact consensus match while
  the gate counts motif match.
- *CFD* reproduces the PAM score on a perfect match, scores a gap and a broken
  PAM as 0, and the lookup table agrees with the scalar function everywhere.
- *The index direction guard.* Position 20 is PAM proximal, so the seed must
  come out as the **harshest** end. If the CFD table were indexed backward the
  seed would look like the most tolerant end, and every score in the project
  would be quietly wrong while still looking reasonable. This is the single
  most valuable test here.

**`test/test_annotation.py`** — the coordinate.

- GFF3 is 1 based inclusive and BED is 0 based half open; both must land on the
  same internal span.
- An alias renames without moving a coordinate.
- A gapped panel maps onto the reference with the fully gapped column anchoring
  nothing and the rest strictly increasing.
- `Element` is ordered narrowest first.
- Orientation is detected both way on a non palindromic sequence.
- A duplicated region is reachable by either copy through `--reference-region`,
  and a coordinate stays in full reference space.

**`test/test_integration.py`** — the whole pipeline on a synthetic panel:
every guide a clean 20 mer, every PAM matching the motif, every score in range,
the table round tripping, two run identical, a reverse complemented panel
giving the same answer, and a run without a reference still scoring everything.
It also parses real RNAfold output text without the network, and checks the
cache key behaviour.

**`test/synthetic.py`** — builds the panel, reference and annotation the
integration test uses.


## 11. Targeting the HIV-1 Provirus

The tool is general, but it was built for one job: finding a guide that cut the
integrated HIV-1 provirus across the diversity of circulating strain. This
section is a complete worked example and the source of the default in the
command of [section 4](#a-complete-example).

The LTR is the target worth attacking. It sits at both end of the provirus, so
one guide cut twice, and it holds the promoter the virus needs to transcribe
itself.

**None of this data is in the repository.** It is large, not ours to
redistribute, and it changes upstream. Fetch it into `input/`:

```bash
# 1. The reference genome: HIV-1 HXB2, GenBank K03455, 9719 nt
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=K03455.1&rettype=gbwithparts&retmode=text" \
    -o input/K03455.gb

# 2. The annotation: the LANL HXB2 landmark spreadsheet, then convert it
curl -L https://www.hiv.lanl.gov/content/sequence/HIV/MAP/HXB2.xlsx \
    -o input/HXB2.xlsx
pip install openpyxl
python tool/lanl_hxb2_gff3.py input/HXB2.xlsx input/HXB2.gff3

# 3. The panel: an LTR alignment from the LANL HIV sequence database,
#    https://www.hiv.lanl.gov/ -> Sequence Database -> Curated Alignments.
#    Any aligned FASTA of the same region works.
```

`--reference-region 1 634` pins the panel to the **5'** LTR, which is where
that number comes from: K03455 annotates `repeat_region 1..634` as the 5' LTR,
and position 1 of the record is its first base. The 3' LTR is the same 634 nt
again at 9086–9719. Without the flag the anchor lands on whichever copy scores
a hair higher and the coordinate come back about 9 kb off; see
[section 7.3](#what---reference-region-is-for).

### Why The Annotation Needs A Separate Source

K03455 is the 1985 submission. It carries the sequence, eight CDS under the ORF
name of the day (`sor 23K protein` is what we now call *vif*), and each LTR as
one undivided 634 nt block. It has no U3, R or U5 split, and no TATA box,
NF-κB site, Sp1 site or TAR. NCBI's own GFF3 export of K03455 carries the same
thin content, and the RefSeq record `NC_001802` begins at the transcription
start site, so it has no U3 and no promoter element at all.

The LANL landmark spreadsheet is the authoritative source for those. It gives
one row per base of HXB2 with explicit start and end marker, and
[`tool/lanl_hxb2_gff3.py`](tool/lanl_hxb2_gff3.py) converts it to GFF3,
recording in each line which sheet column it came from. It yields 77 feature:

| Element | HXB2 position |
| --- | --- |
| U3 (5') | 1–455 |
| TCF-1α | 315–329 |
| NF-κB II, NF-κB I | 350–359, 364–373 |
| Sp1 III, Sp1 II, Sp1 I | 375–386, 388–397, 398–408 |
| TATA box | 427–431 |
| TAR, TAR bulge, TAR loop | 453–513, 476–478, 483–488 |
| +1 transcription start | 456 |
| R (5') | 456–551 |
| poly(A) signal | 527–532 |
| U5 (5'), U5 binding region | 552–634, 560–570 |
| PBS | 636–653 |

plus the Psi packaging stem-loop, the RRE and its stem, the polypurine tract,
the Gag, Pol and Env mature peptide (p17, p24, p7, protease, RT, RNase,
integrase, gp120, gp41), and both Tat and Rev exon.

### What LANL Does Not Publish

A drawn LTR map often also shows AP-1, NFAT, COUP-TF, C/EBP, USF-1,
Ets-1/RBF-1, RBEIII, NRT-1, GATA and NF-1/CREB-ATF across the modulatory region
in the 5' half of U3. **Those are not in the LANL landmark spreadsheet, and not
in the LANL HXB2 Feature Database either** — its full 358 name feature
vocabulary contains NF-κB, Sp1, TATA, TCF-1α and the PBS, and none of the
others. They come from the primary literature as interval relative to the
transcription start, not as a coordinate set anyone distributes.

So they are left out rather than guessed at. If you have a citable coordinate,
add it as an extra GFF3 line and it will be picked up and drawn with everything
else:

```
K03455.1	yourcitation	protein_binding_site	<start>	<end>	.	+	.	Name=AP-1 I
```

### A Provenance Caution

LANL and the GenBank record disagree by one base on where R begins: LANL marks
the transcription start at 456, the GenBank mRNA feature begins at 455. The
converted GFF3 follows LANL throughout, since mixing the two would be worse
than either. Expect a one base offset against a source that follows GenBank.


## 12. Limitation

- **No off-target search.** Every score here comes from the panel alignment
  only. A guide can be perfectly conserved in the target and still cut the host
  genome. Screen a shortlist with Cas-OFFinder or CRISPOR against the host
  genome before ordering anything.
- **The panel is the population.** Every metric inherits whatever bias the
  alignment carries. A panel over-representing one subtype reports conservation
  for that subtype, not for the species.
- **CFD is a model.** It predicts relative cutting from a mismatch table, not
  from your experiment, and it was trained on SpCas9. For another nuclease,
  entropy, PAM conservation, position and annotation stay valid; CFD's mismatch
  term does not, and at a spacer length other than 20 it is not computed at
  all.
- **MFE is the guide alone by default.** The spacer is folded without the
  scaffold unless `--scaffold` is given, and without the target.
- **Annotation is only as good as its source.** The tool reports what the file
  says. It does not validate a coordinate against the sequence.


## 13. Reference

- Doench, J. G. et al. Optimized sgRNA design to maximize activity and minimize
  off-target effects of CRISPR-Cas9. *Nature Biotechnology* **34**, 184–191
  (2016). The CFD model.
- Lorenz, R. et al. ViennaRNA Package 2.0. *Algorithms for Molecular Biology*
  **6**, 26 (2011). The folding model behind RNAfold.
- Gruber, A. R. et al. The Vienna RNA Websuite. *Nucleic Acids Research* **36**,
  W70–W74 (2008). The server used when ViennaRNA is not installed locally.
- Mathews, D. H. et al. Incorporating chemical modification constraints into a
  dynamic programming algorithm for prediction of RNA secondary structure.
  *PNAS* **101**, 7287–7292 (2004). The Turner 2004 parameter.
- Ratner, L. et al. Complete nucleotide sequence of the AIDS virus, HTLV-III.
  *Nature* **313**, 277–284 (1985). GenBank K03455.
- Los Alamos National Laboratory HIV Sequence Database,
  [HXB2 landmark spreadsheet](https://www.hiv.lanl.gov/content/sequence/HIV/MAP/HXB2.xlsx).
  The HIV-1 annotation used in section 11.
