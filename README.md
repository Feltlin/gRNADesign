# HIV-1 LTR gRNA Designer

This tool searches for SpCas9 guide RNAs targeting the HIV-1 LTR (long
terminal repeat, a sequence that appears identically at both ends of every
integrated HIV provirus) that stay effective across a 365-sequence panel of
real HIV-1 subtype-B sequences, rather than being designed against a single
reference genome and hoping it generalizes. It also optionally checks
candidates against the human genome to flag guides that might also cut
somewhere they shouldn't.

No biology background is assumed anywhere in this document — every term is
defined the first time it's used, and section 4 formalizes every score with
an explicit equation and variable list.

## 1. What it needs to run

The tool needs three things on disk before it can do anything:

- **An aligned FASTA panel** — `input/hiv-db.fasta` by default.
  "Aligned" means all 365 sequences have already been stretched/padded (with
  `-` gap characters) so that column number `i` refers to the same
  biological position across every sequence — this alignment step is
  assumed already done; this tool does not do any alignment itself. Each
  sequence is 866 columns long. Besides the four DNA letters and gaps, a
  handful of columns contain IUPAC ambiguity codes (`M, Y, R, N, K, W`) — a
  sequencer's way of saying "couldn't confidently call a single base here."
- **Two small CFD lookup tables** — `data/mismatch_score.pkl` and
  `data/pam_scores.pkl`. CFD ("Cutting Frequency Determination") is a
  published, experimentally-measured model (Doench et al. 2016 — the same
  numbers CRISPOR uses) of how much a mismatch or an imperfect PAM reduces
  Cas9 cutting. These tables are the actual measured numbers, not something
  this project invented — reproduced in full in section 7 so you can see
  exactly where every score comes from. The loader tells the two files
  apart by looking at their key structure rather than trusting the
  filename, so it doesn't matter which order you hand them in.
- **(Optional) a reference genome** — only needed if you want the off-target
  step (section 3) to actually run a search; without it, everything else
  in the tool works fine, and the off-target columns are just left blank.

The FASTA panel lives under `input/`, the two CFD tables live under `data/`,
and every run's output (`gRNA_table.tsv`, `plots/`) lands under `output/` by
default — all three paths are set once in `config.py` and used everywhere
else, so moving any of these directories just means updating `config.py`,
not hunting through the codebase.

## 2. What the tool actually does, step by step

1. **Load the alignment, then mirror it.** `io_fasta.py` reads the FASTA,
   then builds its full reverse complement (every base swapped to its
   Watson–Crick partner — A↔T, C↔G — and the sequence order flipped). Guides
   on the minus strand are found by re-running the *exact same* search logic
   against this mirrored alignment, rather than by trying to read a sequence
   backwards in place, which is a much easier place to introduce a subtle
   bug.
2. **Compute per-column statistics, on both alignments separately.**
   `columns.py` computes, for every column: the base frequencies and gap
   fraction, a majority-vote consensus base, and Shannon entropy — a
   standard measure of how mixed a column is (0 = every sequence agrees,
   higher = more disagreement). Ambiguity codes are excluded from a
   column's denominator rather than being treated as a real base.
3. **Load the CFD tables.** `cfd.py` figures out which pickle is which,
   then exposes the atomic scoring function and a derived "how much does a
   mismatch hurt at each of the 20 spacer positions" weight vector — pulled
   out of the real data, not hand-tuned.
4. **Slide a 23-column window across both strands and filter it down.**
   `windows.py` considers every possible 20-nt spacer + 3-nt PAM window
   (`NGG` — SpCas9 requires this exact 3-letter pattern immediately next to
   any site it cuts; `N` = any base, so only the two `G`s are fixed). A
   window is thrown out immediately, before any scoring, if:
   - **the gap gate fails** — any of its 23 columns has more gaps than
     `g_max` allows (default 10%), meaning too much of the panel doesn't
     actually have a real base there; or
   - **the PAM gate fails** — fewer than `p_pam_min` of the panel (default
     80%) actually carries `GG` at the PAM position, meaning Cas9 couldn't
     cut most of the swarm here even with a perfect spacer.
5. **Score every window that survives.** `score.py` computes a
   seed-weighted conservation score and two different estimates of expected
   cutting activity across the panel (see section 4 for exact formulas).
6. **Assemble the results table.** `table.py` merges both strands' surviving
   windows onto one shared coordinate system (section 6) and writes one row
   per guide.
7. **(Optional) off-target scoring.** `offtarget.py`, wired up through
   `--genome` / `--offinder-results` (section 3), only runs if you ask for
   it.
8. **(Optional) diagnostics on a shortlist.** `diagnostics.py` — meant to be
   called by hand on guides you've already picked out, not run automatically
   on the whole table.
9. **(Optional) plots.** `plotting.py` — see the very last section of this
   document.

Four different questions are kept deliberately separate the whole way
through, and never combined into a single score: how *conserved* a region is
(Shannon entropy), how much a *specific mismatch* is tolerated (CFD), how
*efficiently* Cas9 cuts a perfect match in the first place (not modeled here
— that's a separate tool like Azimuth/DeepSpCas9), and whether a guide
overlaps something the virus can't afford to mutate, like a regulatory
element (also not modeled here). Each answers a different question, so
merging them into one number would just hide information.

## 3. Running it

The entry point is `gRNADesign.py`. It always takes two optional positional
arguments — input FASTA, then output TSV path — defaulting to whatever is
configured in `config.py` if you don't supply them:

```bash
python gRNADesign.py                       # uses config.py's defaults
python gRNADesign.py my_panel.fasta results.tsv
```

That first line is shorthand for spelling every one of those same defaults
out explicitly — useful as a template if you want to start from the current
settings and only change one or two flags:

```bash
python gRNADesign.py input/hiv-db.fasta output/gRNA_table.tsv \
  --mm-pickle data/mismatch_score.pkl \
  --pam-pickle data/pam_scores.pkl \
  --g-max 0.10 \
  --p-pam-min 0.80 \
  --plots output/plots
```

Every run prints its progress, ending with something like:

```
Total candidates: 31
No --genome or --offinder-results given; skipping off-target scoring.
Wrote 31 rows to output/gRNA_table.tsv
```

Below is every flag, explained individually:

- **`fasta`** *(positional, optional)* — path to the aligned input FASTA.
  Defaults to `config.py`'s `fasta_path`.
- **`out`** *(positional, optional)* — path to write the output TSV.
  Defaults to `config.py`'s `table_path`.
- **`--mm-pickle`** — path to the CFD mismatch-score pickle. You'd only
  override this if you have a different/updated CFD table than the bundled
  one. Defaults to `config.py`'s `mm_pickle_path`.
- **`--pam-pickle`** — same idea, for the CFD PAM-score pickle. Defaults to
  `config.py`'s `pam_pickle_path`.
- **`--g-max`** — the gap gate threshold (default `0.10`). Lower = stricter:
  a smaller fraction of the panel is allowed to have a gap anywhere in the
  23-column window before it's rejected. Raise this if you're getting very
  few candidates and are okay with somewhat gappier regions.
- **`--p-pam-min`** — the PAM gate threshold (default `0.80`). Higher =
  stricter: a larger fraction of the panel must actually carry `GG` at the
  PAM position. Lower this if you want to consider PAMs that are only
  conserved in, say, 60% of the panel.
- **`--genome`** — path to a reference genome (e.g. an unmasked GRCh38
  directory) to search for off-targets. Supplying this is what turns the
  off-target step on at all; without it, off-target scoring is skipped
  entirely and those two output columns stay blank. See the dedicated
  walkthrough below for exactly what happens when you use this.
- **`--offinder-query`** — where to write the Cas-OFFinder query file that
  `--genome` generates. Defaults to `<out>.offinder_in.txt`. You'd only set
  this if you want a specific filename/location.
- **`--offinder-results`** — path to a Cas-OFFinder output file that was
  already generated (by this tool or by hand, possibly on a different
  machine). Supplying this skips the search step entirely and jumps
  straight to scoring those results — it's the "I already searched, just
  tell me what it means" flag, as opposed to `--genome`'s "search for me."
  If both `--genome` and `--offinder-results` are given, `--offinder-results`
  wins and `--genome` is ignored.
- **`--offinder-mm`** — how many mismatches Cas-OFFinder is allowed to find
  (default `4`). Higher = broader/slower search, catching more distant
  off-targets; the default is already generous because the LTR resembles
  some human endogenous-retrovirus-derived repeats, so a stingy mismatch
  budget risks silently missing real off-targets.
- **`--plots`** *(optional value)* — off by default (no plots generated).
  Passing `--plots` with no value turns plotting on and writes the three
  diagnostic figures described at the end of this document to
  `config.py`'s `plots_dir` (`output/plots` by default). Passing
  `--plots OUT_DIR` turns plotting on and writes them to `OUT_DIR` instead.

**The off-target step, walked through concretely.** Off-target scoring is
entirely optional and only touches disk/runs anything if you ask for it via
one of the two flags above.

If you pass `--genome`:
```bash
python gRNADesign.py my_panel.fasta results.tsv --genome /path/to/GRCh38_dir
```
this writes a Cas-OFFinder query file, then checks whether the `cas-offinder`
program is installed (found on your system `PATH`). If it is, the tool runs
it automatically and merges the results straight into the table. If it
isn't installed, nothing is searched — the tool just prints the exact
command for you to run yourself later:
```
cas-offinder results.tsv.offinder_in.txt C <cas_offinder_out>
```
Once you have that output file (from this run, or generated separately —
possibly on a different, more capable machine), bring the scores in without
repeating the search:
```bash
python gRNADesign.py my_panel.fasta results.tsv --offinder-results cas_offinder_out.txt
```

If neither flag is given, the whole off-target stage is skipped and
`Offtarget_Count` / `Offtarget_Repeat_Flag` are simply left blank — this is
the default, and is intentional; there's no reason to require a multi-
gigabyte human genome file just to look at the on-target conservation table.

All of these same defaults live in one place, `config.py`, if you'd rather
change them permanently than type flags every run.

## 4. Output columns, and the math behind each one

### 4.1 The table itself

| Column | Meaning |
|---|---|
| `gRNA_5to3` | The guide, written as `SPACER(PAM)` — e.g. `GACAAGATATCCTTGATCTG(TGG)`. The 20 letters before the parenthesis are the actual guide sequence (write it as RNA — U instead of T — to synthesize it); the 3 letters in parentheses are the adjacent PAM, shown for context only, not part of the guide itself. |
| `Start_Pos` | The smallest plus-strand alignment column touched by this guide's 23-column footprint (spacer+PAM together) — one shared ruler for both strands. Section 6 works through exactly what this points at for `+` vs `-` guides. |
| `Derived_Strand` | `+` if the spacer was read directly off the plus strand; `-` if found by reverse-complementing first. This describes *how the guide was found*, not which physical strand the guide RNA binds — see section 6. |
| `PAM_Conservation` | Fraction of the 365 panel variants that actually carry `GG` at this window's PAM. Always ≥ `p_pam_min` by construction. |
| `Weighted_Entropy` | How conserved the spacer is, weighted so that variation near the PAM (where Cas9 is least forgiving) counts more than variation at the far end. **Lower = more conserved where it matters.** Range $[0,1]$. |
| `CFD_Score_Independent` | Predicted cutting activity against the panel, assuming the 20 spacer positions vary independently of each other. Range $[0,1]$; $1.0$ = perfect predicted activity everywhere. |
| `CFD_Score_Exact` | The same idea, computed exactly instead of assuming independence — the real CFD formula run against each of the 365 variants' actual sequence, then averaged. A big gap from `CFD_Score_Independent` means certain mismatches travel together across variants rather than occurring independently (a sub-lineage signal). |
| `Offtarget_Count` | Number of human-genome hits scoring above the CFD off-target threshold (default $0.2$). Blank unless off-target scoring ran. |
| `Offtarget_Repeat_Flag` | `True` if this guide has an unusually large number of total genomic matches (default: more than 50), usually meaning it also matches a repetitive element. Blank unless off-target scoring ran. |

Nothing is pre-sorted — in general you want low `Weighted_Entropy`, high
`CFD_Score_Exact`, and (once available) low `Offtarget_Count`.

### 4.2 Notation used below

| Symbol | Meaning |
|---|---|
| $L$ | Alignment length: $866$ columns. |
| $N$ | Number of panel variants: $365$. |
| $i = 1, \dots, 20$ | Spacer position. $i=1$ is the PAM-distal end (the guide's 5′ end); $i=20$ is PAM-proximal (the "seed", right next to the PAM, empirically the least mismatch-tolerant region). |
| $g_i$ | The guide's own (consensus) base at position $i$. |
| $\text{freq}_i(b)$ | Fraction of the $N$ panel variants with base $b \in \{A,C,G,T\}$ at position $i$'s alignment column (ambiguity codes excluded from the denominator). |
| $H_i$ | Shannon entropy of the alignment column at position $i$ (over all 5 symbols $\{A,C,G,T,-\}$). |
| $w_i$ | The CFD-derived weight for position $i$ (higher = a mismatch here hurts cutting more). |
| $\text{act}(i,x,y)$ | CFD single-mismatch activity: fraction of cutting activity that survives if the guide expects base $x$ at position $i$ but the target actually has $y$ there. $\text{act}(i,x,x) = 1$. |

### 4.3 The equations

**Per-column Shannon entropy** (`columns.entropy`) — using log base 5 means
it's already scaled to $[0,1]$ without a separate normalization step, and
treating gap `-` as a 5th symbol means an all-gap column reads the same as
an all-agreeing column (both "not variable"):

$$
H_i = -\sum_{s \in \{A,C,G,T,-\}} p_i(s)\, \log_5 p_i(s), \qquad 0 \log_5 0 := 0
$$

**Position weights** (`cfd.cfd_position_weights`) — averaged over all 12
possible ordered mismatch pairs among $\{A,C,G,T\}$, then flipped so a
position where mismatches hurt *more* gets a *higher* weight. This is
computed once from the CFD table itself, not hand-tuned:

$$
w_i = 1 - \frac{1}{12}\sum_{\substack{x,y \in \{A,C,G,T\}\\ x \neq y}} \text{act}(i,x,y)
$$

**`Weighted_Entropy`** (`score.weighted_entropy`) — the seed-weighted mean
of the 20 spacer columns' entropy:

$$
\text{Weighted\_Entropy} = \frac{\displaystyle\sum_{i=1}^{20} w_i H_i}{\displaystyle\sum_{i=1}^{20} w_i}
$$

**`CFD_Score_Independent`** (`score.expected_cfd`, first return value) — the
expected activity at each position on its own, multiplied across all 20 as
if they varied independently:

$$
a_i = \sum_{b \in \{A,C,G,T\}} \text{freq}_i(b)\cdot \text{act}(i, g_i, b), \qquad \text{act}(i, g_i, {-}) := 0
$$

$$
\text{CFD\_Score\_Independent} = \prod_{i=1}^{20} a_i
$$

**`CFD_Score_Exact`** (`score.expected_cfd`, second return value) — skip the
independence assumption; average the real CFD score against each variant's
actual window sequence:

$$
\text{CFD\_Score\_Exact} = \frac{1}{N}\sum_{v=1}^{N} \text{CFD}(g, s_v, p_v)
$$

where $s_v, p_v$ are variant $v$'s actual spacer and PAM bases in this
window (a variant contributes $0$ if it has a gap/ambiguity anywhere in the
window, or a PAM that isn't in the PAM table).

**The atomic CFD score** (`cfd.calc_cfd`) — the real Doench 2016 / CRISPOR
formula everything above is built on. $wt$ is the guide's intended 20-mer,
$sg$ is whatever it's being compared against (a variant, or a genomic
off-target hit), $pam$ is the 2-letter key for the PAM actually found
(section 7 has the full table), and $\overline{sg_i}$ denotes the
Watson–Crick complement of $sg_i$:

$$
\text{CFD}(wt, sg, pam) = \text{pam\_score}(pam) \cdot \prod_{i\,:\, wt_i \neq sg_i} \text{mismatch\_score}\big(wt_i, \overline{sg_i}, i\big)
$$

If $wt = sg$ exactly, the product is empty (equal to 1), so
$\text{CFD}(wt, wt, pam) = \text{pam\_score}(pam)$ — this identity is
exactly what the CFD validation test checks (section 9). $\text{act}(i,x,y)$
above is just this same formula evaluated on two 20-mers identical
everywhere except position $i$, with $pam$ fixed to `'GG'` and the result
divided back out by $\text{pam\_score}(\text{'GG'})$ to isolate the pure
per-position mismatch effect from the PAM term.

## 5. Why the guide is shown as `SPACER(PAM)`

`gRNA_5to3` shows the PAM in parentheses purely so it's visible at a glance,
without having to look up `Start_Pos` and count 20 characters over — e.g.
`GACAAGATATCCTTGATCTG(TGG)` immediately tells you the PAM here is `TGG`
(the consensus; individual variants may differ slightly, which is exactly
what `PAM_Conservation` quantifies). It's a display-only convenience: off-
target search and CFD scoring both use the bare 20-nt spacer internally, and
the parenthesized PAM is never treated as part of the guide RNA sequence
itself.

## 6. Coordinate system — what `Start_Pos` actually points at

Both strands are reported on **one shared ruler**: `Start_Pos` is always the
*smallest* plus-strand alignment column touched by the guide's full
23-column footprint (spacer + PAM together). This is what makes `+` and `-`
guides directly comparable and sortable by position — but it means
`Start_Pos` points at a *different part of the guide* depending on strand,
which is easiest to see with real numbers from the table.

**`+` strand example** — `GACAAGATATCCTTGATCTG(TGG)`, `Start_Pos = 69`: this
guide was read directly off the plus strand, so its own 5′ end (the first
spacer letter) sits at plus-column 69, and reading left-to-right from there
gives you the spacer, then the PAM, in the same order the guide itself
reads. `Start_Pos` = the guide's own starting point. This is the intuitive
case.

**`-` strand example** — `CTTTATTGAGGCTTAAGCAG(TGG)`, `Start_Pos = 719`:
this guide was *derived* by reverse-complementing the whole alignment and
running the identical search on the mirrored copy. Mapped back onto the
plus strand, its 23-column footprint spans plus-columns 719–741 — but
because this guide's own 5′→3′ direction runs *right-to-left* along the
plus-strand ruler, its 5′ end is actually at the **right** edge
(plus-column 741), and its PAM sits at the **left** edge (plus-columns
719–721). So for this row, `Start_Pos` (719) marks the *last base of the
PAM*, not the start of the spacer:

| Part of the guide | Its own column (mirrored-alignment frame) | Mapped plus-strand column |
|---|---|---|
| Spacer, 1st base (guide's 5′ end) | 124 | 741 |
| Spacer, 20th base (last base before the PAM) | 143 | 722 |
| PAM, 1st base ("N") | 144 | 721 |
| PAM, 3rd base (last "G") | 146 | **719** ← `Start_Pos` |

**The rule, stated plainly:** for `+` guides, `Start_Pos` is where the guide
itself starts. For `-` guides, `Start_Pos` is where the guide's PAM ends,
because that guide reads "backwards" relative to the plus-strand ruler.
Either way, the single consistent rule is "smallest plus-strand column the
23-nt footprint touches" — that's what makes the ruler shared across both
strands; it just isn't always "the start of the spacer" once you're reading
a minus-strand guide.

## 7. CFD reference tables — the actual numbers

Both tables below are the literal contents of `mismatch_score.pkl` and
`pam_scores.pkl` — every score this tool reports is ultimately built from
these numbers, so they're reproduced here in full rather than just
described.

### 7.1 PAM table — all 16 entries

SpCas9's PAM is `NGG`: the first letter can be anything, so only the other
two positions actually get scored. Read this table as: *if the PAM actually
found is the row/column combination below instead of the ideal `GG`, this
fraction of cutting activity survives* (`1.0` = as good as a perfect PAM;
`0.0` = Cas9 essentially can't cut here).

| 1st PAM letter ↓ \ 2nd PAM letter → | A | C | G | T |
|---|---|---|---|---|
| **A** | 0.000 | 0.000 | 0.259 | 0.000 |
| **C** | 0.000 | 0.000 | 0.107 | 0.000 |
| **G** | 0.069 | 0.022 | **1.000** | 0.016 |
| **T** | 0.000 | 0.000 | 0.039 | 0.000 |

(e.g. row `G`, column `G` = the key `"GG"` = activity `1.000` — the ideal
case. Row `A`, column `G` = the key `"AG"` = activity `0.259`.) Only `GG`
gets full activity; almost every other combination is at or near `0` — this
is exactly why the PAM gate in section 2 exists.

### 7.2 Mismatch table — 240 entries, shown as one matrix

For every one of the 20 spacer positions and every one of the 12 possible
single-base mismatches, this is how much cutting activity survives. Rows =
spacer position (1 = guide's 5′/PAM-distal end, 20 = PAM-proximal seed).
Columns = **guide expects → target actually has**. `1.0` = this mismatch,
here, doesn't hurt at all; near `0.0` = this single mismatch, at this single
position, is enough to essentially kill cutting.

| Pos ↓ \ Mismatch → | A→C | A→G | A→T | C→A | C→G | C→T | G→A | G→C | G→T | T→A | T→C | T→G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.857 | 1.000 | 1.000 | 1.000 | 0.913 | 1.000 | 0.900 | 0.714 | 1.000 | 1.000 | 0.857 | 0.957 |
| 2 | 0.786 | 0.800 | 0.727 | 0.727 | 0.696 | 0.909 | 0.846 | 0.692 | 0.636 | 0.846 | 0.857 | 0.840 |
| 3 | 0.429 | 0.611 | 0.706 | 0.867 | 0.500 | 0.688 | 0.750 | 0.385 | 0.500 | 0.714 | 0.429 | 0.500 |
| 4 | 0.353 | 0.625 | 0.636 | 0.842 | 0.500 | 0.800 | 0.900 | 0.529 | 0.364 | 0.476 | 0.647 | 0.625 |
| 5 | 0.500 | 0.720 | 0.364 | 0.571 | 0.600 | 0.636 | 0.867 | 0.786 | 0.300 | 0.500 | 1.000 | 0.640 |
| 6 | 0.455 | 0.714 | 0.714 | 0.929 | 0.500 | 0.929 | 1.000 | 0.682 | 0.667 | 0.867 | 0.909 | 0.571 |
| 7 | 0.438 | 0.706 | 0.438 | 0.750 | 0.471 | 0.813 | 1.000 | 0.688 | 0.571 | 0.875 | 0.688 | 0.588 |
| 8 | 0.429 | 0.733 | 0.429 | 0.650 | 0.643 | 0.875 | 1.000 | 0.615 | 0.625 | 0.800 | 1.000 | 0.733 |
| 9 | 0.571 | 0.667 | 0.600 | 0.857 | 0.619 | 0.875 | 0.643 | 0.538 | 0.533 | 0.929 | 0.923 | 0.619 |
| 10 | 0.333 | 0.556 | 0.882 | 0.867 | 0.389 | 0.941 | 0.933 | 0.400 | 0.813 | 0.857 | 0.533 | 0.500 |
| 11 | 0.400 | 0.650 | 0.308 | 0.750 | 0.250 | 0.308 | 1.000 | 0.429 | 0.385 | 0.750 | 0.667 | 0.400 |
| 12 | 0.263 | 0.722 | 0.333 | 0.714 | 0.444 | 0.538 | 0.933 | 0.529 | 0.385 | 0.800 | 0.947 | 0.500 |
| 13 | 0.211 | 0.652 | 0.300 | 0.385 | 0.136 | 0.700 | 0.923 | 0.421 | 0.300 | 0.692 | 0.789 | 0.261 |
| 14 | 0.214 | 0.467 | 0.533 | 0.350 | 0.000 | 0.733 | 0.750 | 0.429 | 0.267 | 0.619 | 0.286 | 0.000 |
| 15 | 0.273 | 0.650 | 0.200 | 0.222 | 0.050 | 0.067 | 0.941 | 0.273 | 0.143 | 0.579 | 0.273 | 0.050 |
| 16 | 0.000 | 0.192 | 0.000 | 1.000 | 0.154 | 0.308 | 1.000 | 0.000 | 0.000 | 0.909 | 0.667 | 0.346 |
| 17 | 0.176 | 0.176 | 0.133 | 0.467 | 0.059 | 0.467 | 0.933 | 0.235 | 0.250 | 0.533 | 0.706 | 0.118 |
| 18 | 0.190 | 0.400 | 0.500 | 0.538 | 0.133 | 0.643 | 0.692 | 0.476 | 0.667 | 0.667 | 0.429 | 0.333 |
| 19 | 0.207 | 0.375 | 0.538 | 0.429 | 0.125 | 0.462 | 0.714 | 0.448 | 0.667 | 0.286 | 0.276 | 0.250 |
| 20 | 0.227 | 0.765 | 0.600 | 0.500 | 0.059 | 0.300 | 0.938 | 0.429 | 0.700 | 0.563 | 0.091 | 0.176 |

**How to read a single cell:** row 20, column `A→G` = `0.765` means: if the
guide is designed expecting `A` right next to the PAM (position 20, the
seed) but a panel variant actually has `G` there, about 76.5% of cutting
activity survives that particular mismatch. Compare that to row 20, column
`C→G` = `0.059` — almost total loss of activity, for a *different*
mismatch at the *exact same position*. This is why the tool never treats
"number of mismatches" as a single penalty: which base, mismatched to
which other base, at which position, all matter — and this table (measured
data, not a guess) is where every `CFD_Score_*` value ultimately comes
from. You can regenerate this exact matrix yourself by calling
`cfd.mm_activity(pos0, x, y, mm_scores, pam_scores)` for every position and
base pair.

## 8. Module and function reference

Every file in the project, and what each function in it actually does.

### `config.py`
- **`CONFIG`** (dict) — the single source of truth for every tunable value
  in the project: input file paths (`fasta_path`, `mm_pickle_path`,
  `pam_pickle_path`), the enzyme's parameters (`pam`, `spacer_len`,
  `pam_len` — all fixed for SpCas9), the two window gates (`g_max`,
  `p_pam_min`), off-target parameters (`genome_path`, `offinder_mm`,
  `offtarget_cfd_threshold`, `repeat_hit_threshold`), and output locations
  (`table_path`, `plots_dir`). Nothing else in the project hardcodes these
  values — change them here to change them everywhere.

### `io_fasta.py`
- **`read_alignment(path)`** — reads the aligned FASTA using Biopython's
  `SeqIO`, uppercases every sequence, checks that all rows are the same
  length (raises `ValueError` if not — a real alignment must be rectangular),
  and returns a plain list of strings.
- **`reverse_complement_alignment(rows)`** — returns a new list where every
  row has been reverse-complemented via Biopython's
  `Seq.reverse_complement()`, which correctly maps gaps to gaps and IUPAC
  ambiguity codes to their complementary code (e.g. `M`↔`K`, `N`↔`N`). This
  is what lets the minus strand be searched with identical logic to the
  plus strand.

### `cfd.py`
- **`load_cfd(mm_path, pam_path)`** — unpickles both files, then identifies
  which is the mismatch table and which is the PAM table by inspecting a
  sample key's shape (mismatch keys look like `"rG:dA,14"`; PAM keys are
  bare 2-letter strings like `"AG"`) rather than trusting the filenames.
  Automatically swaps and prints a warning if they were passed in the wrong
  order. Returns `(mm_scores, pam_scores)`.
- **`calc_cfd(wt, sg, pam, mm_scores, pam_scores)`** — the atomic CFD score
  (equation in section 4.3): compares two 20-mers position by position,
  looks up a mismatch penalty for every position where they differ,
  multiplies all of those together with the PAM's own score, and returns
  the product.
- **`mm_activity(pos0, guide_base, target_base, mm_scores, pam_scores)`** —
  the activity for a single mismatch at one position: returns `1.0`
  immediately if the bases match or `0.0` if either base is a gap;
  otherwise builds two 20-mers that are identical everywhere except
  `pos0`, calls `calc_cfd` on them with `pam='GG'`, and divides out
  `pam_scores['GG']` so only the position-specific mismatch penalty
  remains.
- **`cfd_position_weights(mm_scores, pam_scores)`** — computes the weight
  vector $w_i$ from section 4.3, by averaging `mm_activity` over all 12
  possible ordered mismatch pairs at each of the 20 positions.

### `columns.py`
- **`column_frequencies(rows)`** — vectorized with numpy: stacks all rows
  into a character array, counts `A`/`C`/`G`/`T`/`-` per column (ambiguity
  codes are automatically excluded, since they match none of those five
  masks), and normalizes each column's counts into frequencies. Returns
  `(freqs, gapfrac)` — a list of per-column frequency dicts, and a
  parallel list of just the gap fractions.
- **`consensus(freqs)`** — for each column, returns whichever of `A/C/G/T`
  has the highest frequency (ties broken deterministically by `A,C,G,T`
  order); gaps and ambiguity codes are never returned as a consensus base.
- **`entropy(freqs)`** — the Shannon entropy formula from section 4.3,
  vectorized with numpy, returning one float per column.

### `windows.py`
- **`enumerate_windows(consensus, freqs, gapfrac, rows, strand, g_max, p_pam_min, L)`**
  — slides a 23-column window across the given alignment, applies the gap
  gate then the PAM gate (section 2, step 4), and for every window that
  survives both, builds a `Window` record: the spacer sequence, the PAM
  consensus, which alignment columns each occupies, which strand it came
  from, its plus-strand `start_col` (via `coords.plus_coord`), and its PAM
  conservation fraction.

### `coords.py`
- **`plus_coord(start_col, strand, L)`** — the single function responsible
  for the shared coordinate ruler described in section 6. For `+` strand,
  returns the column unchanged. For `-` strand, maps a window found at
  column `j` in the mirrored alignment to plus-strand column
  $L - (j + 23)$.

### `score.py`
- **`weighted_entropy(window, H, w)`** — implements the `Weighted_Entropy`
  equation from section 4.3 for one window.
- **`expected_cfd(window, rows, freqs, mm_scores, pam_scores)`** —
  implements both `CFD_Score_Independent` and `CFD_Score_Exact` (section
  4.3) for one window, returning both as a tuple.

### `table.py`
- **`assemble_table(windows, H, w, rows, freqs, mm_scores, pam_scores)`** —
  scores every window (via `score.py`) and builds the final row dicts,
  including the bare-spacer internal key (`_spacer`, used for off-target
  matching) and the display string (`gRNA_5to3`, `SPACER(PAM)`).
- **`merge_offtarget(out_rows, offtarget_summary)`** — fills in
  `Offtarget_Count` / `Offtarget_Repeat_Flag` on each row by matching on its
  `_spacer` key against `offtarget.score_offinder`'s results.
- **`to_dataframe(out_rows)`** / **`write_table(out_rows, path)`** — build a
  pandas `DataFrame` restricted to the public columns (internal keys like
  `_spacer` never make it into the file) and write it out as a TSV.

### `offtarget.py`
- **`emit_offinder(guides, genome_path, out_path, mm)`** — writes a
  Cas-OFFinder query file: the genome path, the search pattern
  (`N`×20 + `NGG`), and one query line per guide (`"<spacer>NGG <mm>"`).
- **`score_offinder(out_file, guides, mm_scores, pam_scores, cfd_threshold, repeat_hit_threshold)`**
  — parses Cas-OFFinder's raw hit list, matches each hit back to the guide
  that produced it, CFD-scores every hit (via `calc_cfd`) instead of just
  counting raw mismatches, and returns a per-guide summary: total hit
  count, how many scored above `cfd_threshold` (this becomes
  `Offtarget_Count`), and whether the guide has more than
  `repeat_hit_threshold` total hits (`n_repeat_hits`, becomes
  `Offtarget_Repeat_Flag`).

### `diagnostics.py`
*(opt-in — meant to be called by hand on a shortlist, not run automatically
on every candidate)*
- **`weighted_entropy_profile(window, H, w)`** — returns the un-aggregated,
  length-20 list of $w_i H_i$ values for one window, so you can see exactly
  *which* position(s) are driving that guide's `Weighted_Entropy` score
  instead of only the final aggregate number.
- **`per_variant_distance(window, rows, mm_scores, pam_scores)`** — for
  each of the 365 panel variants, sums $1 - \text{act}(i, g_i, \cdot)$
  across the 20 spacer positions (a gap/ambiguity counts as maximal
  distance, `1.0`), giving one CFD-weighted "distance from this guide's
  consensus" number per variant.
- **`summarize_distances(distances)`** — median and 90th percentile of that
  distance list, a quick way to see whether one guide's distances cluster
  tightly or split into sub-groups (possible sub-lineages).

### `plotting.py`
*(all optional — see section 10)*
- **`compute_profiles(fasta_path)`** — loads the alignment, builds its
  reverse complement, and returns per-column entropy for both
  (`H_plus`, `H_minus`), per-column gap fraction for both (`gapfrac_plus`,
  `gapfrac_minus`), and the alignment length `L`.
- **`rolling_window_stat(H, weights, n_starts)`** — a generic sliding
  dot-product helper: for each start position, multiplies a 20-value
  window by a weight vector and sums it (used for the unweighted mean, the
  CFD-weighted mean, and the gap-fraction mean curves alike).
- **`windowed_profiles(H_plus, H_minus, L, w)`** — builds the four entropy
  background curves (plus/minus strand × unweighted/weighted) already
  mapped onto the shared `Start_Pos` ruler from section 6.
- **`gap_windowed_profile(gapfrac_plus, gapfrac_minus, L)`** — the same
  idea, for gap fraction instead of entropy: builds the plus/minus mean
  gap-fraction curves on the same shared ruler.
- **`plot_entropy_by_position`, `plot_window20_unweighted`, `plot_window20_weighted`, `plot_gap_fraction_window20`**
  — draw the four figures described in section 10, including the
  candidate-point overlay (all but the first).
- **`generate_all_plots(out_dir, fasta_path, mm_pickle_path, pam_pickle_path, candidates)`**
  — orchestrates all of the above and saves the four PNGs to `out_dir`.
- **`main()`** — CLI entry point for regenerating plots standalone from an
  existing `table.tsv`, without rerunning the whole pipeline.

### `gRNADesign.py`
- **`run_pipeline(fasta_path, mm_path, pam_path, g_max, p_pam_min)`** — the
  core pipeline, steps 1–6 of section 2, returning the final list of row
  dicts.
- **`run_offtarget(out_rows, args)`** — the optional off-target
  orchestration described in section 3 (`--genome` / `--offinder-results`),
  mutating `out_rows` in place.
- **`main(argv)`** — parses every flag from section 3 and wires everything
  together: runs the pipeline, optionally runs off-target scoring, writes
  the table, and optionally generates plots.

## 9. Tests — what they check, and why

```bash
python -m unittest discover -s tests -t .
```

There are three tiers, each checking a different kind of thing that could
otherwise go silently wrong:

- **`tests/test_unit.py`** — small, fast checks on individual functions
  using hand-built toy data (not the real HIV panel), so each check
  isolates one specific behavior:
  - reverse-complementing a sequence twice returns the original (an
    "involution" — if this ever fails, the minus-strand search is silently
    broken);
  - a column where all 5 symbols (`A,C,G,T,-`) appear equally often has
    entropy exactly `1.0`, and a column that's the same base in every
    sequence has entropy exactly `0.0` (the two edge cases the entropy
    formula must get exactly right);
  - an ambiguity code (`N`, etc.) doesn't get counted as any base when
    computing frequencies;
  - the gap gate and PAM gate each correctly keep/reject a window built
    specifically to sit right on either side of the threshold.
- **`tests/test_cfd_validation.py`** — checks specific to the CFD table,
  and arguably the single most important test file in the project. CFD
  data has a well-known gotcha: different implementations index the 20
  spacer positions in opposite directions, and getting this backwards
  doesn't crash anything — it just silently produces plausible-looking but
  wrong scores everywhere. This file checks:
  - scoring a guide against itself (no mismatches) returns exactly
    `pam_scores['GG']` (the identity mentioned in section 4.3);
  - `pam_scores['GG']` is approximately `1.0` (a perfect PAM should cost
    essentially nothing);
  - `cfd_position_weights()` is heavier toward the PAM-proximal end
    (positions 13–20 average higher than positions 1–8) — if this comes
    out reversed, the position indexing is backwards and every score in
    the table is wrong, even though nothing would have crashed.
- **`tests/test_integration.py`** — runs the *real* pipeline against the
  actual FASTA panel and CFD tables (not toy data), and checks the output
  makes sense as a whole: at least one candidate survives; every
  candidate's `PAM_Conservation` really is above the configured threshold;
  every guide's bare spacer is a clean 20-letter `A/C/G/T` string with no
  leftover gaps or ambiguity codes; the `gRNA_5to3` display string is
  correctly formatted as `SPACER(PAM)`; and running the whole pipeline
  twice produces byte-for-byte identical output (no hidden randomness or
  ordering nondeterminism anywhere in the pipeline).

If you ever change any scoring logic, running this suite is the fast way to
find out whether you broke something — especially the CFD validation tests,
since a CFD indexing bug is exactly the kind of mistake that won't announce
itself any other way.

## 10. Diagnostic plots (`plotting.py`)

Entirely optional. Generate them either as part of a normal run:

```bash
python gRNADesign.py my_panel.fasta results.tsv --plots plots/
```

or standalone, from an existing table, without rerunning the pipeline:

```bash
python plotting.py --table results.tsv --out-dir plots/
```

Four figures are produced, all sharing the same `Start_Pos` ruler as the
main table, with a dashed gridline every 100 columns to make positions easy
to read off:

1. **`entropy_by_position.png`** — raw per-column Shannon entropy across the
   whole alignment, both strands. (The two strands' lines coincide almost
   exactly — entropy doesn't change under complementation, so this is
   expected, not a bug.)
2. **`entropy_window20_unweighted.png`** — plain mean entropy over every
   possible 20-column spacer window (a continuous version of the
   conservation signal, computed everywhere along the genome, not just at
   the surviving gated candidates).
3. **`entropy_window20_weighted.png`** — CFD-seed-weighted mean entropy
   over every 20-column window — this is exactly `Weighted_Entropy`,
   computed continuously instead of only at the candidate windows that
   passed the gates.
4. **`gap_fraction_window20.png`** — mean gap (`-`) fraction over every
   20-column spacer window, i.e. what fraction of the 365 panel variants
   have an actual deletion there (gap count / total sequences). This is the
   same quantity the gap gate (`g_max`, section 2) filters on, computed
   continuously along the whole genome instead of only as a pass/fail
   check — the tall spikes up to `1.0` are regions where most of the panel
   has a deletion (commonly variable-length loop regions), and every
   surviving candidate correctly clusters near `0.0`, since the gap gate
   only lets through windows with almost no gaps.

On plots 2, 3, and 4, every surviving candidate is overlaid as a point:
marker **shape** encodes `Derived_Strand` (circle = `+`, triangle = `-`),
and marker **color** encodes `CFD_Score_Exact` (see the color bar) — so you
can see at a glance both where a guide sits in the conservation/gap
landscape and how good its predicted cutting activity is.

## 11. Credits

The CFD (Cutting Frequency Determination) lookup tables bundled here
(`mismatch_score.pkl`, `pam_scores.pkl`) are not this project's own data —
they are the published, experimentally-measured scores from:

> Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan, K. F., Smith, I., Tothova, Z., Wilen, C., Orchard, R., Virgin, H. W., Listgarten, J., & Root, D. E. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature biotechnology, 34(2), 184–191. https://doi.org/10.1038/nbt.3437

These are the same numbers used by [CRISPOR](https://github.com/maximilianh/crisporWebsite/tree/master/CFD_Scoring), and
are reproduced here in full in section 7 for transparency. All credit for
the underlying measurements belongs to Doench et al.; this project only
consumes the tables (`cfd.py`) to score conservation across the HIV-1 panel.
