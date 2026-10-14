"""The candidate map.

One figure, four panel sharing the reference coordinate as a common x axis:

    1  CFD score per candidate, as a stem whose marker area carries the MFE
    2  running window mean Shannon entropy, as a line
    3  the reference annotation
    4  the guide footprint, with its direction

Drawn in the Swiss International Typographic manner: one grid, flush left
text, a single sans serif in two weight, hairline rule, and a palette of ink
and paper with one signal accent.
"""

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

import annotation

# --- Palette: ink, paper, one accent ---
INK = "#111111"
GRAPHITE = "#6E6E6E"
RULE = "#D7D4CE"
FILL = "#EFEDE8"
PAPER = "#FFFFFF"

# Strand carry the colour, and the same two colour are used by the CFD stem
# and by the footprint arrow so one reads straight onto the other. Both are
# muted rather than primary: a saturated red against a saturated blue
# vibrates, and neither should outshout the ink the rest of the figure is in.
STRAND_PLUS = "#C14B3E"
STRAND_MINUS = "#2F6690"

# --- Type scale ---
TITLE, SUBTITLE, LABEL, TICK, SMALL = 19.0, 10.0, 9.0, 8.5, 7.8

_MARKER_MIN, _MARKER_MAX = 26.0, 210.0


def _style():
    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "text.color": INK,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.6,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": TICK,
        "ytick.labelsize": TICK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "figure.facecolor": PAPER,
        "axes.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "legend.frameon": False,
    })


def _bare(axis, keep_bottom=False, rule=True):
    """Strip every rule that is not carrying information.

    The bottom rule stays on every panel: it is the baseline the whole figure
    is built on, and a panel floating without one reads as unfinished.
    """
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.spines["left"].set_color(INK)
    axis.spines["bottom"].set_visible(True)
    axis.spines["bottom"].set_color(INK)
    axis.spines["bottom"].set_linewidth(0.9)
    if not keep_bottom:
        axis.tick_params(axis="x", length=0, labelbottom=False)
    if rule:
        axis.set_axisbelow(True)
        axis.grid(axis="y", color=RULE, linewidth=0.5)


_MARKER_ABSENT = 58.0


def _area_of(magnitude, extent):
    """Area for one |MFE|, on the scale fixed by `extent`."""
    low, high = extent
    return _MARKER_MIN + (magnitude - low) / (high - low) * (_MARKER_MAX - _MARKER_MIN)


def _marker_area(energy):
    """Map |MFE| to a marker area. A stronger fold reads as a larger marker.

    Returns (area, extent); extent is None when there is nothing to scale by,
    and the caller then draws no size key.
    """
    value = np.abs(np.asarray(energy, dtype=float))
    if not np.isfinite(value).any():
        return np.full(value.shape, _MARKER_ABSENT), None
    low, high = float(np.nanmin(value)), float(np.nanmax(value))
    if high - low < 1e-9:
        return np.where(np.isfinite(value), 96.0, _MARKER_ABSENT), None
    area = _area_of(value, (low, high))
    return np.where(np.isfinite(area), area, _MARKER_ABSENT), (low, high)


def _pack(item, allowance):
    """Assign each (start, end, payload) to the lowest free row."""
    occupied = []
    placed = []
    for start, end, payload in sorted(item, key=lambda entry: entry[0]):
        width = max(end, start + allowance(payload))
        for index, edge in enumerate(occupied):
            if start > edge:
                occupied[index] = width
                placed.append((index, start, end, payload))
                break
        else:
            occupied.append(width)
            placed.append((len(occupied) - 1, start, end, payload))
    return placed, len(occupied)


def _arrow(axis, start, end, level, height, forward, color, span=None):
    """A sleek arrow spanning [start, end], pointing along the strand.

    A slim body with a proportionally long head, so direction reads at a
    glance even where the footprint is only a few pixel wide. The head is
    capped against the axis span so a long feature does not grow a huge one.
    """
    width = end - start
    head = width * 0.44
    if span:
        head = min(head, span * 0.013)
    body = height * 0.46
    top, bottom = level + height / 2, level - height / 2
    if forward:
        tip, neck = end, end - head
        point = [(start, level + body / 2), (neck, level + body / 2), (neck, top),
                 (tip, level), (neck, bottom), (neck, level - body / 2),
                 (start, level - body / 2)]
    else:
        tip, neck = start, start + head
        point = [(end, level + body / 2), (neck, level + body / 2), (neck, top),
                 (tip, level), (neck, bottom), (neck, level - body / 2),
                 (end, level - body / 2)]
    axis.add_patch(Polygon(point, closed=True, facecolor=color,
                           edgecolor="none", joinstyle="miter", zorder=3))


def draw(
    frame,
    entropy_position,
    entropy_value,
    reference,
    path,
    *,
    window=20,
    cfd_line=0.80,
    title="gRNA candidate",
    subtitle="",
    footnote="",
):
    """Render the figure and write it to `path`.

    frame: the output table, already restricted to the candidate to draw.
    entropy_position / entropy_value: the running mean entropy track, in
        reference coordinate.
    """
    _style()
    plotted = frame.dropna(subset=["Ref_Start", "Ref_End"]).copy()

    low = float(np.nanmin(entropy_position)) if len(entropy_position) else 0.0
    high = float(np.nanmax(entropy_position)) if len(entropy_position) else 1.0
    if len(plotted):
        low = min(low, float(plotted["Ref_Start"].min()))
        high = max(high, float(plotted["Ref_End"].max()))
    pad = max(8.0, (high - low) * 0.02)
    limit = (low - pad, high + pad)

    figure = plt.figure(figsize=(16.0, 9.6), dpi=200)
    grid = figure.add_gridspec(
        4, 1, height_ratios=[2.85, 1.45, 1.9, 1.5], hspace=0.17,
        left=0.062, right=0.982, top=0.845, bottom=0.140,
    )
    cfd_axis = figure.add_subplot(grid[0])
    entropy_axis = figure.add_subplot(grid[1], sharex=cfd_axis)
    gene_axis = figure.add_subplot(grid[2], sharex=cfd_axis)
    guide_axis = figure.add_subplot(grid[3], sharex=cfd_axis)

    _title_block(figure, title, subtitle)
    _draw_cfd(cfd_axis, plotted, cfd_line)
    _draw_entropy(entropy_axis, entropy_position, entropy_value, window)
    _draw_gene(gene_axis, reference, limit)
    _draw_guide(guide_axis, plotted, limit[1] - limit[0])

    guide_axis.set_xlim(*limit)
    guide_axis.set_xlabel(f"{reference.identifier} coordinate, nt",
                          fontsize=LABEL, labelpad=6, loc="left")
    if footnote:
        figure.text(0.062, 0.015, footnote, fontsize=SMALL, color=GRAPHITE,
                    ha="left", va="bottom", linespacing=1.7)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def _title_block(figure, title, subtitle):
    figure.text(0.062, 0.955, title, fontsize=TITLE, fontweight="bold",
                ha="left", va="top", color=INK)
    if subtitle:
        figure.text(0.062, 0.912, subtitle, fontsize=SUBTITLE, ha="left",
                    va="top", color=GRAPHITE, linespacing=1.65)
    figure.add_artist(Line2D([0.062, 0.982], [0.873, 0.873], color=INK,
                             linewidth=0.9, transform=figure.transFigure))


def _draw_cfd(axis, frame, cfd_line):
    _bare(axis)
    axis.set_ylabel("CFD score\nagainst the panel", fontsize=LABEL, linespacing=1.5)

    if not len(frame):
        axis.text(0.5, 0.5, "no anchored candidate", transform=axis.transAxes,
                  ha="center", color=GRAPHITE, fontsize=LABEL)
        return

    centre = (frame["Ref_Start"] + frame["Ref_End"]) / 2.0
    value = frame["CFD_Score"].to_numpy(dtype=float)
    area, extent = _marker_area(frame["MFE"].to_numpy(dtype=float))
    forward = (frame["Derived_Strand"] == "+").to_numpy()
    colour = np.where(forward, STRAND_PLUS, STRAND_MINUS)

    axis.set_ylim(0.0, max(1.04, float(np.nanmax(value)) * 1.1))
    axis.vlines(centre, 0.0, value, color=colour, linewidth=0.9, zorder=2)
    axis.scatter(centre, value, s=area, c=colour, edgecolor=PAPER,
                 linewidth=0.7, zorder=3)

    if cfd_line is not None:
        axis.axhline(cfd_line, color=GRAPHITE, linewidth=0.7,
                     linestyle=(0, (4, 3)), zorder=1)
        axis.text(0.0, cfd_line, f"  CFD {cfd_line:.2f}", fontsize=SMALL,
                  color=GRAPHITE, ha="left", va="bottom",
                  transform=axis.get_yaxis_transform())

    handle = [
        Line2D([], [], marker="o", linestyle="none", markersize=7,
               markerfacecolor=STRAND_PLUS, markeredgecolor=PAPER,
               label=f"protospacer on (+)   n = {int(forward.sum())}"),
        Line2D([], [], marker="o", linestyle="none", markersize=7,
               markerfacecolor=STRAND_MINUS, markeredgecolor=PAPER,
               label=f"protospacer on (−)   n = {int((~forward).sum())}"),
    ]
    axis.add_artist(axis.legend(
        handles=handle, loc="upper left", fontsize=SMALL, labelspacing=0.55,
        handletextpad=0.6, borderaxespad=0.6, alignment="left"))

    if extent is not None:
        low, high = extent
        size = [
            Line2D([], [], marker="o", linestyle="none",
                   markersize=np.sqrt(_area_of(magnitude, extent)),
                   markerfacecolor="none", markeredgecolor=INK,
                   markeredgewidth=0.7,
                   label=f"{-magnitude if magnitude else 0.0:.1f}")
            for magnitude in (low, (low + high) / 2.0, high)
        ]
        axis.legend(handles=size, loc="upper right", fontsize=SMALL,
                    title="MFE, kcal/mol\n(marker area)",
                    title_fontsize=SMALL, labelspacing=1.45,
                    handletextpad=1.0, borderaxespad=0.6, alignment="left")


def _draw_entropy(axis, position, value, window):
    _bare(axis)
    axis.set_ylabel(f"Shannon entropy H\n{window} bp running mean",
                    fontsize=LABEL, linespacing=1.5)
    if len(position):
        axis.plot(position, value, color=INK, linewidth=1.15, zorder=3)
        axis.set_ylim(0.0, max(0.05, float(np.nanmax(value)) * 1.25))
    for tick in axis.get_yticks():
        axis.axhline(tick, color=RULE, linewidth=0.5, zorder=0)
    axis.text(0.0, 1.03, "lower = more conserved across the panel",
              transform=axis.transAxes, fontsize=SMALL, color=GRAPHITE,
              ha="left", va="bottom")


def _draw_gene(axis, reference, limit):
    """The annotation track.

    A broad feature (a gene, a region) is drawn as a stacked bar at the
    bottom. A short element (a TATA box, a factor site) is too narrow to carry
    its own label, so it is drawn as a tick on a shared line with a leader up
    to a stacked label, which is what keeps a dense promoter readable.
    """
    axis.set_ylabel("reference\nannotation", fontsize=LABEL, linespacing=1.5)
    for side in ("top", "right", "left", "bottom"):
        axis.spines[side].set_visible(False)
    axis.tick_params(axis="both", length=0, labelbottom=False, labelleft=False)

    low, high = limit
    span = high - low
    visible = [
        (max(feature.start + 1, low), min(feature.end, high), feature)
        for feature in reference.feature
        if feature.kind not in annotation.structural()
        and feature.start + 1 <= high and feature.end >= low
    ]
    if not visible:
        axis.set_ylim(0, 1)
        return

    def label_width(feature):
        return len(feature.name) * span * 0.0050

    # A feature is drawn as a labelled bar only when its own label fits inside
    # it. Otherwise the label would overflow and collide with the neighbouring
    # track, so the feature goes to the tick and leader treatment instead.
    wide = [entry for entry in visible
            if entry[1] - entry[0] >= label_width(entry[2]) * 1.15]
    narrow = [entry for entry in visible if entry not in wide]

    wide.sort(key=lambda entry: (entry[0] - entry[1], entry[0]))
    bar, rows = _pack(wide, lambda feature: 0.0) if wide else ([], 0)
    height, step = 0.44, 0.62
    for row, start, end, feature in bar:
        level = row * step
        axis.add_patch(Rectangle(
            (start, level), max(end - start, span * 0.0012), height,
            facecolor=FILL, edgecolor=INK, linewidth=0.6, zorder=2))
        axis.text(start + span * 0.004, level + height / 2, feature.name,
                  fontsize=SMALL, color=INK, ha="left", va="center", zorder=3)

    top = rows * step
    if not narrow:
        axis.set_ylim(-0.12, max(top, 1.0))
        return

    # The label is centred over its tick, so what has to be packed is the
    # label's own extent either side of the midpoint, not the feature's span.
    item = []
    for start, end, feature in narrow:
        middle = (start + end) / 2.0
        half = label_width(feature) / 2.0
        item.append((middle - half, middle + half, (start, end, feature)))

    tick, label_rows = _pack(item, lambda payload: 0.0)
    tick_base, label_base, label_step = top + 0.12, top + 0.66, 0.46
    for row, _, _, (start, end, feature) in tick:
        level = label_base + (label_rows - 1 - row) * label_step
        middle = (start + end) / 2.0
        axis.add_patch(Rectangle(
            (start, tick_base), max(end - start, span * 0.0016), 0.26,
            facecolor=INK, edgecolor="none", zorder=3))
        axis.plot([middle, middle], [tick_base + 0.26, level], color=RULE,
                  linewidth=0.6, zorder=1)
        axis.text(middle, level, feature.name, fontsize=SMALL, color=INK,
                  ha="center", va="bottom", zorder=3)
    axis.set_ylim(-0.12, label_base + label_rows * label_step + 0.25)


def _draw_guide(axis, frame, span):
    """The footprint track.

    A candidate panel is full of guide that start a base or two apart, and
    drawn on one line those merge into an unreadable block. Each strand is
    therefore packed into as many sub lane as it takes for no two arrow to
    touch, with a clear gap enforced between neighbour.
    """
    _bare(axis, keep_bottom=True)
    axis.set_ylabel("guide\nfootprint", fontsize=LABEL, linespacing=1.5)
    axis.tick_params(axis="y", length=0)

    gap = span * 0.006
    lane = {}
    for sign, forward in ((1, True), (-1, False)):
        entry = [(float(row["Ref_Start"]), float(row["Ref_End"]), row)
                 for _, row in frame.iterrows()
                 if (row["Derived_Strand"] == "+") is forward]
        lane[forward] = _pack(entry, lambda payload: gap) if entry else ([], 0)

    rows = max(lane[True][1], lane[False][1], 1)
    height = min(0.62, 1.30 / rows)
    step = height * 1.34
    axis.set_ylim(-(rows * step + 0.30), rows * step + 0.30)
    axis.set_yticks([step * rows / 2.0, -step * rows / 2.0])
    axis.set_yticklabels(["+", "−"], fontsize=LABEL)
    axis.axhline(0.0, color=RULE, linewidth=0.6, zorder=1)

    for forward in (True, False):
        placed, _ = lane[forward]
        sign = 1 if forward else -1
        colour = STRAND_PLUS if forward else STRAND_MINUS
        for row, start, end, _payload in placed:
            level = sign * (step * (row + 0.5) + 0.06)
            _arrow(axis, start, end, level, height, forward, colour, span=span)

