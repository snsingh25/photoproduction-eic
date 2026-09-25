#!/usr/bin/env python3
"""Shared publication style for the EIC_Observables paper figures.

Design contract
---------------
Typography   Computer Modern via LaTeX, so figure text is the same face and
             the same optical weight as the REVTeX body text.
Geometry     Square panels. Single-column width 3.375in, full width 6.9in.
Colour       Two palettes, each doing exactly one job:

             observables -> nominal categorical. Three fixed slots, assigned
                            in order, never cycled. Backed by line style so
                            identity never rests on hue alone.
             energies    -> ordinal. One hue, monotone light->dark with
                            sqrt(s), so the reader sees the energy ordering
                            in the colour itself. This matters because the
                            paper's claim is about a trend in sqrt(s).

Both palettes were checked with the dataviz validator against a white
surface: the categorical trio passes all-pairs CVD (worst dE 9.2 deutan) and
the normal-vision floor (worst dE 24.0); the ordinal ramp passes
monotonicity, adjacent dL, and light-end contrast. The aqua slot sits at
2.82:1 on white, below the 3:1 mark floor, so every figure that uses it
carries a visible direct label (the AUC value in the legend) as relief.
"""

import matplotlib as mpl

# --- ink ---------------------------------------------------------------
# Chrome is neutral gray, not the design system's warm chrome: the journal
# page is pure white, and the warm steps read muddy beside the cool series
# hues in print. Series colours are unchanged from the validated palette.
INK = "#0b0b0b"      # primary text
MUTED = "#8a8a8a"    # tick marks, chance line, annotations
GRID = "#e9e9e9"     # hairline grid
AXIS = "#bdbdbd"     # frame

# --- observables: nominal categorical, fixed slot order ----------------
# slot 1 blue, slot 2 orange, slot 3 aqua
OBS_COLOR = {
    "jet_psi03":    "#2a78d6",
    "jet_nsd":      "#eb6834",
    "jet_nsubjets": "#1baf7a",
}
# secondary encoding, so identity is never colour-alone
OBS_DASH = {
    "jet_psi03":    (None, None),
    "jet_nsd":      (4.0, 1.6),
    "jet_nsubjets": (1.2, 1.4),
}
OBS_LABEL = {
    "jet_psi03":    r"$\Psi(r=0.3)$",
    "jet_nsd":      r"$n_{\mathrm{SD}}$",
    "jet_nsubjets": r"$n_{\mathrm{subjets}}$",
}
OBS_ORDER = ["jet_psi03", "jet_nsd", "jet_nsubjets"]

# --- energies: ordinal ramp, one hue, light -> dark with sqrt(s) -------
# blue ramp steps 250 / 400 / 500 / 650; light end clears 2:1 on white
ENERGY_COLOR = {
    64:  "#86b6ef",
    105: "#3987e5",
    141: "#256abf",
    300: "#104281",
}

# --- sequential ramp, for density / magnitude (one hue, light -> dark) --
# Same blue family as the ordinal energy ramp, extended to the full
# 100->700 range. The light end is allowed to recede into the white page,
# which is correct for "near zero" in a density map.
SEQ_BLUE = [
    "#ffffff", "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
    "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95",
    "#104281", "#0d366b",
]


def density_cmap(name="paper_blue"):
    """Sequential white->dark-blue colormap for 2D density plots."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, SEQ_BLUE)


# --- figure geometry ---------------------------------------------------
COL_W = 3.375   # REVTeX single column
FULL_W = 6.9    # REVTeX full width


def use():
    """Install the paper rcParams. Call once at the top of a plot script."""
    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",

        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,

        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.6,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,

        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "grid.alpha": 1.0,
        "axes.axisbelow": True,

        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,

        "legend.frameon": False,
        "legend.handlelength": 1.9,
        "legend.handletextpad": 0.6,
        "legend.labelspacing": 0.35,
        "legend.borderaxespad": 0.0,

        "lines.linewidth": 1.4,
        "lines.markersize": 3.0,

        "figure.dpi": 200,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def square(ax):
    """Force a truly square plot box regardless of tick-label widths."""
    ax.set_box_aspect(1.0)


def chance_line(ax):
    """The no-discrimination diagonal, recessive."""
    ax.plot([0, 1], [0, 1], linestyle=(0, (2.5, 2.5)),
            color=MUTED, linewidth=0.6, zorder=0)


def energy_label(sqrts):
    return rf"$\sqrt{{s}} = {sqrts}$ GeV"


def auc_legend_label(obs, auc, width="1.6cm"):
    """Legend entry with the AUC value in an aligned second column.

    The value doubles as the visible direct label that the aqua slot needs
    as contrast relief, so it is not decoration -- do not drop it.
    """
    return rf"\makebox[{width}][l]{{{OBS_LABEL[obs]}}}{auc:.3f}"
