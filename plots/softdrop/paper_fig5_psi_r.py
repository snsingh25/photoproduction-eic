#!/usr/bin/env python3
"""Integrated jet shape Psi(r) — paper Figure 5 reproduction.

Uses the styling and structure of the legacy
ph-new/subprocjets/jets-int/plot_int/trypython/plot_intshapes.py
(Computer Modern via LaTeX, Tufte-style spines off, log-y, smoothed
curves) but reads the underlying numbers from the regenerated dijet
ROOT files in data-jets/, so the values stay in sync with the live
analysis pipeline.

Per-bin output mirrors the original — four standalone PDFs, one per
eta bin — so they can be dropped straight into the paper.

Usage:
    python plots/softdrop/paper_fig5_psi_r.py
    python plots/softdrop/paper_fig5_psi_r.py --no-tex   (mathtext fallback)
    python plots/softdrop/paper_fig5_psi_r.py --out-dir <path>
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import awkward as ak
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import make_interp_spline


# Same r grid baked into jetreco_softdrop.cc.
R_GRID = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

SAMPLES = [
    ("eic64_antikt_dijets",  64),
    ("eic105_antikt_dijets", 105),
    ("eic141_antikt_dijets", 141),
    ("hera300_kt_dijets",    300),
]

ETA_BINS = [
    (-1.0, 0.0, r"$-1 < \eta < 0$",    "injetshapesminus1to0.pdf"),
    ( 0.0, 1.0, r"$0 < \eta < 1$",     "injetshapes0to1.pdf"),
    ( 1.0, 1.5, r"$1 < \eta < 1.5$",   "injetshapes1to1p5.pdf"),
    ( 1.5, 2.0, r"$1.5 < \eta < 2$",   "injetshapes1p5to2.pdf"),
]

QQ_COLOR = "#d62728"   # standard MPL red
GG_COLOR = "#1f77b4"   # standard MPL blue


def find_root(sample_dir):
    for pat in ("dijets_*.root", "alljets_*.root"):
        hits = sorted(glob.glob(f"{sample_dir}/{pat}"))
        if hits:
            return hits[0]
    return None


def read_psi_eta(root_path, category):
    """Return (psi_per_jet, eta_per_jet): shapes (N_jets, 10), (N_jets,)."""
    with uproot.open(root_path) as f:
        key = f"{category}/jets_{category}"
        if key not in f or "jet_psi_curve_flat" not in f[key].keys():
            return None, None
        t = f[key]
        psi_flat = ak.to_numpy(ak.flatten(
            t["jet_psi_curve_flat"].array(), axis=1))
        eta = ak.to_numpy(ak.flatten(
            t["jet_eta"].array(), axis=1))
    if psi_flat.size != 10 * eta.size:
        raise RuntimeError(
            f"psi_flat size {psi_flat.size} != 10 * eta size {eta.size}")
    return psi_flat.reshape(-1, 10), eta


def mean_psi(psi, eta, lo, hi, min_jets=20):
    m = (eta >= lo) & (eta < hi)
    if m.sum() < min_jets:
        return None
    return psi[m].mean(axis=0)


def setup_style(use_tex=True):
    plt.style.use("default")
    plt.rcParams.update({
        "font.size": 24,
        "font.family": "serif",
        "font.serif": ["Computer Modern"] if use_tex else ["DejaVu Serif"],
        "text.usetex": use_tex,
        "axes.linewidth": 1.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.major.size": 6,
        "ytick.minor.size": 3,
        "lines.linewidth": 2,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def create_plot(qq_curves, gg_curves, eta_label, file_path):
    """Plot one eta bin: red QQ + blue GG, one curve per energy.

    qq_curves / gg_curves are dicts {sqrts: numpy(10,)} or with None
    entries for samples without enough statistics in this bin.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    # Smooth interpolation matching the legacy plot.
    r_smooth = np.linspace(R_GRID.min(), R_GRID.max(), 300)
    def smooth(y):
        return make_interp_spline(R_GRID, y, k=3)(r_smooth)

    LINESTYLE = {64: ":", 105: "--", 141: "-.", 300: "-"}

    # QQ family (red) — only the QQ curves carry energy labels for the legend.
    for sqrts in (64, 105, 141, 300):
        y = qq_curves.get(sqrts)
        if y is None:
            continue
        ax.plot(r_smooth, smooth(y), color=QQ_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2,
                label=f"ep {sqrts} GeV")

    # GG family (blue) — same line styles, no separate labels.
    for sqrts in (64, 105, 141, 300):
        y = gg_curves.get(sqrts)
        if y is None:
            continue
        ax.plot(r_smooth, smooth(y), color=GG_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)

    # Axis: linear x, log y, paper-style ticks.
    ax.set_xlim(0.05, 1.0)
    ax.set_yscale("log")
    ax.set_yticks([0.05, 0.1, 1.0])
    ax.set_yticklabels(["0.05", "0.1", "1.0"])
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)

    ax.set_xlabel(r"$r$", fontsize=28)
    ax.set_ylabel(r"$\Psi(r)$", fontsize=28, labelpad=10)
    ax.xaxis.set_label_coords(0.5, -0.08)
    ax.yaxis.set_label_coords(-0.08, 0.5)

    ax.tick_params(axis="both", which="major", labelsize=20)
    ax.tick_params(axis="both", which="minor", labelsize=10)
    ax.minorticks_on()

    # Legend: 4 line-style swatches in black (energy key).
    legend_elements = [
        plt.Line2D([0], [0], color="black", linestyle=":",  linewidth=2, label="64 GeV"),
        plt.Line2D([0], [0], color="black", linestyle="--", linewidth=2, label="105 GeV"),
        plt.Line2D([0], [0], color="black", linestyle="-.", linewidth=2, label="141 GeV"),
        plt.Line2D([0], [0], color="black", linestyle="-",  linewidth=2, label="300 GeV"),
    ]
    ax.legend(handles=legend_elements, loc="lower right",
              bbox_to_anchor=(0.95, 0.15), fontsize=22)

    # QQ / GG colour key, positioned mid-panel as in the paper.
    ax.text(0.60, 0.50, "QQ", transform=ax.transAxes,
            fontsize=24, color=QQ_COLOR, weight="bold")
    ax.text(0.70, 0.50, " / ", transform=ax.transAxes,
            fontsize=24, color="black", weight="bold")
    ax.text(0.75, 0.50, "GG", transform=ax.transAxes,
            fontsize=24, color=GG_COLOR, weight="bold")

    # Eta label (bottom-right corner like the paper).
    ax.text(0.65, 0.10, eta_label, transform=ax.transAxes, fontsize=24)

    plt.tight_layout()
    plt.savefig(file_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  -> {file_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--out-dir",
                    default="data-jets/cross_energy_paperconfig/fig5_panels")
    ap.add_argument("--no-tex", action="store_true",
                    help="Use mathtext instead of LaTeX (for systems without "
                         "Computer Modern installed).")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_style(use_tex=not args.no_tex)

    # --- load ----------------------------------------------------------
    print("Loading dijet samples:")
    sample_data = {}
    for sname, sqrts in SAMPLES:
        rf = find_root(str(Path(args.data_jets) / sname))
        if rf is None:
            print(f"  [skip] {sname}"); continue
        psi_qq, eta_qq = read_psi_eta(rf, "QQ_Events")
        psi_gg, eta_gg = read_psi_eta(rf, "GG_Events")
        if psi_qq is None or psi_gg is None:
            print(f"  [skip] {sname}: missing jet_psi_curve_flat"); continue
        sample_data[sqrts] = ((psi_qq, eta_qq), (psi_gg, eta_gg))
        print(f"  [load] {sname} sqrt(s)={sqrts}: "
              f"N_QQ={psi_qq.shape[0]:,} N_GG={psi_gg.shape[0]:,}")

    if not sample_data:
        print("nothing to plot"); sys.exit(1)

    # --- plot per eta bin ----------------------------------------------
    print(f"\nWriting per-eta-bin plots to {out_dir.resolve()}/")
    for lo, hi, label, fname in ETA_BINS:
        qq_curves = {}
        gg_curves = {}
        for sqrts, ((psi_qq, eta_qq), (psi_gg, eta_gg)) in sample_data.items():
            qq_curves[sqrts] = mean_psi(psi_qq, eta_qq, lo, hi)
            gg_curves[sqrts] = mean_psi(psi_gg, eta_gg, lo, hi)
        create_plot(qq_curves, gg_curves, label, str(out_dir / fname))

    print("\nAll panels written.")


if __name__ == "__main__":
    main()
