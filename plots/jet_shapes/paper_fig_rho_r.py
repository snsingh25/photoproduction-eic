#!/usr/bin/env python3
"""Differential jet shape rho(r) — live-data version, QQ vs GG split.

Sibling of paper_fig5_psi_r.py. The integrated jet shape Psi(r) is
already stored per-jet as jet_psi_curve_flat (10 r-points,
r = 0.1, 0.2, ..., 1.0). The differential jet shape is the derivative
of Psi(r):

    rho(r_i_center) = ( Psi(r_i) - Psi(r_{i-1}) ) / Delta r ,    Delta r = 0.1

Psi(0) is taken to be 0, so the first rho point corresponds to the
annulus 0 < r < 0.1 and is centred at r = 0.05. The subsequent points
are centred at r = 0.15, 0.25, ..., 0.95.

Styling matches paper_fig5_psi_r.py: Computer Modern via LaTeX,
Tufte-style spines off, four PDFs (one per eta bin), QQ red / GG blue,
energy-coded line styles.

Usage:
    python plots/jet_shapes/paper_fig_rho_r.py
    python plots/jet_shapes/paper_fig_rho_r.py --no-tex     (mathtext fallback)
    python plots/jet_shapes/paper_fig_rho_r.py --out-dir <path>
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
from scipy.interpolate import make_interp_spline


# Same r grid baked into jetreco_softdrop.cc — Psi(r) stored at these r values.
R_GRID    = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
DELTA_R   = 0.1
# Annulus centres: rho(r) is reported at the midpoint of each (r-DR, r] bin,
# so the 10 centres span r = 0.05, 0.15, ..., 0.95. The first point uses
# Psi(0)=0 as the lower-edge boundary — matching the paper convention so
# the GG peak structure near r ~ 0.2 is captured.
RHO_R     = R_GRID                       # paper convention: report rho at the UPPER edge of each annulus (r = 0.1, 0.2, ..., 1.0)

SAMPLES = [
    ("eic64_antikt_dijets",   64),
    ("eic105_antikt_dijets",  105),
    ("eic141_antikt_dijets",  141),
    # HERA uses anti-kT EtMin=10 (matching the EIC samples).
    ("hera300_antikt_dijets", 300),
]

ETA_BINS = [
    (-1.0, 0.0, r"$-1 < \eta < 0$",    "diffjetshapesminus1to0.pdf"),
    ( 0.0, 1.0, r"$0 < \eta < 1$",     "diffjetshapes0to1.pdf"),
    ( 1.0, 1.5, r"$1 < \eta < 1.5$",   "diffjetshapes1to1p5.pdf"),
    ( 1.5, 2.0, r"$1.5 < \eta < 2$",   "diffjetshapes1p5to2.pdf"),
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


def mean_rho(psi, eta, lo, hi, min_jets=20):
    """Per-eta-bin mean rho(r) computed from the Psi(r) curve.

    rho_i = ( Psi(r_i) - Psi(r_{i-1}) ) / Delta r , i = 1..10
    with Psi(r_0) = 0 for the first bin, matching the paper convention.
    """
    m = (eta >= lo) & (eta < hi)
    if m.sum() < min_jets:
        return None
    psi_mean = psi[m].mean(axis=0)                          # shape (10,)
    return np.diff(np.concatenate(([0.0], psi_mean))) / DELTA_R   # shape (10,)


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

    # Smooth interpolation matching paper_fig5_psi_r.py.
    r_smooth = np.linspace(RHO_R.min(), RHO_R.max(), 300)
    def smooth(y):
        return make_interp_spline(RHO_R, y, k=3)(r_smooth)

    LINESTYLE = {64: ":", 105: "--", 141: "-.", 300: "-"}

    # QQ family (red).
    for sqrts in (64, 105, 141, 300):
        y = qq_curves.get(sqrts)
        if y is None:
            continue
        ax.plot(r_smooth, smooth(y), color=QQ_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)

    # GG family (blue).
    for sqrts in (64, 105, 141, 300):
        y = gg_curves.get(sqrts)
        if y is None:
            continue
        ax.plot(r_smooth, smooth(y), color=GG_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)

    # Axis: linear x, linear y for differential shape (paper convention).
    # rho is reported starting at r = 0.15, but the x-axis extends to
    # r = 0 with a visible gap so the missing low-r region is explicit.
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)

    ax.set_xlabel(r"$r$", fontsize=28)
    ax.set_ylabel(r"$\rho(r)$", fontsize=28, labelpad=10)
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
    ax.legend(handles=legend_elements, loc="upper right",
              bbox_to_anchor=(0.95, 0.95), fontsize=22)

    # QQ / GG colour key, positioned mid-panel.
    ax.text(0.60, 0.60, "QQ", transform=ax.transAxes,
            fontsize=24, color=QQ_COLOR, weight="bold")
    ax.text(0.70, 0.60, " / ", transform=ax.transAxes,
            fontsize=24, color="black", weight="bold")
    ax.text(0.75, 0.60, "GG", transform=ax.transAxes,
            fontsize=24, color=GG_COLOR, weight="bold")

    # Eta label — placed in the lower-right empty region (after the
    # peak has decayed) so it doesn't overlap the QQ peak at r ~ 0.1.
    ax.text(0.55, 0.55, eta_label, transform=ax.transAxes, fontsize=24)

    plt.tight_layout()
    plt.savefig(file_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  -> {file_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--out-dir",
                    default=str(Path(__file__).resolve().parent
                                / "output/fig_rho_panels"))
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
            qq_curves[sqrts] = mean_rho(psi_qq, eta_qq, lo, hi)
            gg_curves[sqrts] = mean_rho(psi_gg, eta_gg, lo, hi)
        create_plot(qq_curves, gg_curves, label, str(out_dir / fname))

    print("\nAll panels written.")


if __name__ == "__main__":
    main()
