#!/usr/bin/env python3
"""Side-by-side comparison of two rho(r) discretizations.

Upper-edge convention (what the active scripts currently use):
    rho_i = (Psi(r_i) - Psi(r_{i-1})) / Delta r ,   r_i = 0.1, 0.2, ..., 1.0
    -- 10 points; first bin = solid disk (0, 0.1] using Psi(0)=0.
    -- This is the paper's reported convention.

Centered-bin convention:
    rho_k = (Psi(r_k + 0.05) - Psi(r_k - 0.05)) / Delta r ,
             r_k = 0.1, 0.2, ..., 0.9     (9 points)
    -- Each rho value sits at the centre of a Delta r-wide annulus.
    -- Requires Psi at half-integer points, obtained by linear
       interpolation between our stored Psi values at 0.1, 0.2, ..., 1.0
       and Psi(0)=0; equivalent to the standard centered finite
       difference (Psi(r_k + 0.1) - Psi(r_k - 0.1)) / (2 Delta r).

For each (sqrt(s), eta-bin) the script overlays QQ and GG in both
conventions; the change of binning slightly shifts and smooths the
curve but the QQ-vs-GG separation is unchanged.

Output: plots/jet_shapes/output/compare_rho_conventions/
        compare_<sample>_eta_<bin>.pdf
"""

import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline


R_GRID  = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
DELTA_R = 0.1

UPPER_EDGE_R = R_GRID                              # 10 points at upper edges
CENTERED_R   = R_GRID[:-1]                         # 9 points at 0.1..0.9

ETA_BINS = [
    (-1.0, 0.0, r"$-1 < \eta < 0$",  "minus1to0"),
    ( 0.0, 1.0, r"$0 < \eta < 1$",   "0to1"),
    ( 1.0, 1.5, r"$1 < \eta < 1.5$", "1to1p5"),
    ( 1.5, 2.0, r"$1.5 < \eta < 2$", "1p5to2"),
]

QQ_COLOR = "#c0392b"
GG_COLOR = "#2980b9"

REPO = Path(__file__).resolve().parents[2]


def find_dijet_root(sample_dir):
    base = REPO / "data-jets" / sample_dir
    hits = sorted(glob.glob(str(base / "dijets_*.root")))
    return hits[0] if hits else None


def mean_psi(root_path, category, lo, hi, min_jets=20):
    with uproot.open(root_path) as f:
        t = f[f"{category}/jets_{category}"]
        eta = ak.to_numpy(ak.flatten(t["jet_eta"].array(), axis=1))
        psi_flat = ak.to_numpy(ak.flatten(
            t["jet_psi_curve_flat"].array(), axis=1))
    psi = psi_flat.reshape(-1, 10)
    m = (eta >= lo) & (eta < hi)
    if m.sum() < min_jets:
        return None
    return psi[m].mean(axis=0)


def rho_upper_edge(psi_mean):
    """Upper-edge convention: rho at r = 0.1, 0.2, ..., 1.0."""
    return np.diff(np.concatenate(([0.0], psi_mean))) / DELTA_R   # (10,)


def rho_centered(psi_mean):
    """Centered-bin convention: rho at r = 0.1, 0.2, ..., 0.9.
    Uses centered finite difference, equivalent to linear-interpolating
    Psi to half-integer points."""
    psi_full = np.concatenate(([0.0], psi_mean))  # (11,) at r=0, 0.1, ..., 1.0
    # centered diff at r=r_k uses psi at r_k+/-0.1, divided by 2*Delta r
    return (psi_full[2:11] - psi_full[0:9]) / (2 * DELTA_R)        # (9,)


def main():
    plt.rcParams.update({
        "font.size": 16, "font.family": "serif",
        "font.serif": ["Computer Modern"], "text.usetex": True,
        "axes.linewidth": 1.2, "axes.spines.top": False,
        "axes.spines.right": False, "lines.linewidth": 1.8,
        "legend.frameon": False, "figure.facecolor": "white",
    })

    out_dir = REPO / "plots" / "jet_shapes" / "output" / "compare_rho_conventions"
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_dir = "hera300_antikt_dijets"
    jet_root = find_dijet_root(sample_dir)
    if jet_root is None:
        print(f"missing ROOT for {sample_dir}"); sys.exit(1)
    print(f"Loading {sample_dir}")

    for lo, hi, eta_lbl, slug in ETA_BINS:
        psi_qq = mean_psi(jet_root, "QQ_Events", lo, hi)
        psi_gg = mean_psi(jet_root, "GG_Events", lo, hi)
        if psi_qq is None or psi_gg is None:
            print(f"  [skip] eta={slug}: too few jets")
            continue

        fig, ax = plt.subplots(figsize=(9, 7))

        # Smooth splines for visual comparison
        r_smooth_ue = np.linspace(UPPER_EDGE_R.min(), UPPER_EDGE_R.max(), 300)
        r_smooth_ct = np.linspace(CENTERED_R.min(),  CENTERED_R.max(),  300)

        for psi, color, label in [(psi_qq, QQ_COLOR, "QQ"),
                                   (psi_gg, GG_COLOR, "GG")]:
            ue = rho_upper_edge(psi)
            ct = rho_centered(psi)
            ax.plot(r_smooth_ue, make_interp_spline(UPPER_EDGE_R, ue, k=3)(r_smooth_ue),
                    color=color, linestyle="-", linewidth=2,
                    label=f"{label}: upper-edge (current)")
            ax.plot(r_smooth_ct, make_interp_spline(CENTERED_R, ct, k=3)(r_smooth_ct),
                    color=color, linestyle="--", linewidth=2,
                    label=f"{label}: centered (proposed)")
            # also drop the actual data points to make the discretization explicit
            ax.plot(UPPER_EDGE_R, ue, "o", color=color, markersize=6,
                    markeredgecolor="black", markeredgewidth=0.5, zorder=5)
            ax.plot(CENTERED_R, ct, "s", color=color, markersize=6,
                    markeredgecolor="black", markeredgewidth=0.5, zorder=5,
                    markerfacecolor="none")

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel(r"$r$", fontsize=20)
        ax.set_ylabel(r"$\rho(r)$", fontsize=20)
        ax.set_title(rf"HERA 300, anti-$k_T$, dijets ; {eta_lbl}", fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=14)
        ax.minorticks_on()
        ax.legend(loc="upper right", fontsize=12)

        out = out_dir / f"compare_hera300_eta_{slug}.pdf"
        plt.tight_layout()
        plt.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"  -> {out}")


if __name__ == "__main__":
    main()
