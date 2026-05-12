#!/usr/bin/env python3
"""Diagnostic: rho(r) at 1.5 < eta < 2 for HERA-300 reco'd two ways
(kT EtMin=17 — the canonical paper-config, vs anti-kT EtMin=10 — matched
to the EIC samples) plus the three EIC samples at anti-kT EtMin=10.

Tests the hypothesis that the offset between HERA and EIC in the
Direct vs Resolved rho(r) figure at forward eta is driven by the
algorithm + EtMin choice rather than sqrt(s). If HERA-antikt-EtMin10
overlays the EIC curves, hypothesis confirmed.

Usage:
    python plots/jet_shapes/compare_hera_algo_rho.py
"""

import glob
from pathlib import Path

import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline


R_GRID  = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
DELTA_R = 0.1
RHO_R   = R_GRID - 0.5 * DELTA_R           # 10 centres: 0.05, 0.15, ..., 0.95

ETA_LO, ETA_HI = 1.5, 2.0      # forward eta bin where the HERA-EIC gap was largest

REPO = Path(__file__).resolve().parents[2]

SAMPLES = [
    ("HERA 300, kT, EtMin=17 (canonical)",     "hera300_kt_dijets",      "#000000", "-"),
    ("HERA 300, anti-kT, EtMin=10 (new)",      "hera300_antikt_dijets",  "#7f7f7f", "--"),
    ("EIC 141, anti-kT, EtMin=10",             "eic141_antikt_dijets",   "#1f77b4", "-"),
    ("EIC 105, anti-kT, EtMin=10",             "eic105_antikt_dijets",   "#2ca02c", "-"),
    ("EIC 64,  anti-kT, EtMin=10",             "eic64_antikt_dijets",    "#d62728", "-"),
]


def find_dijet_root(sample_dir):
    base = REPO / "data-jets" / sample_dir
    hits = sorted(glob.glob(str(base / "dijets_*.root")))
    return hits[0] if hits else None


def mean_rho_in_eta_bin(root_path, lo, hi):
    """Average Psi(r) over jets in eta bin (across all 3 categories),
    then take finite-difference rho(r)."""
    psi_all, eta_all = [], []
    with uproot.open(root_path) as f:
        for cat in ("QQ_Events", "GG_Events", "GQ_Events"):
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            psi_flat = ak.to_numpy(ak.flatten(
                t["jet_psi_curve_flat"].array(), axis=1))
            eta = ak.to_numpy(ak.flatten(t["jet_eta"].array(), axis=1))
            if psi_flat.size == 0:
                continue
            psi_all.append(psi_flat.reshape(-1, 10))
            eta_all.append(eta)
    psi = np.vstack(psi_all)
    eta = np.concatenate(eta_all)
    m = (eta >= lo) & (eta < hi)
    if m.sum() < 50:
        return None, m.sum()
    psi_mean = psi[m].mean(axis=0)
    rho = np.diff(np.concatenate(([0.0], psi_mean))) / DELTA_R
    return rho, m.sum()


def main():
    plt.rcParams.update({
        "font.size": 16, "font.family": "serif",
        "font.serif": ["Computer Modern"], "text.usetex": True,
        "axes.linewidth": 1.2, "axes.spines.top": False,
        "axes.spines.right": False, "lines.linewidth": 2,
        "legend.frameon": False, "figure.facecolor": "white",
    })

    fig, ax = plt.subplots(figsize=(8, 8))
    r_smooth = np.linspace(RHO_R.min(), RHO_R.max(), 300)

    for label, sample_dir, color, ls in SAMPLES:
        rf = find_dijet_root(sample_dir)
        if rf is None:
            print(f"  [skip] {sample_dir}: no ROOT")
            continue
        rho, n = mean_rho_in_eta_bin(rf, ETA_LO, ETA_HI)
        if rho is None:
            print(f"  [skip] {sample_dir}: too few jets ({n})")
            continue
        sp = make_interp_spline(RHO_R, rho, k=3)(r_smooth)
        ax.plot(r_smooth, sp, color=color, linestyle=ls, linewidth=2.5,
                label=f"{label}  (N={n:,})")
        print(f"  [load] {sample_dir}: N={n:,}, rho peak={rho.max():.3f}")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$r$", fontsize=22)
    ax.set_ylabel(r"$\rho(r)$", fontsize=22, labelpad=10)
    ax.set_title(r"$\rho(r)$ for $1.5 < \eta < 2$", fontsize=18, pad=10)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.legend(loc="upper right", fontsize=12, bbox_to_anchor=(0.98, 0.95))
    plt.tight_layout()

    out = Path(__file__).resolve().parent / "output" / "compare_hera_algo_rho_forward.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nPlot saved: {out}")


if __name__ == "__main__":
    main()
