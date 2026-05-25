#!/usr/bin/env python3
"""Thick-jet subprocess composition vs eta — live-data version (HERA 300).

Thick jets are defined as Psi(r=0.3) < 0.6 (gluon-like, broad energy
profile), matching src/thinthick/thick_thin_analysis.cc. For each eta bin,
this script counts the fraction of thick jets that come from QQ_Events,
GG_Events, and GQ_Events trees of the regenerated dijet ROOT.

Inputs come from <repo>/data-jets/hera300_kt_dijets/dijets_*.root.
Output PDF goes to plots/efficiency_purity/output/.
"""

import glob
from pathlib import Path

import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.interpolate import make_interp_spline


REPO_ROOT = Path(__file__).resolve().parents[2]
PSI_THICK_CUT = 0.6       # Psi(r=0.3) < PSI_THICK_CUT  -> thick (gluon-like)
SAMPLE_DIR = "hera300_antikt_dijets"


def find_dijet_root(sample_dir):
    base = REPO_ROOT / "data-jets" / sample_dir
    hits = sorted(glob.glob(str(base / "dijets_*.root")))
    return hits[0] if hits else None


def thick_composition(root_path, eta_edges, psi_cut=PSI_THICK_CUT):
    """Return (qq_pct, gg_pct, gq_pct) arrays of length N for thick jets."""
    counts = {c: np.zeros(len(eta_edges) - 1, dtype=int)
              for c in ("QQ_Events", "GG_Events", "GQ_Events")}
    with uproot.open(root_path) as f:
        for cat in counts:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            eta = ak.to_numpy(ak.flatten(t["jet_eta"].array(), axis=1))
            psi = ak.to_numpy(ak.flatten(t["jet_psi03"].array(), axis=1))
            thick = psi < psi_cut
            for i in range(len(eta_edges) - 1):
                m = thick & (eta >= eta_edges[i]) & (eta < eta_edges[i + 1])
                counts[cat][i] = m.sum()
    total = counts["QQ_Events"] + counts["GG_Events"] + counts["GQ_Events"]
    safe = np.where(total > 0, total, 1)
    return (
        100.0 * counts["QQ_Events"] / safe,
        100.0 * counts["GG_Events"] / safe,
        100.0 * counts["GQ_Events"] / safe,
    )

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'
rcParams['font.size'] = 18
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 22
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
# rcParams['legend.fontsize'] = 18
# rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 300
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.width'] = 2.0
rcParams['ytick.major.width'] = 2.0
rcParams['xtick.minor.width'] = 1.0
rcParams['ytick.minor.width'] = 1.0

def create_purity_plot():

    eta_edges   = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])
    eta_centers = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    # ----- legacy paper values (HERA 300, kept for reference) --------------
    # thick_QQ = np.array([42.5, 28.4, 19.0, 14.5])
    # thick_GG = np.array([ 8.2, 20.5, 33.2, 38.6])
    # thick_GQ = np.array([49.4, 51.2, 47.8, 46.9])
    # -----------------------------------------------------------------------

    # Live-data computation from the regenerated HERA 300 dijet sample.
    rf = find_dijet_root(SAMPLE_DIR)
    if rf is None:
        raise SystemExit(f"missing dijet ROOT under data-jets/{SAMPLE_DIR}/")
    thick_QQ, thick_GG, thick_GQ = thick_composition(rf, eta_edges)
    print(f"thick (Psi(0.3) < {PSI_THICK_CUT}) composition per eta bin:")
    print(f"  QQ : {np.round(thick_QQ, 2)}")
    print(f"  GG : {np.round(thick_GG, 2)}")
    print(f"  GQ : {np.round(thick_GQ, 2)}")


    # ==========================================================================
    # CREATE FIGURE AND AXES
    # ==========================================================================
    
    fig, ax = plt.subplots(figsize=(9, 9))
    
    # Set minimalistic style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # ==========================================================================
    # CREATE SMOOTH INTERPOLATION
    # ==========================================================================

    eta_smooth = np.linspace(eta_centers.min(), eta_centers.max(), 300)  # Stay within data bounds

    # Create spline interpolations (k=3 for cubic splines)
    spline_QQ = make_interp_spline(eta_centers, thick_QQ, k=3)
    spline_GG = make_interp_spline(eta_centers, thick_GG, k=3)
    spline_GQ = make_interp_spline(eta_centers, thick_GQ, k=3)

    # Generate smooth curves
    thick_QQ_smooth = spline_QQ(eta_smooth)
    thick_GG_smooth = spline_GG(eta_smooth)
    thick_GQ_smooth = spline_GQ(eta_smooth)
    
    # ==========================================================================
    # PLOT DATA
    # ==========================================================================
    
    # Define colors matching ROOT style
    color_QQ = '#CC0000'  # Red
    color_GG = '#00AA00'  # Green  
    color_GQ = '#0066CC'  # Blue
    
    # Plot smooth curves
    ax.plot(eta_smooth, thick_QQ_smooth, color=color_QQ, linewidth=4, label=r'$\mathrm{QQ}$',zorder=2)
    ax.plot(eta_smooth, thick_GG_smooth, color=color_GG, linewidth=4, label=r'$\mathrm{GG}$',zorder=2)
    ax.plot(eta_smooth, thick_GQ_smooth, color=color_GQ, linewidth=4, label=r'$\mathrm{GQ}$',zorder=2)
    
    # Add the original data points as markers
    ax.scatter(eta_centers, thick_QQ, 
              color=color_QQ, 
              s=150, 
              marker='o',
              edgecolors='black',
              linewidth=1,
              zorder=3)
    
    ax.scatter(eta_centers, thick_GG, 
              color=color_GG, 
              s=150, 
              marker='s',
              edgecolors='black',
              linewidth=1,
              zorder=3)
    
    ax.scatter(eta_centers, thick_GQ, 
              color=color_GQ, 
              s=150, 
              marker='^',
              edgecolors='black',
              linewidth=1,
              zorder=3)
    
    # ==========================================================================
    # AXIS CONFIGURATION
    # ==========================================================================
    
    # Set axis labels with LaTeX
    ax.set_xlabel(r'\textbf{Pseudorapidity} $\boldsymbol{\eta}$', fontsize=26)
    ax.set_ylabel(r'\textbf{Percentage (\%)}', fontsize=26)
    
    # Set axis ranges
    # ax.set_xlim(-1.2, 3.2)
    ax.set_xlim(-1.0, 3.0)
    ax.set_ylim(0, 100)
    
    # Set custom x-ticks to match eta bins
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.set_xticklabels([r'$-1$', r'$0$', r'$1$', r'$2$', r'$3$'])
    
    # Set y-ticks
    ax.set_yticks(np.arange(0, 101, 10))
    
    # Add minor ticks for better readability
    ax.minorticks_on()
    ax.tick_params(which='minor', length=5, width=1.0)
    ax.tick_params(which='major', length=8, width=1.5)
    
    # Add grid for better readability (optional - comment out for more minimal look)
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5, which='major')
    
    # ==========================================================================
    # ADD LEGEND (MINIMALISTIC)
    # ==========================================================================
    
    legend = ax.legend(loc='upper right', 
                      title=r'\textbf{Thick Jets' + '\n' + r'\textbf{ep 300 GeV}',
                      frameon=False,
                      numpoints=1,
                      handlelength=2.5,
                      columnspacing=1,
                      fontsize=26,
                      title_fontsize=26)
    
    # ==========================================================================
    # FINAL ADJUSTMENTS
    # ==========================================================================
    
    plt.tight_layout(pad=1.5)
    
    # Save figures
    out_dir = Path(__file__).resolve().parent / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / "purity_thick_hera.pdf", bbox_inches="tight")
    print(f"Plot saved: {out_dir / 'purity_thick_hera.pdf'}")
    
    return fig, ax

def main():
    fig, ax = create_purity_plot()
    plt.show()

if __name__ == "__main__":
    main()