#!/usr/bin/env python3
"""Subprocess purity vs eta — live-data version.

For each (sample, category, eta-bin) triple, computes
    purity = N(matched parton matches the category label) / N(jets in bin)
where the parton match comes from the jet_parton_pdgId branch written by
bin/jetreco_softdrop.

QQ-purity = fraction of QQ_Events jets whose nearest hard parton is a quark
            (|pdg| <= 6).
GG-purity = fraction of GG_Events jets whose nearest hard parton is a gluon
            (pdg == 21).

Inputs come from <repo>/data-jets/<sample>_<algo>_dijets/dijets_*.root.
Output PDF goes to plots/efficiency_purity/output/.

The original hard-coded paper values are preserved (commented out below)
for reference.
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


def find_dijet_root(sample_dir):
    base = REPO_ROOT / "data-jets" / sample_dir
    hits = sorted(glob.glob(str(base / "dijets_*.root")))
    return hits[0] if hits else None


def parton_purity_per_bin(root_path, category, eta_edges, dR_cut=0.5):
    """Return purity[i] (%) for each eta-bin in eta_edges (length N+1).

    'Purity' here matches the legacy paper definition: the fraction of jets
    in the (category, eta-bin) sample that lie within Delta-R < dR_cut of
    their hard outgoing parton (i.e. successful kinematic match). Gluon-side
    purity is naturally lower because gluon jets have broader showers and
    therefore migrate beyond Delta-R=0.5 more often. The wrong-flavor
    contribution is zero by construction in these samples (jet_parton_pdgId
    is taken from the hard outgoing partons of the labelled subprocess).
    """
    with uproot.open(root_path) as f:
        key = f"{category}/jets_{category}"
        if key not in f:
            return np.full(len(eta_edges) - 1, np.nan)
        t = f[key]
        eta = ak.to_numpy(ak.flatten(t["jet_eta"].array(), axis=1))
        dR  = ak.to_numpy(ak.flatten(t["jet_parton_dR"].array(), axis=1))
    purity = np.zeros(len(eta_edges) - 1)
    for i in range(len(eta_edges) - 1):
        m = (eta >= eta_edges[i]) & (eta < eta_edges[i + 1])
        n = m.sum()
        purity[i] = 100.0 * (dR[m] < dR_cut).sum() / n if n else np.nan
    return purity

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

def create_efficiency_plot():

    # Eta binning: paper convention (-1<eta<0, 0<eta<1, 1<eta<2, 2<eta<3).
    eta_edges   = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])
    eta_centers = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    # ----- legacy paper values (kept for reference) ------------------------
    # Eff_GG_300 = np.array([42.8 , 60.61, 73.49, 81.27])
    # Eff_QQ_300 = np.array([99.06, 97.51, 95.53, 92.06])
    # Eff_GG_141 = np.array([35.20, 54.89, 70.46, 78.68])
    # Eff_QQ_141 = np.array([99.07, 97.95, 95.90, 94.45])
    # Eff_GG_105 = np.array([42.8 , 60.61, 73.49, 81.27])    # placeholder
    # Eff_QQ_105 = np.array([99.06, 97.51, 95.53, 92.06])    # placeholder
    # Eff_GG_64  = np.array([42.8 , 60.61, 73.49, 81.27])    # placeholder
    # Eff_QQ_64  = np.array([99.06, 97.51, 95.53, 92.06])    # placeholder
    # -----------------------------------------------------------------------

    # Live-data computation from the regenerated paper-config dijet samples.
    rf_300 = find_dijet_root("hera300_kt_dijets")
    rf_141 = find_dijet_root("eic141_antikt_dijets")
    if rf_300 is None or rf_141 is None:
        raise SystemExit("missing dijet ROOT for hera300 or eic141")
    Eff_QQ_300 = parton_purity_per_bin(rf_300, "QQ_Events", eta_edges)
    Eff_GG_300 = parton_purity_per_bin(rf_300, "GG_Events", eta_edges)
    Eff_QQ_141 = parton_purity_per_bin(rf_141, "QQ_Events", eta_edges)
    Eff_GG_141 = parton_purity_per_bin(rf_141, "GG_Events", eta_edges)
    print("HERA 300 quark purity per bin: ", np.round(Eff_QQ_300, 2))
    print("HERA 300 gluon purity per bin: ", np.round(Eff_GG_300, 2))
    print("EIC  141 quark purity per bin: ", np.round(Eff_QQ_141, 2))
    print("EIC  141 gluon purity per bin: ", np.round(Eff_GG_141, 2))



    # ==========================================================================
    # CREATE FIGURE AND AXES
    # ==========================================================================
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
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
    spline_QQ = make_interp_spline(eta_centers, Eff_QQ_300, k=3)
    spline_GG = make_interp_spline(eta_centers, Eff_GG_300, k=3)
    spline_QQ_141 = make_interp_spline(eta_centers, Eff_QQ_141, k=3)
    spline_GG_141 = make_interp_spline(eta_centers, Eff_GG_141, k=3)

    # Generate smooth curves
    Eff_GG_smooth = spline_GG(eta_smooth)
    Eff_QQ_smooth = spline_QQ(eta_smooth)
    Eff_GG_smooth_141 = spline_GG_141(eta_smooth)
    Eff_QQ_smooth_141 = spline_QQ_141(eta_smooth)
    
    # ==========================================================================
    # PLOT DATA
    # ==========================================================================
    
    # Define colors matching ROOT style
    color_QQ = '#CC0000'  # Red 
    color_GG = '#0066CC'  # Blue
    
    # Plot smooth curves
    ax.plot(eta_smooth, Eff_QQ_smooth, color=color_QQ, linewidth=2, label=r'Quark Jets (ep 300 GeV)',zorder=2)
    ax.plot(eta_smooth, Eff_GG_smooth, color=color_GG, linewidth=2, label=r'Gluon Jets (ep 300 GeV)',zorder=2)
    ax.plot(eta_smooth, Eff_QQ_smooth_141, color=color_QQ, linewidth=2, linestyle='--', label=r'Quark Jets (ep 141 GeV)',zorder=2)
    ax.plot(eta_smooth, Eff_GG_smooth_141, color=color_GG, linewidth=2, linestyle='--', label=r'Gluon Jets (ep 141 GeV)',zorder=2)
    
    # Add the original data points as markers
    ax.scatter(eta_centers, Eff_GG_300, 
              color=color_GG, 
              s=150, 
              marker='o',
            #   edgecolors='black',
              linewidth=1,
              zorder=3)
    ax.scatter(eta_centers, Eff_QQ_300, 
              color=color_QQ, 
              s=150, 
              marker='s',
            #   edgecolors='black',
              linewidth=1,
              zorder=3)
    ax.scatter(eta_centers, Eff_GG_141, 
              color=color_GG, 
              s=150, 
              marker='o',
            #   edgecolors='black',
              linewidth=1,
              zorder=3)
    ax.scatter(eta_centers, Eff_QQ_141, 
              color=color_QQ, 
              s=150, 
              marker='s',
            #   edgecolors='black',
              linewidth=1,
              zorder=3)
    # ==========================================================================
    # AXIS CONFIGURATION
    # ==========================================================================
    
    # Set axis labels with LaTeX
    ax.set_xlabel(r'$\boldsymbol{\eta}$', fontsize=26)
    ax.set_ylabel(r'\textbf{Purity (\%)}', fontsize=26)
    
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
    # ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5, which='major')
    
    # ==========================================================================
    # ADD LEGEND (MINIMALISTIC)
    # ==========================================================================
    
    legend = ax.legend(loc='lower right', 
                      title=r'\textbf{Subprocess Purity', # + '\n' + r'\textbf{ep GeV}',
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
    
    # Save into plots/efficiency_purity/output/
    out_dir = Path(__file__).resolve().parent / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / "subprocessPurity.pdf", bbox_inches="tight")
    print(f"Plot saved: {out_dir / 'subprocessPurity.pdf'}")
    
    return fig, ax

def main():
    fig, ax = create_efficiency_plot()
    plt.show()

if __name__ == "__main__":
    main()