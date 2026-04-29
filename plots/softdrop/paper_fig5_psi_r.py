#!/usr/bin/env python3
"""Recreate the paper's Figure 5: integrated jet shape Psi(r) vs r
for QQ (red) and GG (blue) at four ep energies, in four eta bins.

Layout (2x2):
    top-left:  -1 < eta < 0
    top-right:  0 < eta < 1
    bot-left:   1 < eta < 1.5
    bot-right:  1.5 < eta < 2

Line styles by sqrt(s):
    64 GeV  -> dotted
    105 GeV -> dashed
    141 GeV -> dash-dot
    300 GeV -> solid

Y-axis is log from 0.05 to 1.0. Computer Modern serif via LaTeX.

Usage:
    python plots/softdrop/paper_fig5_psi_r.py
    python plots/softdrop/paper_fig5_psi_r.py --no-tex   (matplotlib mathtext fallback)
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib
import matplotlib.pyplot as plt


SAMPLES = [
    ("eic64_antikt_dijets",  64),
    ("eic105_antikt_dijets", 105),
    ("eic141_antikt_dijets", 141),
    ("hera300_kt_dijets",    300),
]

LINESTYLE = {64: ":", 105: "--", 141: "-.", 300: "-"}
COLOR_QQ  = "#c1273e"   # paper's red
COLOR_GG  = "#2c5fa8"   # paper's blue

ETA_BINS = [(-1.0, 0.0), (0.0, 1.0), (1.0, 1.5), (1.5, 2.0)]
ETA_TITLES = [r"$-1 < \eta < 0$", r"$0 < \eta < 1$",
              r"$1 < \eta < 1.5$", r"$1.5 < \eta < 2$"]

# Same r grid baked into jetreco_softdrop.cc
R_GRID = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])


def find_root(sample_dir):
    hits = sorted(glob.glob(f"{sample_dir}/dijets_*.root"))
    if not hits:
        hits = sorted(glob.glob(f"{sample_dir}/alljets_*.root"))
    return hits[0] if hits else None


def read_psi_eta(root_path, category):
    """Return (psi_per_jet, eta_per_jet) as numpy arrays of shape
    (N_jets, 10) and (N_jets,). The flat ROOT branch packs 10 floats
    per saved jet in r=0.1..1.0 order; we reshape to (-1, 10)."""
    import awkward as ak
    with uproot.open(root_path) as f:
        key = f"{category}/jets_{category}"
        if key not in f:
            return None, None
        t = f[key]
        if "jet_psi_curve_flat" not in t.keys():
            return None, None
        psi_flat_jagged = t["jet_psi_curve_flat"].array()   # one vector<float> per event
        eta_jagged      = t["jet_eta"].array()              # one vector<float> per event

        # Flatten the event axis to one big 1D awkward, then to numpy.
        psi_all = ak.to_numpy(ak.flatten(psi_flat_jagged, axis=1))
        eta_all = ak.to_numpy(ak.flatten(eta_jagged,      axis=1))

        # Sanity: psi_all has 10 entries per saved jet.
        if psi_all.size != 10 * eta_all.size:
            raise RuntimeError(
                f"psi_flat size {psi_all.size} != 10 * eta size {eta_all.size}")
        psi_per_jet = psi_all.reshape(-1, 10)
    return psi_per_jet, eta_all


def mean_psi_per_bin(psi, eta, lo, hi, min_jets=20):
    """Returns the per-r mean Psi(r) for jets with eta in [lo, hi),
    or None if too few jets."""
    mask = (eta >= lo) & (eta < hi)
    if mask.sum() < min_jets:
        return None
    return psi[mask].mean(axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--out", default="data-jets/cross_energy_paperconfig/paper_fig5_psi_r.pdf")
    ap.add_argument("--no-tex", action="store_true",
                    help="Disable LaTeX rendering (use matplotlib mathtext)")
    args = ap.parse_args()

    if not args.no_tex:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif":  ["Computer Modern Roman"],
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 12,
        })
    else:
        plt.rcParams.update({"font.size": 14})

    # --- load ----------------------------------------------------------
    samples_data = {}
    for sname, sqrts in SAMPLES:
        sdir = Path(args.data_jets) / sname
        rf = find_root(str(sdir))
        if rf is None:
            print(f"  [skip] {sname}: no ROOT"); continue
        psi_qq, eta_qq = read_psi_eta(rf, "QQ_Events")
        psi_gg, eta_gg = read_psi_eta(rf, "GG_Events")
        if psi_qq is None or psi_gg is None:
            print(f"  [skip] {sname}: missing jet_psi_curve branch")
            continue
        samples_data[sqrts] = {
            "QQ": (psi_qq, eta_qq),
            "GG": (psi_gg, eta_gg),
        }
        print(f"  [load] {sname} (sqrt(s)={sqrts}): "
              f"N_QQ={psi_qq.shape[0]:,}, N_GG={psi_gg.shape[0]:,}")

    if not samples_data:
        print("nothing to plot"); sys.exit(1)

    # --- plot ----------------------------------------------------------
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, (lo, hi), title in zip(axes.flat, ETA_BINS, ETA_TITLES):
        for sqrts, dat in samples_data.items():
            ls = LINESTYLE[sqrts]
            for cat, color in [("QQ", COLOR_QQ), ("GG", COLOR_GG)]:
                psi, eta = dat[cat]
                m = mean_psi_per_bin(psi, eta, lo, hi)
                if m is None:
                    continue
                ax.plot(R_GRID, m, color=color, linestyle=ls, linewidth=1.6)

        ax.set_xscale("linear")
        ax.set_yscale("log")
        ax.set_ylim(0.05, 1.10)
        ax.set_xlim(0.05, 1.05)
        ax.set_xlabel(r"$r$")
        ax.set_ylabel(r"$\Psi(r)$")
        ax.grid(False)
        # eta label in bottom-right corner of each panel (paper convention)
        ax.text(0.96, 0.04, title, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=14)

        # Single combined legend per panel — paper-style:
        # red "QQ" / blue "GG" header, then 4 energies (black, linestyle by sqrts)
        header_handles = [
            Line2D([0], [0], color=COLOR_QQ, lw=1.8),
            Line2D([0], [0], color=COLOR_GG, lw=1.8),
        ]
        # An invisible spacer between header and energy block
        spacer = Line2D([0], [0], color="none", lw=0)
        energy_handles = [
            Line2D([0], [0], color="black", linestyle=LINESTYLE[s], lw=1.6)
            for s in (64, 105, 141, 300)
        ]
        all_handles = header_handles + [spacer] + energy_handles
        all_labels  = ["QQ", "GG", "", "64 GeV", "105 GeV", "141 GeV", "300 GeV"]
        leg = ax.legend(all_handles, all_labels, loc="lower right",
                        bbox_to_anchor=(0.99, 0.18), frameon=False,
                        handlelength=2.0, fontsize=11,
                        title="QQ / GG", title_fontsize=12)
        # Color the header title to match the paper's red/blue split
        leg.get_title().set_color("black")
        # Manually colour the QQ / GG entries' text
        for txt, col in zip(leg.get_texts()[:2], [COLOR_QQ, COLOR_GG]):
            txt.set_color(col)

    fig.tight_layout()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    print(f"\nWrote {out_path.resolve()}")


if __name__ == "__main__":
    main()
