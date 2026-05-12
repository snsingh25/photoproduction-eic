#!/usr/bin/env python3
"""Integrated jet shape Psi(r) and differential jet shape rho(r),
split by Direct vs Resolved photoproduction events.

Sibling of paper_fig5_psi_r.py (Psi(r) QQ vs GG) and
paper_fig_rho_r.py (rho(r) QQ vs GG). Instead of splitting by
outgoing-parton flavour (QQ / GG / GQ), this script splits jets
by photon-production mode:

    Direct   :  processCode in [271..284]   (gamma-parton vertex)
    Resolved :  processCode in [111..124]   (partonic gamma vertex)

processCode is per-event in the source ROOT (data/<sample>/*.root,
QQ_Events/GG_Events/GQ_Events trees, each row has eventID +
processCode); we build a global eventID -> processCode dict and
look it up per jet via the eventID field on each jet-level row.
No re-reco needed.

Output: 8 PDFs into
    plots/jet_shapes/output/fig_psi_rho_direct_resolved/
    psi_directresolved_<etabin>.pdf   (4 panels)
    rho_directresolved_<etabin>.pdf   (4 panels)

Each plot overlays Direct (red) and Resolved (blue) for all 4 sqrt(s)
values (line style encodes energy: 64 GeV ":", 105 GeV "--",
141 GeV "-.", 300 GeV "-").

Usage:
    python plots/jet_shapes/paper_fig_psi_rho_direct_resolved.py
    python plots/jet_shapes/paper_fig_psi_rho_direct_resolved.py --no-tex
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


R_GRID  = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
DELTA_R = 0.1
# rho(r) at annulus centres. We drop the first bin (which would rely on
# Psi(0)=0 with no measurement of Psi at sub-0.1 radii), so the lowest r
# reported is r = 0.15.
RHO_R   = R_GRID[1:] - 0.5 * DELTA_R       # 9 centres: 0.15, 0.25, ..., 0.95

SAMPLES = [
    ("eic64_antikt_dijets",  "eic64",   64),
    ("eic105_antikt_dijets", "eic105",  105),
    ("eic141_antikt_dijets", "eic141",  141),
    # HERA uses anti-kT EtMin=10 (matching the EIC samples) so the four
    # sqrt(s) values are compared with the same algorithm and per-jet
    # ET threshold. The canonical kT EtMin=17 HERA sample lives in
    # data-jets/hera300_kt_dijets/ and is what Tables I/II/III in the
    # paper are reconstructed from -- it is intentionally NOT used here.
    ("hera300_antikt_dijets", "hera300", 300),
]

ETA_BINS = [
    (-1.0, 0.0, r"$-1 < \eta < 0$",    "minus1to0"),
    ( 0.0, 1.0, r"$0 < \eta < 1$",     "0to1"),
    ( 1.0, 1.5, r"$1 < \eta < 1.5$",   "1to1p5"),
    ( 1.5, 2.0, r"$1.5 < \eta < 2$",   "1p5to2"),
]

DIRECT_COLOR   = "#d62728"   # red
RESOLVED_COLOR = "#1f77b4"   # blue

DIRECT_CODES   = set(range(271, 275)) | set(range(281, 285))
RESOLVED_CODES = set(range(111, 117)) | set(range(121, 125))


def find_root(sample_dir, base):
    for pat in ("dijets_*.root", "alljets_*.root"):
        hits = sorted(glob.glob(str(Path(base) / sample_dir / pat)))
        if hits:
            return hits[0]
    return None


def find_event_root(sample_short, data_dir):
    hits = sorted(glob.glob(f"{data_dir}/{sample_short}/*.root"))
    return hits[0] if hits else None


def build_event2pc(event_root_path):
    """Global eventID -> processCode dict, summed over all event-level
    category trees (QQ_Events, GG_Events, GQ_Events)."""
    out = {}
    with uproot.open(event_root_path) as f:
        for cat in ("QQ_Events", "GG_Events", "GQ_Events"):
            key = f"{cat}/{cat}"
            if key not in f:
                continue
            t = f[key]
            for e, p in zip(t["eventID"].array(library="np"),
                            t["processCode"].array(library="np")):
                out[int(e)] = int(p)
    return out


def load_jets_tagged(jet_root_path, ev2pc):
    """Return (psi, eta, tag) for every jet across QQ/GG/GQ trees.

    psi  -- (N_jets, 10)
    eta  -- (N_jets,)
    tag  -- (N_jets,) string array, 'direct' or 'resolved' or 'other'
    """
    psi_all, eta_all, tag_all = [], [], []
    with uproot.open(jet_root_path) as f:
        for cat in ("QQ_Events", "GG_Events", "GQ_Events"):
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            ev_ids = t["eventID"].array(library="np")          # per-row scalar
            n_jets = ak.to_numpy(ak.num(t["jet_eta"].array(), axis=1))
            psi_flat = ak.to_numpy(ak.flatten(
                t["jet_psi_curve_flat"].array(), axis=1))
            eta = ak.to_numpy(ak.flatten(t["jet_eta"].array(), axis=1))
            if psi_flat.size != 10 * eta.size:
                raise RuntimeError(f"psi/eta size mismatch in {key}")
            psi = psi_flat.reshape(-1, 10)
            # Look up processCode per event row, then broadcast across jets.
            pc_per_event = np.array([ev2pc.get(int(e), -1) for e in ev_ids],
                                    dtype=np.int32)
            pc_per_jet = np.repeat(pc_per_event, n_jets)
            tags = np.full(pc_per_jet.shape, "other", dtype=object)
            tags[np.isin(pc_per_jet, list(DIRECT_CODES))]   = "direct"
            tags[np.isin(pc_per_jet, list(RESOLVED_CODES))] = "resolved"
            psi_all.append(psi); eta_all.append(eta); tag_all.append(tags)
    if not psi_all:
        return None
    return (np.vstack(psi_all),
            np.concatenate(eta_all),
            np.concatenate(tag_all))


def mean_psi(psi, eta, tags, target_tag, lo, hi, min_jets=20):
    m = (eta >= lo) & (eta < hi) & (tags == target_tag)
    if m.sum() < min_jets:
        return None, m.sum()
    return psi[m].mean(axis=0), m.sum()


def psi_to_rho(psi_mean):
    """Convert one Psi(r) curve into a rho(r) curve.

    rho_i = (Psi(r_i) - Psi(r_{i-1})) / Delta r , for i = 2..10.
    The first bin is dropped so the lowest r reported is r = 0.15.
    """
    return np.diff(psi_mean) / DELTA_R     # shape (9,)


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


LINESTYLE = {64: ":", 105: "--", 141: "-.", 300: "-"}


def _combined_legend():
    """Single legend listing colour (Direct/Resolved) and line style (sqrt(s))."""
    return [
        plt.Line2D([0], [0], color=DIRECT_COLOR,   linewidth=3, label="Direct"),
        plt.Line2D([0], [0], color=RESOLVED_COLOR, linewidth=3, label="Resolved"),
        plt.Line2D([0], [0], color="black", linestyle=":",  linewidth=2, label="64 GeV"),
        plt.Line2D([0], [0], color="black", linestyle="--", linewidth=2, label="105 GeV"),
        plt.Line2D([0], [0], color="black", linestyle="-.", linewidth=2, label="141 GeV"),
        plt.Line2D([0], [0], color="black", linestyle="-",  linewidth=2, label="300 GeV"),
    ]


def create_psi_plot(direct_curves, resolved_curves, eta_label, file_path):
    fig, ax = plt.subplots(figsize=(8, 8))
    r_smooth = np.linspace(R_GRID.min(), R_GRID.max(), 300)
    def smooth(y): return make_interp_spline(R_GRID, y, k=3)(r_smooth)

    for sqrts in (64, 105, 141, 300):
        y = direct_curves.get(sqrts)
        if y is None: continue
        ax.plot(r_smooth, smooth(y), color=DIRECT_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)
    for sqrts in (64, 105, 141, 300):
        y = resolved_curves.get(sqrts)
        if y is None: continue
        ax.plot(r_smooth, smooth(y), color=RESOLVED_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)

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

    ax.legend(handles=_combined_legend(), loc="lower right",
              bbox_to_anchor=(0.97, 0.03), fontsize=18, ncol=1)
    ax.text(0.04, 0.92, eta_label, transform=ax.transAxes, fontsize=24)

    plt.tight_layout()
    plt.savefig(file_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  -> {file_path}")


def create_rho_plot(direct_curves, resolved_curves, eta_label, file_path):
    fig, ax = plt.subplots(figsize=(8, 8))
    r_smooth = np.linspace(RHO_R.min(), RHO_R.max(), 300)
    def smooth(y): return make_interp_spline(RHO_R, y, k=3)(r_smooth)

    for sqrts in (64, 105, 141, 300):
        y = direct_curves.get(sqrts)
        if y is None: continue
        ax.plot(r_smooth, smooth(psi_to_rho(y)), color=DIRECT_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)
    for sqrts in (64, 105, 141, 300):
        y = resolved_curves.get(sqrts)
        if y is None: continue
        ax.plot(r_smooth, smooth(psi_to_rho(y)), color=RESOLVED_COLOR,
                linestyle=LINESTYLE[sqrts], linewidth=2)

    ax.set_xlim(RHO_R.min(), 1.0)   # rho starts at r = 0.15
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$r$", fontsize=28)
    ax.set_ylabel(r"$\rho(r)$", fontsize=28, labelpad=10)
    ax.xaxis.set_label_coords(0.5, -0.08)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    ax.tick_params(axis="both", which="major", labelsize=20)
    ax.tick_params(axis="both", which="minor", labelsize=10)
    ax.minorticks_on()

    ax.legend(handles=_combined_legend(), loc="upper right",
              bbox_to_anchor=(0.97, 0.97), fontsize=18, ncol=1)
    ax.text(0.40, 0.88, eta_label, transform=ax.transAxes, fontsize=24)

    plt.tight_layout()
    plt.savefig(file_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  -> {file_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data",      default="data")
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--out-dir",
                    default=str(Path(__file__).resolve().parent
                                / "output/fig_psi_rho_direct_resolved"))
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_style(use_tex=not args.no_tex)

    print("Loading samples (jet ROOT + event-level processCode lookup):")
    sample_data = {}
    for jdir, ev_short, sqrts in SAMPLES:
        jet_root = find_root(jdir, args.data_jets)
        ev_root  = find_event_root(ev_short, args.data)
        if jet_root is None or ev_root is None:
            print(f"  [skip] {jdir}"); continue
        ev2pc = build_event2pc(ev_root)
        result = load_jets_tagged(jet_root, ev2pc)
        if result is None:
            print(f"  [skip] {jdir}: no jet trees"); continue
        psi, eta, tags = result
        n_dir = (tags == "direct").sum()
        n_res = (tags == "resolved").sum()
        sample_data[sqrts] = (psi, eta, tags)
        print(f"  [load] {jdir} sqrt(s)={sqrts}: "
              f"N_direct={n_dir:,} N_resolved={n_res:,}")

    if not sample_data:
        print("nothing to plot"); sys.exit(1)

    print(f"\nWriting per-eta-bin Psi(r) + rho(r) plots to {out_dir.resolve()}/")
    for lo, hi, label, slug in ETA_BINS:
        direct_psi, resolved_psi = {}, {}
        for sqrts, (psi, eta, tags) in sample_data.items():
            dp, n_d = mean_psi(psi, eta, tags, "direct",   lo, hi)
            rp, n_r = mean_psi(psi, eta, tags, "resolved", lo, hi)
            if dp is not None: direct_psi[sqrts]   = dp
            if rp is not None: resolved_psi[sqrts] = rp
        create_psi_plot(direct_psi, resolved_psi, label,
                        str(out_dir / f"psi_directresolved_{slug}.pdf"))
        create_rho_plot(direct_psi, resolved_psi, label,
                        str(out_dir / f"rho_directresolved_{slug}.pdf"))

    print("\nAll panels written.")


if __name__ == "__main__":
    main()
