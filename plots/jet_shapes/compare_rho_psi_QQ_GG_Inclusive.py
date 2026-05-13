#!/usr/bin/env python3
"""Combined comparison: Psi(r) and rho(r) with three classifications
overlaid in one panel -- QQ, GG, and Direct+Resolved (Dir+Res, the
inclusive photoproduction sample) -- across all four sqrt(s) values.

Sibling of compare_rho_psi_all.py but with Direct and Resolved
merged into one inclusive curve so the QQ/GG/inclusive structure
is the headline.

For each of the 4 eta bins, produces:
    psi_qqggincl_<etabin>.pdf   (Psi(r), 12 curves: 3 classes x 4 energies)
    rho_qqggincl_<etabin>.pdf   (rho(r), same layout)

Class -> colour:
    QQ       red       (#c0392b)
    GG       blue      (#2980b9)
    Dir+Res  black     (#000000)

sqrt(s) -> line style:  64 ":"  /  105 "--"  /  141 "-."  /  300 "-"

Inputs:
    data-jets/<sample>_antikt_dijets/   (HERA also at anti-kT EtMin=10
                                         so all four samples are
                                         compared apples-to-apples)
    data/<sample>/*.root                (event-level processCode lookup)

Output: plots/jet_shapes/output/compare_rho_psi_QQ_GG_Inclusive/
"""

import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import make_interp_spline


R_GRID  = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
DELTA_R = 0.1
RHO_R   = R_GRID                       # paper convention: report rho at the UPPER edge of each annulus (r = 0.1, 0.2, ..., 1.0)

SAMPLES = [
    ("eic64_antikt_dijets",   "eic64",   64),
    ("eic105_antikt_dijets",  "eic105",  105),
    ("eic141_antikt_dijets",  "eic141",  141),
    ("hera300_antikt_dijets", "hera300", 300),
]

ETA_BINS = [
    (-1.0, 0.0, r"$-1 < \eta < 0$",  "minus1to0"),
    ( 0.0, 1.0, r"$0 < \eta < 1$",   "0to1"),
    ( 1.0, 1.5, r"$1 < \eta < 1.5$", "1to1p5"),
    ( 1.5, 2.0, r"$1.5 < \eta < 2$", "1p5to2"),
]

CLASS_COLOUR = {
    "QQ":       "#c0392b",   # red
    "GG":       "#2980b9",   # blue
    "Dir+Res":  "#000000",   # black (inclusive)
}
LINESTYLE = {64: ":", 105: "--", 141: "-.", 300: "-"}

DIRECT_CODES   = set(range(271, 275)) | set(range(281, 285))
RESOLVED_CODES = set(range(111, 117)) | set(range(121, 125))
INCLUSIVE_CODES = DIRECT_CODES | RESOLVED_CODES

REPO = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def find_dijet_root(sample_dir):
    base = REPO / "data-jets" / sample_dir
    hits = sorted(glob.glob(str(base / "dijets_*.root")))
    return hits[0] if hits else None


def find_event_root(sample_short):
    hits = sorted(glob.glob(str(REPO / "data" / sample_short / "*.root")))
    return hits[0] if hits else None


def build_event2pc(event_root):
    out = {}
    with uproot.open(event_root) as f:
        for cat in ("QQ_Events", "GG_Events", "GQ_Events"):
            key = f"{cat}/{cat}"
            if key not in f:
                continue
            t = f[key]
            for e, p in zip(t["eventID"].array(library="np"),
                            t["processCode"].array(library="np")):
                out[int(e)] = int(p)
    return out


def load_sample(jet_root, ev2pc):
    """Return dict with keys 'QQ', 'GG', 'Dir+Res' -> (psi: (N, 10), eta: (N,))."""
    out = {"QQ": ([], []), "GG": ([], []), "Dir+Res": ([], [])}

    with uproot.open(jet_root) as f:
        for cat in ("QQ_Events", "GG_Events", "GQ_Events"):
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            ev_ids   = t["eventID"].array(library="np")
            n_jets   = ak.to_numpy(ak.num(t["jet_eta"].array(), axis=1))
            psi_flat = ak.to_numpy(ak.flatten(
                t["jet_psi_curve_flat"].array(), axis=1))
            eta      = ak.to_numpy(ak.flatten(t["jet_eta"].array(), axis=1))
            if psi_flat.size == 0:
                continue
            psi = psi_flat.reshape(-1, 10)

            # Parton-flavour class
            if cat == "QQ_Events":
                out["QQ"][0].append(psi); out["QQ"][1].append(eta)
            elif cat == "GG_Events":
                out["GG"][0].append(psi); out["GG"][1].append(eta)

            # Inclusive photon-mode class (Direct OR Resolved)
            pc_per_event = np.array([ev2pc.get(int(e), -1) for e in ev_ids],
                                    dtype=np.int32)
            pc_per_jet = np.repeat(pc_per_event, n_jets)
            incl_mask = np.isin(pc_per_jet, list(INCLUSIVE_CODES))
            if incl_mask.any():
                out["Dir+Res"][0].append(psi[incl_mask])
                out["Dir+Res"][1].append(eta[incl_mask])

    result = {}
    for k, (psis, etas) in out.items():
        if not psis:
            result[k] = (np.empty((0, 10)), np.empty((0,)))
        else:
            result[k] = (np.vstack(psis), np.concatenate(etas))
    return result


def mean_psi_in_eta(psi, eta, lo, hi, min_jets=10):
    """Lowered min_jets from 20 to 10 so the low-stats EIC samples
    (EIC 64 has only ~10 GG jets at -1<eta<0, EIC 105 has ~16)
    still produce a curve, albeit a noisier one."""
    m = (eta >= lo) & (eta < hi)
    if m.sum() < min_jets:
        return None
    return psi[m].mean(axis=0)


def psi_to_rho(psi_curve):
    return np.diff(np.concatenate(([0.0], psi_curve))) / DELTA_R


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def setup_style():
    plt.style.use("default")
    plt.rcParams.update({
        "font.size": 18, "font.family": "serif",
        "font.serif": ["Computer Modern"], "text.usetex": True,
        "axes.linewidth": 1.2, "axes.spines.top": False,
        "axes.spines.right": False, "lines.linewidth": 1.8,
        "legend.frameon": False, "figure.facecolor": "white",
    })


def _legend_handles():
    """Layout with ncol=2 fills column-major. To get
        col 1 = {QQ, GG, Dir+Res, [empty]}
        col 2 = {64, 105, 141, 300 GeV}
    we pad the class column with a transparent spacer so each
    column holds 4 entries."""
    handles = []
    for label, col in CLASS_COLOUR.items():
        handles.append(plt.Line2D([0], [0], color=col, linewidth=3, label=label))
    handles.append(plt.Line2D([0], [0], color="none", label=" "))   # spacer
    for sqrts in (64, 105, 141, 300):
        handles.append(plt.Line2D([0], [0], color="black",
                                  linestyle=LINESTYLE[sqrts],
                                  linewidth=2, label=f"{sqrts} GeV"))
    return handles


def plot_panel(curves_by_class_sqrts, eta_label, out_pdf, *, kind="psi"):
    fig, ax = plt.subplots(figsize=(9, 8))
    r_smooth = np.linspace(RHO_R.min(), RHO_R.max(), 300) if kind == "rho" \
               else np.linspace(R_GRID.min(), R_GRID.max(), 300)

    for cls in ("QQ", "GG", "Dir+Res"):
        per_sqrts = curves_by_class_sqrts.get(cls, {})
        col = CLASS_COLOUR[cls]
        for sqrts in (64, 105, 141, 300):
            psi = per_sqrts.get(sqrts)
            if psi is None:
                continue
            if kind == "rho":
                y = psi_to_rho(psi)
                x_pts = RHO_R
            else:
                y = psi
                x_pts = R_GRID
            sp = make_interp_spline(x_pts, y, k=3)(r_smooth)
            ax.plot(r_smooth, sp, color=col,
                    linestyle=LINESTYLE[sqrts], linewidth=1.8)

    if kind == "psi":
        ax.set_xlim(0.05, 1.0)
        ax.set_yscale("log")
        ax.set_yticks([0.05, 0.1, 1.0])
        ax.set_yticklabels(["0.05", "0.1", "1.0"])
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.get_major_formatter().set_scientific(False)
        ax.set_ylabel(r"$\Psi(r)$", fontsize=22)
    else:
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(bottom=0.0)
        ax.set_ylabel(r"$\rho(r)$", fontsize=22)

    ax.set_xlabel(r"$r$", fontsize=22)
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.minorticks_on()

    eta_x = 0.55 if kind == "rho" else 0.05
    eta_y = 0.55 if kind == "rho" else 0.10
    ax.text(eta_x, eta_y, eta_label, transform=ax.transAxes,
            fontsize=20, va="center")
    ax.legend(handles=_legend_handles(),
              loc="upper right" if kind == "rho" else "lower right",
              fontsize=12, ncol=2)
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  -> {out_pdf}")


def main():
    setup_style()
    out_dir = REPO / "plots" / "jet_shapes" / "output" / "compare_rho_psi_QQ_GG_Inclusive"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading samples (jet ROOT + processCode lookup):")
    sample_data = {}
    for jdir, ev_short, sqrts in SAMPLES:
        jet_root = find_dijet_root(jdir)
        ev_root  = find_event_root(ev_short)
        if jet_root is None or ev_root is None:
            print(f"  [skip] {jdir}: missing ROOT"); continue
        ev2pc = build_event2pc(ev_root)
        sample_data[sqrts] = load_sample(jet_root, ev2pc)
        sizes = {k: v[0].shape[0] for k, v in sample_data[sqrts].items()}
        print(f"  [load] sqrt(s)={sqrts}: {sizes}")

    if not sample_data:
        print("nothing to plot"); sys.exit(1)

    print(f"\nWriting plots to {out_dir.resolve()}/")
    for lo, hi, label, slug in ETA_BINS:
        curves = {cls: {} for cls in ("QQ", "GG", "Dir+Res")}
        for sqrts, by_class in sample_data.items():
            for cls, (psi, eta) in by_class.items():
                m_psi = mean_psi_in_eta(psi, eta, lo, hi)
                if m_psi is not None:
                    curves[cls][sqrts] = m_psi
        plot_panel(curves, label, str(out_dir / f"psi_qqggincl_{slug}.pdf"),
                   kind="psi")
        plot_panel(curves, label, str(out_dir / f"rho_qqggincl_{slug}.pdf"),
                   kind="rho")

    print("\nAll panels written.")


if __name__ == "__main__":
    main()
