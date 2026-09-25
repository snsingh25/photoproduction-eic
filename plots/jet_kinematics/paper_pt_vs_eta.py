#!/usr/bin/env python3
"""Jet pT vs eta density for the paper-config dijet samples.

Shows where the accepted jets actually sit in the (eta, pT) plane, one
square panel per centre-of-mass energy. The overlaid line is the mean pT
in each eta slice, which is the quantity that drives the scale
contamination of n_subjets discussed in the ROC paper: n_subjets counts
emissions above an absolute pT^2 threshold, so wherever <pT> moves with
eta, the observable moves with it.

QQ and GG jets are combined -- the question here is where jets land, not
how the two classes differ.

Usage:
    python plots/jet_kinematics/paper_pt_vs_eta.py
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "style"))
import paper_style as ps  # noqa: E402

SAMPLES = [
    (64,  "eic64_antikt_dijets"),
    (105, "eic105_antikt_dijets"),
    (141, "eic141_antikt_dijets"),
    (300, "hera300_antikt_dijets"),
]

ETA_RANGE = (-2.0, 4.0)
PT_RANGE = (5.0, 32.0)
NBINS = (90, 68)


def find_root(data_jets, sample):
    hits = sorted(glob.glob(f"{data_jets}/{sample}/dijets_*.root"))
    return hits[0] if hits else None


def load_pt_eta(root_path):
    """Flat per-jet pT and eta, QQ and GG combined."""
    pt, eta = [], []
    with uproot.open(root_path) as f:
        for cat in ("QQ_Events", "GG_Events"):
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            px = np.concatenate([np.asarray(a) for a in
                                 t["jet_px"].array(library="np") if len(a)])
            py = np.concatenate([np.asarray(a) for a in
                                 t["jet_py"].array(library="np") if len(a)])
            e = np.concatenate([np.asarray(a) for a in
                                t["jet_eta"].array(library="np") if len(a)])
            pt.append(np.hypot(px, py))
            eta.append(e)
    if not pt:
        return None, None
    return np.concatenate(pt), np.concatenate(eta)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--out-dir",
                    default=str(HERE / "output/paper_pt_vs_eta"))
    args = ap.parse_args()

    data_jets = Path(args.data_jets).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ps.use()
    cmap = ps.density_cmap()

    fig, axes = plt.subplots(2, 2, figsize=(ps.FULL_W, ps.FULL_W * 1.02))
    xedges = np.linspace(*ETA_RANGE, NBINS[0] + 1)
    yedges = np.linspace(*PT_RANGE, NBINS[1] + 1)

    meshes = []
    summary = ["Where the jets sit: peak of the (eta, pT) density and the",
               "spread of <pT> across eta. anti-kT R=1.0, dijets, ET>10/7 GeV.",
               "=" * 72,
               f"{'sqrt(s)':>8}{'jets':>9}{'peak eta':>10}{'peak pT':>9}"
               f"{'<pT> min':>10}{'<pT> max':>10}{'spread':>9}"]
    for idx, (sqrts, sample) in enumerate(SAMPLES):
        ax = axes[idx // 2][idx % 2]
        ps.square(ax)
        rf = find_root(data_jets, sample)
        if rf is None:
            print(f"  [skip] {sample}")
            ax.set_visible(False)
            continue
        pt, eta = load_pt_eta(rf)
        print(f"  [load] {sqrts:>3} GeV  {len(pt)} jets  "
              f"pT {pt.min():.1f}-{pt.max():.1f}  eta {eta.min():.1f}-{eta.max():.1f}")

        H, _, _ = np.histogram2d(eta, pt, bins=[xedges, yedges])
        H = np.ma.masked_where(H == 0, H)
        m = ax.pcolormesh(xedges, yedges, H.T, cmap=cmap,
                          norm=LogNorm(vmin=1, vmax=H.max()),
                          rasterized=True, shading="flat")
        meshes.append(m)

        # mean pT per eta slice: the trend that drives the n_subjets effect
        centres = 0.5 * (xedges[:-1] + xedges[1:])
        mean_pt = np.array([
            pt[(eta >= lo) & (eta < hi)].mean()
            if ((eta >= lo) & (eta < hi)).sum() > 30 else np.nan
            for lo, hi in zip(xedges[:-1], xedges[1:])])
        ax.plot(centres, mean_pt, color=ps.OBS_COLOR["jet_nsd"],
                linewidth=1.6, zorder=5, label=r"$\langle p_T \rangle$")

        iy, ix = np.unravel_index(np.argmax(H.filled(0)), H.shape)
        ycent = 0.5 * (yedges[:-1] + yedges[1:])
        summary.append(
            f"{sqrts:>8}{len(pt):>9}{centres[iy]:>10.2f}{ycent[ix]:>9.1f}"
            f"{np.nanmin(mean_pt):>10.1f}{np.nanmax(mean_pt):>10.1f}"
            f"{np.nanmax(mean_pt)-np.nanmin(mean_pt):>9.1f}")

        ax.set_xlim(*ETA_RANGE)
        ax.set_ylim(*PT_RANGE)
        ax.set_title(ps.energy_label(sqrts), pad=6)
        ax.grid(False)
        if idx // 2 == 1:
            ax.set_xlabel(r"jet $\eta$")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel(r"jet $p_T$ [GeV]")
        else:
            ax.set_yticklabels([])
        if idx == 0:
            ax.legend(loc="upper right", labelcolor=ps.INK)

    fig.subplots_adjust(wspace=0.06, hspace=0.14, right=0.88)
    cax = fig.add_axes([0.90, 0.13, 0.018, 0.74])
    cb = fig.colorbar(meshes[0], cax=cax)
    cb.set_label("jets / bin")
    cb.outline.set_linewidth(0.6)
    cb.outline.set_edgecolor(ps.AXIS)

    out = out_dir / "fig_pt_vs_eta.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"Wrote {out}")

    log_path = out_dir / "pt_vs_eta_summary.log"
    log_path.write_text("\n".join(summary) + "\n")
    print()
    print("\n".join(summary))


if __name__ == "__main__":
    main()
