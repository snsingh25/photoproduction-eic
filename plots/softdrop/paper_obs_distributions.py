#!/usr/bin/env python3
"""QQ vs GG distributions of each observable, per pseudorapidity bin.

One panel per eta bin plus a fifth panel pooling the whole analysis window,
so the per-bin shapes and the pooled shape can be read side by side. This is
the direct picture behind Sec. V: n_subjets shifts bodily to lower values as
eta grows, because its effective threshold rises as cosh(eta), and pooling
therefore stacks several differently-shifted distributions on top of each
other. n_SD and Psi(0.3) do not move.

Usage:
    python plots/softdrop/paper_obs_distributions.py
    python plots/softdrop/paper_obs_distributions.py --sample eic141_antikt_dijets
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "style"))
import paper_style as ps  # noqa: E402

import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location(
    "roc_cross_energy_with_psi", HERE / "roc_cross_energy_with_psi.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

ETA_BINS = [(-1.0, 0.0), (0.0, 1.0), (1.0, 1.5), (1.5, 2.0)]
WINDOW = (-1.0, 2.0)

# QQ / GG are two classes, not an ordered scale: two fixed categorical slots.
# Red for quarks and blue for gluons follows the usual convention.
CLS = {"QQ_Events": ("QQ", "#e34948", "///"),
       "GG_Events": ("GG", "#2a78d6", "\\\\\\")}

# observable -> (branch, axis label, integer bin edges or None for continuous)
OBS = {
    "nsubjets": ("jet_nsubjets", r"$n_{\mathrm{sbj}}$",
                 np.arange(-0.5, 11.5, 1.0)),
    "nsd":      ("jet_nsd", r"$n_{\mathrm{SD}}$",
                 np.arange(-0.5, 11.5, 1.0)),
    "psi03":    ("jet_psi03", r"$\Psi(r=0.3)$",
                 np.linspace(0.0, 1.0, 26)),
}

SQRTS = {"eic64_antikt_dijets": 64, "eic105_antikt_dijets": 105,
         "eic141_antikt_dijets": 141, "hera300_antikt_dijets": 300}


def panel(ax, d, branch, edges, title, sqrts, show_legend):
    hts = []
    for cat, (lbl, col, hatch) in CLS.items():
        if cat not in d:
            continue
        v = d[cat][branch]
        h, _ = np.histogram(v, bins=edges, density=True)
        hts.append(h)
        ax.hist(v, bins=edges, density=True, histtype="stepfilled",
                facecolor=col, alpha=0.25, hatch=hatch, edgecolor=col,
                linewidth=0.0, label=lbl)
        ax.hist(v, bins=edges, density=True, histtype="step",
                edgecolor=col, linewidth=1.3)
    # headroom so the corner labels never sit on top of the histogram
    top = max(h.max() for h in hts) if hts else 1.0
    ax.set_ylim(0, top * 1.42)
    ax.text(0.04, 0.95, title, transform=ax.transAxes,
            va="top", ha="left", fontsize=7.5)
    ax.text(0.04, 0.845, rf"$ep$, {sqrts}~GeV", transform=ax.transAxes,
            va="top", ha="left", fontsize=7.5, color=ps.MUTED)
    if show_legend:
        ax.legend(loc="upper right", fontsize=7, borderpad=0.35)
    ax.grid(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--sample", default="hera300_antikt_dijets")
    ap.add_argument("--out-dir",
                    default=str(HERE / "output/obs_distributions"))
    args = ap.parse_args()

    data_jets = Path(args.data_jets).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    sqrts = SQRTS.get(args.sample, "?")

    rf = rc.find_root(data_jets, args.sample)
    if rf is None:
        print(f"no ROOT file for {args.sample}")
        sys.exit(1)
    d = rc.load(rf)
    print(f"[load] {Path(rf).name}")

    ps.use()
    slices = [(rf"${lo:g} < \eta_{{\mathrm{{jet}}}} < {hi:g}$",
               rc.slice_eta(d, lo, hi)) for lo, hi in ETA_BINS]
    slices.append((rf"${WINDOW[0]:g} < \eta_{{\mathrm{{jet}}}} < {WINDOW[1]:g}$"
                   r"  (all)", rc.slice_eta(d, *WINDOW)))

    for key, (branch, xlabel, edges) in OBS.items():
        fig, axes = plt.subplots(2, 3, figsize=(ps.FULL_W, ps.FULL_W * 0.62))
        flat = [axes[r][c] for r in range(2) for c in range(3)]
        for i, (title, sub) in enumerate(slices):
            ax = flat[i]
            if "QQ_Events" not in sub or "GG_Events" not in sub:
                ax.set_visible(False)
                continue
            panel(ax, sub, branch, edges, title, sqrts, show_legend=(i == 0))
            if i >= 2:
                ax.set_xlabel(xlabel)
            if i % 3 == 0:
                ax.set_ylabel(r"$\frac{1}{N_{\mathrm{jet}}}\,"
                              r"\frac{dN_{\mathrm{jet}}}{d" +
                              xlabel.strip("$") + r"}$", fontsize=8)
        flat[5].set_visible(False)
        flat[2].set_xlabel(xlabel)
        fig.subplots_adjust(wspace=0.26, hspace=0.22)
        out = out_dir / f"dist_{key}_{args.sample}.pdf"
        fig.savefig(out)
        plt.close(fig)
        print(f"Wrote {out}")

    # numeric companion: where each distribution sits, bin by bin
    lines = [f"Mean and width per eta bin — {args.sample}", "=" * 72,
             f"{'bin':>16}{'N_QQ':>9}{'N_GG':>8}"
             f"{'<QQ>':>9}{'<GG>':>9}{'shift':>8}"]
    for key, (branch, _x, _e) in OBS.items():
        lines.append(f"--- {branch} ---")
        for title, sub in slices:
            if "QQ_Events" not in sub or "GG_Events" not in sub:
                continue
            q, g = sub["QQ_Events"][branch], sub["GG_Events"][branch]
            t = title.replace(r"\eta_{\mathrm{jet}}", "eta").replace("$", "")
            lines.append(f"{t:>16}{q.size:>9}{g.size:>8}"
                         f"{q.mean():>9.3f}{g.mean():>9.3f}"
                         f"{g.mean()-q.mean():>8.3f}")
    log = out_dir / f"dist_summary_{args.sample}.log"
    log.write_text("\n".join(lines) + "\n")
    print()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
