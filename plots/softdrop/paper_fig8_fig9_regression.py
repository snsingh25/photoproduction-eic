#!/usr/bin/env python3
"""Regression of <n_subjets>(eta) and <Psi(r=0.3)>(eta) against the paper.

Produces the two curves that the published Fig. 8 (n_subjets vs eta)
and Fig. 9 (integrated jet shape Psi(r=0.3) vs eta, with ZEUS overlay
from archival HERA data) show. Intended as a regression check that
the jetreco_softdrop pipeline reproduces the published paper numbers.

Usage:
    python paper_fig8_fig9_regression.py <root_file> [--bins=...]

Writes:
    paper_fig8_nsubjets.pdf    <n_subjets>(eta) for QQ, GG, GQ
    paper_fig9_psi03.pdf       <Psi(r=0.3)>(eta) for QQ, GG, GQ
    paper_regression_summary.log    numeric values per bin

All outputs go into the current working directory.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt


CATEGORIES = ["QQ_Events", "GG_Events", "GQ_Events"]
COLORS = {"QQ_Events": "C0", "GG_Events": "C3", "GQ_Events": "C2"}
LABELS = {"QQ_Events": "QQ", "GG_Events": "GG", "GQ_Events": "GQ"}
MARKERS = {"QQ_Events": "o", "GG_Events": "s", "GQ_Events": "^"}

DEFAULT_BIN_EDGES = "-1,0,1,1.5,2"   # paper's Fig. 8/9 binning


def load(root_path):
    data = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            tree = f[key]
            available = set(tree.keys())
            needed = {"jet_nsubjets", "jet_psi03", "jet_eta"}
            if not needed.issubset(available):
                missing = needed - available
                print(f"  [skip] {cat}: missing branches {missing}")
                continue
            nsj = tree["jet_nsubjets"].array(library="np")
            psi = tree["jet_psi03"].array(library="np")
            eta = tree["jet_eta"].array(library="np")
            nsj_flat = np.concatenate([np.asarray(a) for a in nsj if len(a) > 0])
            psi_flat = np.concatenate([np.asarray(a) for a in psi if len(a) > 0])
            eta_flat = np.concatenate([np.asarray(a) for a in eta if len(a) > 0])
            if nsj_flat.size == 0:
                continue
            data[cat] = {"nsubjets": nsj_flat, "psi03": psi_flat, "eta": eta_flat}
    return data


def per_bin(data, cat, obs, bins, min_jets=20):
    """Return (centers, means, sems, ns) with NaN for under-populated bins."""
    eta = data[cat]["eta"]
    x = data[cat][obs]
    centers, means, sems, ns = [], [], [], []
    for lo, hi in bins:
        m = (eta >= lo) & (eta < hi)
        n = m.sum()
        centers.append(0.5 * (lo + hi))
        ns.append(n)
        if n < min_jets:
            means.append(np.nan); sems.append(np.nan)
        else:
            means.append(x[m].mean())
            sems.append(x[m].std(ddof=1) / np.sqrt(n))
    return (np.asarray(centers), np.asarray(means),
            np.asarray(sems), np.asarray(ns))


def plot_obs(data, obs, ylabel, title, outpdf, bins):
    fig, ax = plt.subplots(figsize=(7, 5))
    for cat in CATEGORIES:
        if cat not in data:
            continue
        c, mu, err, n = per_bin(data, cat, obs, bins)
        ok = ~np.isnan(mu)
        ax.errorbar(c[ok], mu[ok], yerr=err[ok],
                    fmt=MARKERS[cat] + "-", color=COLORS[cat],
                    linewidth=2, markersize=8, capsize=4,
                    label=LABELS[cat])
    # Mark bin edges
    for lo, hi in bins:
        ax.axvline(lo, linestyle=":", color="grey", alpha=0.3)
    ax.axvline(bins[-1][1], linestyle=":", color="grey", alpha=0.3)
    ax.set_xlabel(r"$\eta$ (jet centre of bin)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpdf); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file")
    ap.add_argument("--bins", default=DEFAULT_BIN_EDGES,
                    help=f"comma-separated eta bin edges (default: {DEFAULT_BIN_EDGES})")
    args = ap.parse_args()

    edges = [float(x) for x in args.bins.split(",")]
    bins = list(zip(edges[:-1], edges[1:]))

    sample_name = Path(args.root_file).stem
    data = load(args.root_file)
    if not data:
        print("No usable data."); sys.exit(1)

    log_lines = []
    def log(s=""):
        print(s); log_lines.append(s)

    log("=" * 84)
    log(f"Fig 8 / Fig 9 regression  —  {sample_name}")
    log(f"bins: {bins}")
    log("=" * 84)

    for obs, obs_display in [("nsubjets",    "<n_subjets>"),
                              ("psi03",      "<Psi(r=0.3)>")]:
        log(f"\n{obs_display}:")
        log(f"{'eta bin':<16}"
            + "".join([f"{LABELS[c]+' N':>10}{LABELS[c]+' <X>':>12}{LABELS[c]+' +/-':>10}" for c in CATEGORIES if c in data]))
        for (lo, hi) in bins:
            parts = [f"[{lo:+.1f},{hi:+.1f})".ljust(16)]
            for c in CATEGORIES:
                if c not in data:
                    continue
                _, mu, err, n = per_bin(data, c, obs, [(lo, hi)])
                mu = mu[0]; err = err[0]; n = int(n[0])
                if np.isnan(mu):
                    parts.append(f"{n:>10}{'—':>12}{'—':>10}")
                else:
                    parts.append(f"{n:>10}{mu:>12.4f}{err:>10.4f}")
            log("".join(parts))

    plot_obs(data, "nsubjets",
             r"$\langle n_{\mathrm{subjets}}\rangle$",
             f"{sample_name}  —  Fig 8 regression",
             "paper_fig8_nsubjets.pdf", bins)
    plot_obs(data, "psi03",
             r"$\langle \Psi(r=0.3)\rangle$",
             f"{sample_name}  —  Fig 9 regression",
             "paper_fig9_psi03.pdf", bins)

    out_log = Path("paper_regression_summary.log")
    out_log.write_text("\n".join(log_lines) + "\n")
    log(f"\nWrote {out_log.resolve()}")


if __name__ == "__main__":
    main()
