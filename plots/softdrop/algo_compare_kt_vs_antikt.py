#!/usr/bin/env python3
"""Side-by-side comparison of kT vs anti-kT jet reconstruction.

Reads two jetreco_softdrop output ROOT files (one produced with kT, one
with anti-kT — same input events, same selection) and tabulates the
key observables per eta bin for each labelled subprocess:

  - N jets that pass E_T > 17 GeV, |eta| < 4
  - <Psi(r=0.3)>
  - <n_subjets>
  - <n_SD>    (baseline (z_cut, beta) = (0.1, 0.0))

and prints the kT - anti-kT difference in each cell. Also writes a
two-panel PDF overlaying the anti-kT and kT curves for <Psi(0.3)> and
<n_subjets> on the paper's binning.

Usage:
    python algo_compare_kt_vs_antikt.py <antikt_root> <kt_root>
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
PAPER_BINS = [(-1.0, 0.0), (0.0, 1.0), (1.0, 1.5), (1.5, 2.0)]


def load(root_path):
    data = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            branches = {b: t[b].array(library="np") for b in
                        ("jet_nsd", "jet_nsubjets", "jet_psi03", "jet_eta")}
            data[cat] = {
                k: np.concatenate([np.asarray(a) for a in v if len(a) > 0])
                for k, v in branches.items()
            }
    return data


def per_bin(data, cat, obs, bins, min_jets=20):
    eta = data[cat]["jet_eta"]
    x = data[cat][obs]
    out = []
    for lo, hi in bins:
        m = (eta >= lo) & (eta < hi)
        n = int(m.sum())
        if n < min_jets:
            out.append((n, np.nan, np.nan))
        else:
            out.append((n, x[m].mean(), x[m].std(ddof=1) / np.sqrt(n)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("antikt_root")
    ap.add_argument("kt_root")
    args = ap.parse_args()

    d_ak = load(args.antikt_root)
    d_kt = load(args.kt_root)
    if not d_ak or not d_kt:
        print("Missing data"); sys.exit(1)

    log_lines = []
    def log(s=""):
        print(s); log_lines.append(s)

    log("=" * 102)
    log("kT vs anti-kT comparison")
    log(f"  anti-kT : {args.antikt_root}")
    log(f"  kT      : {args.kt_root}")
    log("=" * 102)

    # Inclusive sample size per category
    log("\n--- Inclusive jet counts (|eta| < 4, ET > 17 GeV) ---")
    log(f"{'Category':<12}{'anti-kT N':>12}{'kT N':>12}{'kT - anti-kT':>16}{'% diff':>10}")
    for cat in CATEGORIES:
        if cat not in d_ak or cat not in d_kt:
            continue
        n_a = d_ak[cat]["jet_nsd"].size
        n_k = d_kt[cat]["jet_nsd"].size
        pct = 100 * (n_k - n_a) / n_a if n_a else 0.0
        log(f"{LABELS[cat]:<12}{n_a:>12,}{n_k:>12,}{n_k - n_a:>+16,}{pct:>+9.2f}%")

    for obs, display in [("jet_psi03",    "<Psi(r=0.3)>"),
                         ("jet_nsubjets", "<n_subjets>"),
                         ("jet_nsd",      "<n_SD>")]:
        log("\n" + "-" * 102)
        log(f"{display}  per eta bin, per category")
        log("-" * 102)
        for cat in CATEGORIES:
            if cat not in d_ak or cat not in d_kt:
                continue
            log(f"\n  {LABELS[cat]}")
            log(f"  {'eta bin':<16}{'ak N':>8}{'ak <X>':>12}{'kt N':>8}{'kt <X>':>12}"
                f"{'kt - ak':>12}{'% diff':>9}")
            ak = per_bin(d_ak, cat, obs, PAPER_BINS)
            kt = per_bin(d_kt, cat, obs, PAPER_BINS)
            for (lo, hi), (na, ma, ea), (nk, mk, ek) in zip(PAPER_BINS, ak, kt):
                bin_str = f"[{lo:+.1f},{hi:+.1f})"
                if np.isnan(ma) or np.isnan(mk):
                    log(f"  {bin_str:<16}{na:>8}{'—':>12}{nk:>8}{'—':>12}{'—':>12}{'—':>9}")
                else:
                    diff = mk - ma
                    pct = 100 * diff / ma if ma != 0 else 0
                    log(f"  {bin_str:<16}{na:>8}{ma:>12.4f}{nk:>8}{mk:>12.4f}"
                        f"{diff:>+12.4f}{pct:>+8.2f}%")

    # --- plots ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, obs, ylabel in zip(axes,
                               ("jet_psi03", "jet_nsubjets"),
                               (r"$\langle \Psi(r=0.3)\rangle$",
                                r"$\langle n_{\mathrm{subjets}}\rangle$")):
        centers = [0.5*(lo+hi) for lo, hi in PAPER_BINS]
        for cat in ("QQ_Events", "GG_Events", "GQ_Events"):
            if cat not in d_ak or cat not in d_kt:
                continue
            ak = per_bin(d_ak, cat, obs, PAPER_BINS)
            kt = per_bin(d_kt, cat, obs, PAPER_BINS)
            ak_m = [v[1] for v in ak];  ak_e = [v[2] for v in ak]
            kt_m = [v[1] for v in kt];  kt_e = [v[2] for v in kt]
            ax.errorbar(centers, ak_m, yerr=ak_e, fmt="o-",
                        color=COLORS[cat], linewidth=2, markersize=8,
                        label=f"{LABELS[cat]} anti-kT")
            ax.errorbar(centers, kt_m, yerr=kt_e, fmt="s--",
                        color=COLORS[cat], linewidth=2, markersize=8,
                        markerfacecolor="white",
                        label=f"{LABELS[cat]} kT")
        ax.set_xlabel(r"$\eta$ (bin centre)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(frameon=False, fontsize=9, ncol=2)
    fig.suptitle("anti-kT (solid circles) vs kT (open squares)", fontsize=13, y=1.00)
    fig.tight_layout()
    fig.savefig("algo_compare_kt_vs_antikt.pdf"); plt.close(fig)
    log(f"\nWrote algo_compare_kt_vs_antikt.pdf")

    out_log = Path("algo_compare_kt_vs_antikt.log")
    out_log.write_text("\n".join(log_lines) + "\n")
    print(f"Wrote {out_log.resolve()}")


if __name__ == "__main__":
    main()
