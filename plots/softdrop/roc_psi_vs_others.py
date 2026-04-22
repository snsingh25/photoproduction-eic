#!/usr/bin/env python3
"""ROC for Psi(r=0.3) alongside n_SD and n_subjets on the same jets.

Psi(r=0.3) is the paper's primary substructure discriminant (thin/thick
jet classification). This script treats it as a continuous tagger and
compares it head-to-head with the two multiplicity tagggers used for
the N4 referee-response.

Sign conventions:
    - n_SD, n_subjets : larger for gluon jets (signal = GG, X >= c cut).
    - Psi(r=0.3)      : larger for quark jets (signal = GG, X <= c cut).
      Implemented by running the same ROC code on (-Psi) so a single
      "signal-larger" convention applies to all three observables.

Usage:
    python plots/softdrop/roc_psi_vs_others.py <root_file>
    (default: data-jets/hera300_pT7/alljets_*.root)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt


CATEGORIES = ["QQ_Events", "GG_Events"]


def load(root_path):
    data = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            for br in ("jet_nsd", "jet_nsubjets", "jet_psi03"):
                if br not in t.keys():
                    print(f"  [warn] {cat}: missing {br}")
                    return None
            nsd = t["jet_nsd"].array(library="np")
            nsj = t["jet_nsubjets"].array(library="np")
            psi = t["jet_psi03"].array(library="np")
            data[cat] = {
                "nsd":      np.concatenate([np.asarray(a) for a in nsd if len(a) > 0]),
                "nsubjets": np.concatenate([np.asarray(a) for a in nsj if len(a) > 0]),
                "psi03":    np.concatenate([np.asarray(a) for a in psi if len(a) > 0]),
            }
    return data


def roc_continuous(sig, bkg, n_points=200):
    """Sweep a fine threshold grid for a continuous-valued discriminant.
    Convention: larger X favours signal; cut is X >= c."""
    lo = float(min(sig.min(), bkg.min()))
    hi = float(max(sig.max(), bkg.max()))
    # include one step past the ends so the curve reaches (0,0) and (1,1)
    grid = np.linspace(lo - 1e-6, hi + 1e-6, n_points)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    auc = np.trapezoid(eff_sig[order], eff_bkg[order])
    return eff_sig[order], eff_bkg[order], auc


def roc_integer(sig, bkg):
    vmin = int(min(sig.min(), bkg.min()))
    vmax = int(max(sig.max(), bkg.max()))
    grid = np.arange(vmin, vmax + 2)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    auc = np.trapezoid(eff_sig[order], eff_bkg[order])
    return eff_sig[order], eff_bkg[order], auc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file",
                    default="data-jets/hera300_pT7/alljets_hera300_pT7_R10_EtMin17.root",
                    nargs="?")
    args = ap.parse_args()

    data = load(args.root_file)
    if not data:
        print(f"Missing data in {args.root_file}"); sys.exit(1)

    sample_name = Path(args.root_file).stem
    qq = data["QQ_Events"]
    gg = data["GG_Events"]

    log_lines = []
    def log(s=""):
        print(s); log_lines.append(s)

    log("=" * 88)
    log(f"Psi(r=0.3) vs n_SD vs n_subjets — ROC comparison  —  {sample_name}")
    log("=" * 88)
    log(f"  N_QQ = {qq['nsd'].size:,}   N_GG = {gg['nsd'].size:,}\n")

    # Means and std
    log(f"{'Observable':<18}{'<QQ>':>10}{'sig_QQ':>10}{'<GG>':>10}{'sig_GG':>10}"
        f"{'sep':>10}{'AUC':>10}")
    log("-" * 78)

    results = []
    # n_SD  (gluon-larger, integer)
    es, eb, auc = roc_integer(gg["nsd"], qq["nsd"])
    results.append(("n_SD",         gg["nsd"], qq["nsd"], es, eb, auc, "C1", "integer"))
    # n_subjets  (gluon-larger, integer)
    es, eb, auc = roc_integer(gg["nsubjets"], qq["nsubjets"])
    results.append(("n_subjets",    gg["nsubjets"], qq["nsubjets"], es, eb, auc, "C4", "integer"))
    # Psi(r=0.3)  (QUARK-larger -> flip sign so "signal-larger" convention holds)
    sig_psi = -gg["psi03"]
    bkg_psi = -qq["psi03"]
    es, eb, auc = roc_continuous(sig_psi, bkg_psi, n_points=400)
    results.append((r"Psi(r=0.3)", gg["psi03"], qq["psi03"], es, eb, auc, "C2", "continuous"))

    for name, gg_arr, qq_arr, _es, _eb, auc, _col, _kind in results:
        sep = (gg_arr.mean() - qq_arr.mean()) / np.sqrt(
            gg_arr.var(ddof=1) + qq_arr.var(ddof=1))
        log(f"{name:<18}"
            f"{qq_arr.mean():>10.3f}{qq_arr.std(ddof=1):>10.3f}"
            f"{gg_arr.mean():>10.3f}{gg_arr.std(ddof=1):>10.3f}"
            f"{sep:>+10.3f}{auc:>10.3f}")

    # --- plot -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    for name, _gg, _qq, es, eb, auc, col, _kind in results:
        ax.plot(eb, es, "o-" if _kind == "integer" else "-",
                color=col, linewidth=2, markersize=6 if _kind == "integer" else 0,
                label=f"{name}  AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="grey", alpha=0.5, label="random")
    ax.set_xlabel(r"$\varepsilon$(QQ)")
    ax.set_ylabel(r"$\varepsilon$(GG)")
    ax.set_title(f"{sample_name}  —  ROC (GG signal vs QQ background)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.3)
    ax.text(0.02, 0.98,
            f"$N_{{QQ}}={qq['nsd'].size:,}$\n$N_{{GG}}={gg['nsd'].size:,}$",
            transform=ax.transAxes, ha="left", va="top", fontsize=11,
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "0.7"})

    fig.tight_layout()
    out_pdf = Path("roc_psi_vs_others.pdf")
    fig.savefig(out_pdf); plt.close(fig)
    log(f"\nWrote {out_pdf.resolve()}")
    Path("roc_psi_vs_others.log").write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
