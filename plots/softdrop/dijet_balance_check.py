#!/usr/bin/env python3
"""Dijet balance sanity check for QQ / GG / GQ samples.

At LO the two outgoing hard partons are back-to-back in phi and
balanced in pT. Initial- and final-state radiation smears these,
but the bulk of the distribution should stay Born-like. A broad
Delta-phi peak or a strongly asymmetric pT-ratio distribution is a
warning that the "LO-parton -> jet" identification (on which our
QQ/GG/GQ labelling rests) is being diluted by hard extra emissions.

For each jetreco_softdrop output ROOT, read the already-computed
dijet branches (events with >=2 jets only) and produce:

    dijet_balance.pdf            two-panel figure: Delta-phi and pT ratio
    dijet_balance_summary.log    numeric summary with Born-cut fractions

Outputs into the current working directory. Intended to be run inside
data-jets/<sample>_pT7/ next to alljets_*.root.

Usage:
    python dijet_balance_check.py <root_file>
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

# Born-like selection: back-to-back and reasonably balanced.
DPHI_BORN = 2.5      # rad; standard HERA dijet azimuthal cut
PTRATIO_BORN = 0.5   # ET_subleading / ET_leading


def load(root_path):
    """Return dict: cat -> {dphi, pt_ratio} as flat 1-D numpy arrays.
    Only events with n_jets >= 2 contribute (the dijet branches are
    filled only for those)."""
    data = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            tree = f[key]
            dphi = tree["dijet_delta_phi"].array(library="np")
            et_lead = tree["leading_jet_et"].array(library="np")
            et_sub = tree["subleading_jet_et"].array(library="np")
            # The dijet branches are filled with -999 for events with
            # <2 jets; filter those out.
            mask = (et_lead > 0) & (et_sub > 0) & (dphi >= 0)
            if mask.sum() == 0:
                continue
            data[cat] = {
                "dphi":     dphi[mask],
                "pt_ratio": et_sub[mask] / et_lead[mask],
                "et_lead":  et_lead[mask],
                "et_sub":   et_sub[mask],
            }
    return data


def fmt_pct(n, tot):
    return f"{100*n/tot:5.2f}%" if tot else "  —   "


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: plots/softdrop/output/<sample>/)")
    args = ap.parse_args()

    data = load(args.root_file)
    if not data:
        print("No usable dijet data in", args.root_file)
        sys.exit(1)

    sample = Path(args.root_file).resolve().parent.name
    out_dir = Path(args.out_dir) if args.out_dir else Path(__file__).resolve().parent / "output" / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_name = Path(args.root_file).stem
    log_lines = []
    def log(s=""):
        print(s); log_lines.append(s)

    log("=" * 88)
    log(f"Dijet balance check  —  {sample_name}")
    log("=" * 88)
    log(f"Input: {args.root_file}")
    log(f"Born-like cuts used for fraction: Delta-phi > {DPHI_BORN} rad "
        f"AND pT_ratio > {PTRATIO_BORN}")
    log("")

    log(f"{'Category':<14}{'N dijets':>12}"
        f"{'<Delta-phi>':>14}{'<pT_ratio>':>14}"
        f"{'Born-like':>14}{'Delta-phi>2.5':>16}{'pTr>0.5':>12}")
    log("-" * 96)

    for cat in CATEGORIES:
        if cat not in data:
            log(f"{cat:<14}  (not found)")
            continue
        d = data[cat]
        n = d["dphi"].size
        dphi_mean = d["dphi"].mean()
        ptr_mean = d["pt_ratio"].mean()
        born = ((d["dphi"] > DPHI_BORN) & (d["pt_ratio"] > PTRATIO_BORN)).sum()
        n_dphi = (d["dphi"] > DPHI_BORN).sum()
        n_ptr = (d["pt_ratio"] > PTRATIO_BORN).sum()
        log(f"{cat:<14}{n:>12}"
            f"{dphi_mean:>14.3f}{ptr_mean:>14.3f}"
            f"{fmt_pct(born, n):>14}"
            f"{fmt_pct(n_dphi, n):>16}"
            f"{fmt_pct(n_ptr, n):>12}")

    log("")
    log("Interpretation:")
    log("  - Delta-phi should peak sharply near pi (LO back-to-back).")
    log("  - pT_ratio should peak near 1 (LO momentum balance).")
    log("  - Born-like fraction >~70% is expected with pT_hat > 7 GeV.")
    log("  - A strongly skewed pT-ratio or broad Delta-phi can dilute the")
    log("    'hard-parton -> jet' identification used for QQ/GG/GQ labels.")

    # --- plots ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    # Delta-phi
    ax = axes[0]
    for cat in CATEGORIES:
        if cat not in data:
            continue
        ax.hist(data[cat]["dphi"], bins=50, range=(0, np.pi),
                density=True, histtype="step", linewidth=2,
                color=COLORS[cat], label=LABELS[cat])
    ax.axvline(DPHI_BORN, linestyle="--", color="grey", alpha=0.5,
               label=f"$\\Delta\\phi$ = {DPHI_BORN} cut")
    ax.axvline(np.pi, linestyle=":", color="black", alpha=0.4,
               label=r"$\pi$ (LO)")
    ax.set_xlabel(r"$\Delta\phi$(jet$_1$, jet$_2$)  [rad]")
    ax.set_ylabel("fraction of dijet events")
    ax.set_title(f"{sample_name}")
    ax.set_xlim(0, np.pi)
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.3)

    # pT ratio
    ax = axes[1]
    for cat in CATEGORIES:
        if cat not in data:
            continue
        ax.hist(data[cat]["pt_ratio"], bins=50, range=(0, 1),
                density=True, histtype="step", linewidth=2,
                color=COLORS[cat], label=LABELS[cat])
    ax.axvline(PTRATIO_BORN, linestyle="--", color="grey", alpha=0.5,
               label=f"$p_T$-ratio > {PTRATIO_BORN} cut")
    ax.axvline(1.0, linestyle=":", color="black", alpha=0.4,
               label="1 (LO)")
    ax.set_xlabel(r"$E_T$(jet$_2$) / $E_T$(jet$_1$)")
    ax.set_ylabel("fraction of dijet events")
    ax.set_title(f"{sample_name}")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_pdf = out_dir / "dijet_balance.pdf"
    fig.savefig(out_pdf)
    plt.close(fig)

    log(f"\nWrote {out_pdf.resolve()}")

    out_log = out_dir / "dijet_balance_summary.log"
    out_log.write_text("\n".join(log_lines) + "\n")
    print(f"Wrote {out_log.resolve()}")


if __name__ == "__main__":
    main()
