#!/usr/bin/env python3
"""Compare softdrop multiplicity (n_sd) vs kT exclusive n_subjets as
quark/gluon discriminants on the same jets.

Addresses reviewer point N4 on Photoproduction_v1.pdf:
"it would be easy to check other observables well-known to be good
discriminants like softdrop multiplicity, which are specifically known
to be better than n_subjets."

Usage:
    python compare_nsd_vs_nsubjets.py <root_file> [--eta-binned]

Writes:
    compare_nsd_vs_nsubjets.pdf          main plot
    compare_nsd_vs_nsubjets_summary.txt  numeric summary for pasting
"""

import argparse
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt


CATEGORIES = ["QQ_Events", "GG_Events", "GQ_Events"]
LABELS = {"QQ_Events": "QQ", "GG_Events": "GG", "GQ_Events": "GQ"}
COLORS = {"QQ_Events": "C0", "GG_Events": "C3", "GQ_Events": "C2"}

ETA_BINS = [(-1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]


# -----------------------------------------------------------------------------
def load(path):
    """Return dict: category -> dict{'nsd','nsubjets','eta'} of flat arrays."""
    data = {}
    with uproot.open(path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            tree = f[key]
            nsd = tree["jet_nsd"].array(library="np")
            nsj = tree["jet_nsubjets"].array(library="np")
            eta = tree["jet_eta"].array(library="np")
            if len(nsd) == 0:
                continue
            nsd_flat = np.concatenate([np.asarray(a) for a in nsd if len(a) > 0])
            nsj_flat = np.concatenate([np.asarray(a) for a in nsj if len(a) > 0])
            eta_flat = np.concatenate([np.asarray(a) for a in eta if len(a) > 0])
            data[cat] = {"nsd": nsd_flat, "nsubjets": nsj_flat, "eta": eta_flat}
    return data


# -----------------------------------------------------------------------------
def separation(a, b):
    """Mean separation in units of pooled sigma: (<b>-<a>) / sqrt(s_a^2 + s_b^2)."""
    da = a.mean() - b.mean()
    return (b.mean() - a.mean()) / np.sqrt(a.var(ddof=1) + b.var(ddof=1))


def roc(sig, bkg, values):
    """ROC for integer-valued discriminants.

    Signal = higher-value class (we treat GG as signal for Q/G tagging).
    For each possible cut c, efficiency is P(X >= c).
    Returns (eff_sig, eff_bkg, auc).
    """
    grid = np.arange(values.min(), values.max() + 2)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    # Both start at 1 (cut = min) and end at 0 (cut > max). Integrate:
    # AUC = integral of eff_sig d(eff_bkg) from 1 to 0, flipped to positive.
    # Use trapezoid on sorted eff_bkg ascending.
    order = np.argsort(eff_bkg)
    auc = np.trapz(eff_sig[order], eff_bkg[order])
    return eff_sig, eff_bkg, auc


# -----------------------------------------------------------------------------
def summarize(data, summary_lines):
    def line(s):
        print(s)
        summary_lines.append(s)

    line("=" * 74)
    line("HEAD-TO-HEAD: n_sd  vs  n_subjets")
    line("=" * 74)
    line(f"{'Category':<12}{'N_jets':>10}"
         f"{'<n_sd>':>10}{'sd_std':>10}"
         f"{'<n_subj>':>10}{'sj_std':>10}")
    for cat in CATEGORIES:
        if cat not in data:
            continue
        d = data[cat]
        line(f"{LABELS[cat]:<12}{d['nsd'].size:>10d}"
             f"{d['nsd'].mean():>10.2f}{d['nsd'].std(ddof=1):>10.2f}"
             f"{d['nsubjets'].mean():>10.2f}{d['nsubjets'].std(ddof=1):>10.2f}")

    if "QQ_Events" in data and "GG_Events" in data:
        qq, gg = data["QQ_Events"], data["GG_Events"]

        sep_nsd = separation(qq["nsd"], gg["nsd"])
        sep_nsj = separation(qq["nsubjets"], gg["nsubjets"])
        line("")
        line(f"Mean separation  (<GG> - <QQ>) / sqrt(sig_QQ^2 + sig_GG^2):")
        line(f"    n_sd       = {sep_nsd:+.3f}")
        line(f"    n_subjets  = {sep_nsj:+.3f}")
        if sep_nsd > sep_nsj:
            line(f"    -> n_sd separation is {sep_nsd / sep_nsj:.2f}x that of n_subjets")
        else:
            line(f"    -> n_subjets separation is {sep_nsj / sep_nsd:.2f}x that of n_sd")

        all_vals = np.concatenate([qq["nsd"], gg["nsd"], qq["nsubjets"], gg["nsubjets"]])
        _, _, auc_nsd = roc(gg["nsd"], qq["nsd"], all_vals)
        _, _, auc_nsj = roc(gg["nsubjets"], qq["nsubjets"], all_vals)
        line("")
        line("AUC (area under ROC, GG as signal vs QQ as background):")
        line(f"    n_sd       = {auc_nsd:.4f}")
        line(f"    n_subjets  = {auc_nsj:.4f}")
        line(f"    (random = 0.5, perfect = 1.0; higher is better)")


# -----------------------------------------------------------------------------
def plot_main(data, out_pdf):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # (1) n_sd distributions
    ax = axes[0]
    for cat in ["QQ_Events", "GG_Events", "GQ_Events"]:
        if cat not in data:
            continue
        ax.hist(data[cat]["nsd"],
                bins=np.arange(-0.5, 12.5, 1),
                density=True, histtype="step", linewidth=2,
                label=LABELS[cat], color=COLORS[cat])
    ax.set_xlabel(r"$n_{\mathrm{SD}}$")
    ax.set_ylabel("Fraction of jets")
    ax.set_title(r"Softdrop multiplicity ($z_{\mathrm{cut}}=0.1,\ \beta=0$)")
    ax.legend(frameon=False)
    ax.set_xlim(-0.5, 12)

    # (2) n_subjets distributions
    ax = axes[1]
    for cat in ["QQ_Events", "GG_Events", "GQ_Events"]:
        if cat not in data:
            continue
        ax.hist(data[cat]["nsubjets"],
                bins=np.arange(-0.5, 12.5, 1),
                density=True, histtype="step", linewidth=2,
                label=LABELS[cat], color=COLORS[cat])
    ax.set_xlabel(r"$n_{\mathrm{subjets}}$")
    ax.set_ylabel("Fraction of jets")
    ax.set_title(r"$k_T$ exclusive subjets ($y_{\mathrm{cut}}=5\times10^{-4}$)")
    ax.legend(frameon=False)
    ax.set_xlim(-0.5, 12)

    # (3) ROC curves
    ax = axes[2]
    if "QQ_Events" in data and "GG_Events" in data:
        qq, gg = data["QQ_Events"], data["GG_Events"]
        for obs, label, col in [("nsd", r"$n_{\mathrm{SD}}$", "C1"),
                                ("nsubjets", r"$n_{\mathrm{subjets}}$", "C4")]:
            all_vals = np.concatenate([qq[obs], gg[obs]])
            eff_sig, eff_bkg, auc = roc(gg[obs], qq[obs], all_vals)
            ax.plot(eff_bkg, eff_sig, "o-", label=f"{label} (AUC={auc:.3f})",
                    color=col, markersize=5)
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="random")
        ax.set_xlabel(r"$\varepsilon$(QQ)")
        ax.set_ylabel(r"$\varepsilon$(GG)")
        ax.set_title("ROC: GG as signal, QQ as background")
        ax.legend(frameon=False, loc="lower right")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)

    fig.tight_layout()
    fig.savefig(out_pdf)
    print(f"\nWrote {out_pdf}")


# -----------------------------------------------------------------------------
def eta_binned_summary(data, summary_lines):
    """Per-eta-bin separation for both observables."""
    if "QQ_Events" not in data or "GG_Events" not in data:
        return

    def line(s):
        print(s)
        summary_lines.append(s)

    line("")
    line("=" * 74)
    line("Per-eta-bin separation (<GG> - <QQ>) / pooled_sigma")
    line("=" * 74)
    line(f"{'eta bin':<14}{'n_QQ':>8}{'n_GG':>8}"
         f"{'sep(n_sd)':>12}{'sep(n_subj)':>14}")
    for lo, hi in ETA_BINS:
        mqq = (data["QQ_Events"]["eta"] >= lo) & (data["QQ_Events"]["eta"] < hi)
        mgg = (data["GG_Events"]["eta"] >= lo) & (data["GG_Events"]["eta"] < hi)
        if mqq.sum() < 20 or mgg.sum() < 20:
            line(f"{lo:>5.1f},{hi:>4.1f}   (too few jets — skipped)")
            continue
        s_nsd = separation(data["QQ_Events"]["nsd"][mqq],
                           data["GG_Events"]["nsd"][mgg])
        s_nsj = separation(data["QQ_Events"]["nsubjets"][mqq],
                           data["GG_Events"]["nsubjets"][mgg])
        line(f"[{lo:>4.1f},{hi:>4.1f})"
             f"{mqq.sum():>8d}{mgg.sum():>8d}"
             f"{s_nsd:>12.3f}{s_nsj:>14.3f}")


# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root_file")
    ap.add_argument("--eta-binned", action="store_true",
                    help="Also print per-eta-bin comparison table")
    args = ap.parse_args()

    data = load(args.root_file)
    if not data:
        print(f"No usable trees in {args.root_file}")
        sys.exit(1)

    summary_lines = []
    summarize(data, summary_lines)
    if args.eta_binned:
        eta_binned_summary(data, summary_lines)

    with open("compare_nsd_vs_nsubjets_summary.txt", "w") as fh:
        fh.write("\n".join(summary_lines) + "\n")
    print("\nWrote compare_nsd_vs_nsubjets_summary.txt")

    plot_main(data, "compare_nsd_vs_nsubjets.pdf")


if __name__ == "__main__":
    main()
