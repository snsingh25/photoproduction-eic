#!/usr/bin/env python3
"""Sanity check for the per-eta-bin ROC results.

Expected behaviour for a well-behaved q/g discriminant X:
  Inclusive: <X>_QQ and <X>_GG differ substantially -> high AUC.
  Per-bin:   they differ less (binning removes kinematic variation).

The current ROC scan shows the OPPOSITE for n_subjets: per-bin
separation ~0.7-0.84, inclusive only 0.21. This script explicitly
checks whether the eta slicing is doing the right thing by printing:

  1. Per-bin <X>, sigma, N for QQ and GG (both observables).
  2. The variance decomposition
         sigma^2_inclusive = < sigma^2_bin >_bins + Var_bins(<X>_bin)
     which should match the inclusive sigma to within a few %.
  3. An independent check that the eta slicing picks the same jets
     the inclusive pool contains (sum of per-bin N = inclusive N).
  4. Plot: <X>_QQ(eta) and <X>_GG(eta) as step functions.

Usage:
    python sanity_eta_means.py <alljets_*.root>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt


CATEGORIES = ["QQ_Events", "GG_Events"]
BINS = [(-1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]


def load(root_path):
    d = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            tree = f[key]
            nsd = np.concatenate(
                [np.asarray(a) for a in tree["jet_nsd"].array(library="np") if len(a) > 0]
            )
            nsj = np.concatenate(
                [np.asarray(a) for a in tree["jet_nsubjets"].array(library="np") if len(a) > 0]
            )
            eta = np.concatenate(
                [np.asarray(a) for a in tree["jet_eta"].array(library="np") if len(a) > 0]
            )
            d[cat] = {"nsd": nsd, "nsubjets": nsj, "eta": eta}
    return d


def sep(a, b):
    return (b.mean() - a.mean()) / np.sqrt(a.var(ddof=1) + b.var(ddof=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root_file")
    args = ap.parse_args()

    data = load(args.root_file)
    if "QQ_Events" not in data or "GG_Events" not in data:
        print("missing QQ or GG tree"); sys.exit(1)

    sample_name = Path(args.root_file).stem

    print("=" * 92)
    print(f"Sanity check for per-eta binning  —  {sample_name}")
    print("=" * 92)

    for obs in ("nsubjets", "nsd"):
        print(f"\n--- observable: {obs} ---")

        # inclusive
        qq_all = data["QQ_Events"][obs]
        gg_all = data["GG_Events"][obs]
        qq_eta_all = data["QQ_Events"]["eta"]
        gg_eta_all = data["GG_Events"]["eta"]

        # Restrict to |eta| < 4 so per-bin and inclusive use the same pool.
        mask_qq_all = (qq_eta_all >= -4) & (qq_eta_all < 4)
        mask_gg_all = (gg_eta_all >= -4) & (gg_eta_all < 4)
        qq_all = qq_all[mask_qq_all]
        gg_all = gg_all[mask_gg_all]
        qq_eta_all = qq_eta_all[mask_qq_all]
        gg_eta_all = gg_eta_all[mask_gg_all]

        mu_qq_incl = qq_all.mean();  var_qq_incl = qq_all.var(ddof=1)
        mu_gg_incl = gg_all.mean();  var_gg_incl = gg_all.var(ddof=1)

        print(f"INCLUSIVE (|eta|<4):")
        print(f"  QQ: N={qq_all.size:6d}  <X>={mu_qq_incl:.4f}  sigma={np.sqrt(var_qq_incl):.4f}")
        print(f"  GG: N={gg_all.size:6d}  <X>={mu_gg_incl:.4f}  sigma={np.sqrt(var_gg_incl):.4f}")
        print(f"  Delta_mu = {mu_gg_incl-mu_qq_incl:+.4f}")
        print(f"  sep      = {sep(qq_all, gg_all):+.4f}")

        print(f"\n{'bin':<14}{'N_QQ':>8}{'<QQ>':>10}{'sig_QQ':>10}"
              f"{'N_GG':>8}{'<GG>':>10}{'sig_GG':>10}{'Delta':>10}{'sep':>10}")

        # keep per-bin stats for the decomposition check
        per_bin_qq = []
        per_bin_gg = []
        for lo, hi in BINS:
            mqq = (qq_eta_all >= lo) & (qq_eta_all < hi)
            mgg = (gg_eta_all >= lo) & (gg_eta_all < hi)
            qb = qq_all[mqq]
            gb = gg_all[mgg]
            if qb.size < 2 or gb.size < 2:
                print(f"[{lo:+4.1f},{hi:+4.1f})  —  too few jets")
                continue
            per_bin_qq.append((qb.size, qb.mean(), qb.var(ddof=1)))
            per_bin_gg.append((gb.size, gb.mean(), gb.var(ddof=1)))
            print(f"[{lo:+4.1f},{hi:+4.1f})   "
                  f"{qb.size:>5}{qb.mean():>10.4f}{np.sqrt(qb.var(ddof=1)):>10.4f}"
                  f"{gb.size:>8}{gb.mean():>10.4f}{np.sqrt(gb.var(ddof=1)):>10.4f}"
                  f"{gb.mean()-qb.mean():>+10.4f}{sep(qb,gb):>+10.4f}")

        # Variance decomposition: sigma^2_total = E[sigma^2_bin] + Var(mu_bin)
        def decompose(bin_stats, all_arr):
            Ns = np.array([b[0] for b in bin_stats], dtype=float)
            mus = np.array([b[1] for b in bin_stats])
            vars_ = np.array([b[2] for b in bin_stats])
            tot_N = Ns.sum()
            mu_pool = (Ns * mus).sum() / tot_N
            within = (Ns * vars_).sum() / tot_N
            between = (Ns * (mus - mu_pool)**2).sum() / tot_N
            reconstructed = within + between
            observed = all_arr.var(ddof=1) * (all_arr.size - 1) / all_arr.size  # population variance
            return mu_pool, within, between, reconstructed, observed

        mu_p, wQ, bQ, rQ, obsQ = decompose(per_bin_qq, qq_all)
        mu_pG, wG, bG, rG, obsG = decompose(per_bin_gg, gg_all)

        print(f"\nvariance decomposition (sigma^2 = within + between):")
        print(f"  QQ:  within^2 = {wQ:.4f}   between^2 = {bQ:.4f}"
              f"   sum = {rQ:.4f}   observed^2 = {obsQ:.4f}   ratio = {rQ/obsQ:.3f}")
        print(f"  GG:  within^2 = {wG:.4f}   between^2 = {bG:.4f}"
              f"   sum = {rG:.4f}   observed^2 = {obsG:.4f}   ratio = {rG/obsG:.3f}")
        print(f"  -> if ratio ~= 1.0, per-bin + between-bin variance reconstructs inclusive.")
        print(f"  -> fraction of QQ variance that is between-bin: {bQ/rQ*100:.1f}%")
        print(f"  -> fraction of GG variance that is between-bin: {bG/rG*100:.1f}%")

        # Sanity: N should sum correctly
        print(f"\nN sanity: sum per-bin N_QQ = {sum(b[0] for b in per_bin_qq)}  vs inclusive N_QQ = {qq_all.size}  (diff = out-of-bin jets)")
        print(f"          sum per-bin N_GG = {sum(b[0] for b in per_bin_gg)}  vs inclusive N_GG = {gg_all.size}")

    # --- plot <X>(eta) for QQ and GG ---------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    eta_fine = np.linspace(-4, 4, 33)
    eta_centers = 0.5 * (eta_fine[:-1] + eta_fine[1:])

    for ax, obs, ylabel in zip(axes,
                               ("nsubjets", "nsd"),
                               (r"$\langle n_{\mathrm{subjets}}\rangle$",
                                r"$\langle n_{\mathrm{SD}}\rangle$")):
        for cat, label, col in [("QQ_Events", "QQ", "C0"),
                                ("GG_Events", "GG", "C3")]:
            eta = data[cat]["eta"]
            x = data[cat][obs]
            means = []
            counts = []
            for lo, hi in zip(eta_fine[:-1], eta_fine[1:]):
                m = (eta >= lo) & (eta < hi)
                if m.sum() < 5:
                    means.append(np.nan); counts.append(0)
                else:
                    means.append(x[m].mean()); counts.append(m.sum())
            means = np.array(means); counts = np.array(counts)
            ax.plot(eta_centers, means, "o-", color=col, label=f"{label}")
        ax.set_xlabel(r"$\eta$")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{sample_name}")
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)

    fig.tight_layout()
    out_pdf = Path("sanity_eta_means.pdf")
    fig.savefig(out_pdf)
    print(f"\nWrote {out_pdf.resolve()}")


if __name__ == "__main__":
    main()
