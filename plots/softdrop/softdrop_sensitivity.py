#!/usr/bin/env python3
"""Softdrop-parameter sensitivity of the q/g tagging performance.

For each (z_cut, beta) configuration precomputed by jetreco_softdrop.cc
(jet_nsd, jet_nsd_beta1, jet_nsd_loose), reads HERA 300 output and
tabulates <n_SD>_QQ, <n_SD>_GG, AUC (inclusive), AUC in eta=[0,1)
and AUC in eta=[1.5, 2).

Usage:
    python softdrop_sensitivity.py <root_file>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot


CATEGORIES = ["QQ_Events", "GG_Events"]

# (branch_name, display_label, (z_cut, beta), short_tag)
CONFIGS = [
    ("jet_nsd",        "(0.1, 0.0)  baseline",    (0.1,  0.0)),
    ("jet_nsd_beta1",  "(0.1, 1.0)  original SD", (0.1,  1.0)),
    ("jet_nsd_loose",  "(0.05, 0.0) loose",       (0.05, 0.0)),
]


def load(root_path):
    data = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            nsd   = t["jet_nsd"].array(library="np")
            nsdb  = t["jet_nsd_beta1"].array(library="np")
            nsdl  = t["jet_nsd_loose"].array(library="np")
            eta   = t["jet_eta"].array(library="np")
            data[cat] = {
                "jet_nsd":       np.concatenate([np.asarray(a) for a in nsd  if len(a) > 0]),
                "jet_nsd_beta1": np.concatenate([np.asarray(a) for a in nsdb if len(a) > 0]),
                "jet_nsd_loose": np.concatenate([np.asarray(a) for a in nsdl if len(a) > 0]),
                "eta":           np.concatenate([np.asarray(a) for a in eta  if len(a) > 0]),
            }
    return data


def roc_auc(sig, bkg):
    vmin = int(min(sig.min(), bkg.min()))
    vmax = int(max(sig.max(), bkg.max()))
    grid = np.arange(vmin, vmax + 2)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    return np.trapezoid(eff_sig[order], eff_bkg[order])


def slice_eta(data, cat, lo, hi, branch):
    m = (data[cat]["eta"] >= lo) & (data[cat]["eta"] < hi)
    return data[cat][branch][m]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: plots/softdrop/output/<sample>/)")
    args = ap.parse_args()

    data = load(args.root_file)
    if "QQ_Events" not in data or "GG_Events" not in data:
        print("missing QQ or GG tree"); sys.exit(1)

    sample = Path(args.root_file).resolve().parent.name
    out_dir = Path(args.out_dir) if args.out_dir else Path(__file__).resolve().parent / "output" / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_name = Path(args.root_file).stem
    log_lines = []
    def log(s=""):
        print(s); log_lines.append(s)

    log("=" * 102)
    log(f"Softdrop-parameter sensitivity  —  {sample_name}")
    log("=" * 102)
    log(f"{'Config':<24}{'<nSD>_QQ':>10}{'sig_QQ':>10}{'<nSD>_GG':>10}{'sig_GG':>10}"
        f"{'AUC incl':>12}{'AUC [0,1)':>12}{'AUC [1.5,2)':>14}")
    log("-" * 102)

    for branch, label, (_zc, _b) in CONFIGS:
        qq_all = data["QQ_Events"][branch]
        gg_all = data["GG_Events"][branch]
        mu_qq = qq_all.mean(); sig_qq = qq_all.std(ddof=1)
        mu_gg = gg_all.mean(); sig_gg = gg_all.std(ddof=1)
        auc_incl = roc_auc(gg_all, qq_all)

        qq_01 = slice_eta(data, "QQ_Events", 0.0, 1.0, branch)
        gg_01 = slice_eta(data, "GG_Events", 0.0, 1.0, branch)
        auc_01 = roc_auc(gg_01, qq_01) if qq_01.size >= 20 and gg_01.size >= 20 else np.nan

        qq_15 = slice_eta(data, "QQ_Events", 1.5, 2.0, branch)
        gg_15 = slice_eta(data, "GG_Events", 1.5, 2.0, branch)
        auc_15 = roc_auc(gg_15, qq_15) if qq_15.size >= 20 and gg_15.size >= 20 else np.nan

        def fmt(x): return f"{x:>12.4f}" if not np.isnan(x) else f"{'—':>12}"
        def fmt2(x): return f"{x:>14.4f}" if not np.isnan(x) else f"{'—':>14}"

        log(f"{label:<24}{mu_qq:>10.3f}{sig_qq:>10.3f}{mu_gg:>10.3f}{sig_gg:>10.3f}"
            f"{auc_incl:>12.4f}{fmt(auc_01)}{fmt2(auc_15)}")

    log("")
    log("Physics interpretation:")
    log("  - Stricter cut (lower z_cut) counts more splittings, raises <n_SD>")
    log("    for both classes, but also lowers per-splitting gluon-vs-quark")
    log("    Casimir leverage -> AUC may not improve much.")
    log("  - Non-zero beta (angle-dependent grooming) tightens the cut at")
    log("    small angles and loosens it at large ones -> counts fewer")
    log("    splittings and can approach the Cambridge/Aachen ungroomed")
    log("    limit as beta -> large.")

    out_log = out_dir / "softdrop_sensitivity.log"
    out_log.write_text("\n".join(log_lines) + "\n")
    print(f"\nWrote {out_log.resolve()}")


if __name__ == "__main__":
    main()
