#!/usr/bin/env python3
"""ROC curves for n_SD and n_subjets: per-eta-bin and inclusive.

Reads the output of `jetreco_softdrop` (jet_nsd / jet_nsubjets / jet_eta
branches in per-category TTrees) and produces:

    roc_inclusive.pdf             overall ROC (eta-pooled)
    roc_eta_<lo>_<hi>.pdf         one PDF per eta bin
    roc_all_bins.pdf              multi-panel grid (inclusive + all bins)
    roc_summary.log               numerical details and per-bin counts

All outputs go into the current working directory. Intended to be run
inside `data-jets/<sample>_pT7/` so the outputs sit next to the
`alljets_*.root` file.

Usage:
    python roc_eta_scan.py <root_file> [--bins <lo>,<hi>,...]
    python roc_eta_scan.py alljets_hera300_pT7_R10_EtMin17.root

Options:
    --bins  comma-separated eta bin edges, e.g. "-2,-1,0,1,2,3,4"
            (default covers full |eta|<4 acceptance in 1-unit steps)
    --min-jets  minimum QQ or GG jet count to include a bin (default 20)
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt


CATEGORIES = ["QQ_Events", "GG_Events"]
DEFAULT_BIN_EDGES = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]

COLOR_NSD = "C1"        # orange
COLOR_NSUBJETS = "C4"   # purple


# -----------------------------------------------------------------------------
def load(root_path):
    """Return dict: cat -> {nsd, nsubjets, eta} flat numpy arrays."""
    data = {}
    with uproot.open(root_path) as f:
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
            data[cat] = {
                "nsd":      np.concatenate([np.asarray(a) for a in nsd if len(a) > 0]),
                "nsubjets": np.concatenate([np.asarray(a) for a in nsj if len(a) > 0]),
                "eta":      np.concatenate([np.asarray(a) for a in eta if len(a) > 0]),
            }
    return data


def slice_by_eta(data, lo, hi):
    """Return subset of data with eta in [lo, hi). None if empty."""
    out = {}
    for cat, d in data.items():
        m = (d["eta"] >= lo) & (d["eta"] < hi)
        if m.sum() == 0:
            continue
        out[cat] = {k: v[m] for k, v in d.items()}
    return out


# -----------------------------------------------------------------------------
def roc_integer(sig, bkg):
    """ROC for integer-valued discriminants.

    Threshold c ranges over all reachable integer cut values; efficiency
    is the fraction of jets at or above c. GG (sig) should have larger
    typical values than QQ (bkg), so higher cuts cull QQ first.

    Returns (eff_sig, eff_bkg, auc) with arrays sorted by ascending
    eff_bkg, ready to plot on (QQ efficiency, GG efficiency) axes.
    """
    vmin = int(min(sig.min(), bkg.min()))
    vmax = int(max(sig.max(), bkg.max()))
    grid = np.arange(vmin, vmax + 2)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    eff_sig_s = eff_sig[order]
    eff_bkg_s = eff_bkg[order]
    auc = np.trapezoid(eff_sig_s, eff_bkg_s)
    return eff_sig_s, eff_bkg_s, auc


def separation(qq, gg):
    """(<gg> - <qq>) / sqrt(var_qq + var_gg). Positive if gg is higher."""
    return (gg.mean() - qq.mean()) / np.sqrt(
        qq.var(ddof=1) + gg.var(ddof=1)
    )


# -----------------------------------------------------------------------------
def plot_roc_single(sub, title, out_pdf, stats):
    """Draw a single ROC panel for n_SD vs n_subjets into `out_pdf`."""
    fig, ax = plt.subplots(figsize=(6, 5.5))

    qq = sub["QQ_Events"]
    gg = sub["GG_Events"]

    for obs, label, col in [
        ("nsd",      r"$n_{\mathrm{SD}}$",     COLOR_NSD),
        ("nsubjets", r"$n_{\mathrm{subjets}}$", COLOR_NSUBJETS),
    ]:
        es, eb, auc = roc_integer(gg[obs], qq[obs])
        ax.plot(eb, es, "o-", color=col, linewidth=2, markersize=6,
                label=f"{label}  AUC = {auc:.3f}")
        stats[obs] = {"auc": auc, "sep": separation(qq[obs], gg[obs])}

    ax.plot([0, 1], [0, 1], "--", color="grey", alpha=0.5, label="random")
    ax.set_xlabel(r"$\varepsilon$(QQ)  — background acceptance")
    ax.set_ylabel(r"$\varepsilon$(GG)  — signal efficiency")
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.3)

    nqq = qq["nsd"].size
    ngg = gg["nsd"].size
    ax.text(0.02, 0.98, f"$N_{{QQ}} = {nqq:,}$\n$N_{{GG}} = {ngg:,}$",
            transform=ax.transAxes, ha="left", va="top", fontsize=11,
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "0.7"})

    stats["n_qq"] = nqq
    stats["n_gg"] = ngg

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_grid(panels, out_pdf, sample_name):
    """Multi-panel grid. `panels` is list of (title, sub) with the
    inclusive panel first, then per-eta-bin panels."""
    n = len(panels)
    if n == 0:
        return
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.3 * nrows))
    axes = np.atleast_2d(axes).ravel()

    for ax, (title, sub) in zip(axes, panels):
        qq = sub["QQ_Events"]
        gg = sub["GG_Events"]
        for obs, label, col in [
            ("nsd",      r"$n_{\mathrm{SD}}$",     COLOR_NSD),
            ("nsubjets", r"$n_{\mathrm{subjets}}$", COLOR_NSUBJETS),
        ]:
            es, eb, auc = roc_integer(gg[obs], qq[obs])
            ax.plot(eb, es, "o-", color=col, linewidth=1.8, markersize=4,
                    label=f"{label}  AUC={auc:.3f}")
        ax.plot([0, 1], [0, 1], "--", color="grey", alpha=0.4)
        ax.set_xlabel(r"$\varepsilon$(QQ)")
        ax.set_ylabel(r"$\varepsilon$(GG)")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
        ax.set_title(title, fontsize=11)
        ax.legend(frameon=False, loc="lower right", fontsize=9)
        ax.grid(alpha=0.3)
        ax.text(0.02, 0.98,
                f"$N_{{QQ}}={qq['nsd'].size:,}$\n$N_{{GG}}={gg['nsd'].size:,}$",
                transform=ax.transAxes, ha="left", va="top", fontsize=9,
                bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "0.7"})

    # Hide unused axes
    for ax in axes[len(panels):]:
        ax.set_visible(False)

    fig.suptitle(f"{sample_name} — ROC curves for $n_{{SD}}$ vs $n_{{subjets}}$",
                 fontsize=13, y=1.00)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file", help="jetreco_softdrop output ROOT")
    ap.add_argument("--bins", default=",".join(str(x) for x in DEFAULT_BIN_EDGES),
                    help="comma-separated eta bin edges")
    ap.add_argument("--min-jets", type=int, default=20,
                    help="skip bins with fewer QQ or GG jets than this")
    args = ap.parse_args()

    bin_edges = [float(x) for x in args.bins.split(",")]
    bins = list(zip(bin_edges[:-1], bin_edges[1:]))

    sample_name = Path(args.root_file).stem
    log_path = Path("roc_summary.log")
    log_lines = []

    def log(s=""):
        print(s)
        log_lines.append(s)

    log("=" * 78)
    log(f"ROC scan  —  sample: {sample_name}")
    log(f"input:   {args.root_file}")
    log(f"bins:    {bins}  (min {args.min_jets} jets per class to include)")
    log("=" * 78)

    data = load(args.root_file)
    if not all(c in data for c in CATEGORIES):
        log("ERROR: missing QQ_Events or GG_Events tree — nothing to do.")
        log_path.write_text("\n".join(log_lines) + "\n")
        sys.exit(1)

    n_qq_total = data["QQ_Events"]["nsd"].size
    n_gg_total = data["GG_Events"]["nsd"].size
    log(f"\nLoaded: N_QQ = {n_qq_total:,},  N_GG = {n_gg_total:,}")

    panels = []   # list of (title, sub) for the grid
    rows = []     # list of dicts for summary table

    # --- inclusive ---------------------------------------------------------
    log("\n" + "-" * 78)
    log("Inclusive (eta-pooled)")
    log("-" * 78)
    inc_stats = {}
    plot_roc_single(data, f"{sample_name}  —  inclusive",
                    "roc_inclusive.pdf", inc_stats)
    log(f"  N_QQ={inc_stats['n_qq']:,}  N_GG={inc_stats['n_gg']:,}")
    log(f"  n_SD       : AUC = {inc_stats['nsd']['auc']:.4f}   sep = {inc_stats['nsd']['sep']:+.3f}")
    log(f"  n_subjets  : AUC = {inc_stats['nsubjets']['auc']:.4f}   sep = {inc_stats['nsubjets']['sep']:+.3f}")
    log(f"  -> wrote roc_inclusive.pdf")
    panels.append((f"inclusive  $(|\\eta|<4)$", data))
    rows.append({
        "bin": "inclusive", "n_qq": inc_stats["n_qq"], "n_gg": inc_stats["n_gg"],
        "auc_nsd": inc_stats["nsd"]["auc"], "sep_nsd": inc_stats["nsd"]["sep"],
        "auc_nsj": inc_stats["nsubjets"]["auc"], "sep_nsj": inc_stats["nsubjets"]["sep"],
    })

    # --- per-eta-bin -------------------------------------------------------
    for lo, hi in bins:
        log("\n" + "-" * 78)
        log(f"eta in [{lo:+.1f}, {hi:+.1f})")
        log("-" * 78)
        sub = slice_by_eta(data, lo, hi)
        if "QQ_Events" not in sub or "GG_Events" not in sub:
            log("  (no jets in this bin — skipped)")
            rows.append({"bin": f"[{lo:+.1f},{hi:+.1f})", "n_qq": 0, "n_gg": 0,
                         "auc_nsd": np.nan, "sep_nsd": np.nan,
                         "auc_nsj": np.nan, "sep_nsj": np.nan})
            continue
        nqq = sub["QQ_Events"]["nsd"].size
        ngg = sub["GG_Events"]["nsd"].size
        if nqq < args.min_jets or ngg < args.min_jets:
            log(f"  N_QQ={nqq}  N_GG={ngg}  — below min-jets={args.min_jets}, skipping PDF")
            rows.append({"bin": f"[{lo:+.1f},{hi:+.1f})", "n_qq": nqq, "n_gg": ngg,
                         "auc_nsd": np.nan, "sep_nsd": np.nan,
                         "auc_nsj": np.nan, "sep_nsj": np.nan})
            continue

        pdf_name = f"roc_eta_{lo:+.1f}_{hi:+.1f}.pdf".replace("+", "p").replace("-", "m")
        st = {}
        plot_roc_single(sub,
                        f"{sample_name}  —  $\\eta \\in [{lo:.1f}, {hi:.1f})$",
                        pdf_name, st)
        log(f"  N_QQ={st['n_qq']:,}  N_GG={st['n_gg']:,}")
        log(f"  n_SD       : AUC = {st['nsd']['auc']:.4f}   sep = {st['nsd']['sep']:+.3f}")
        log(f"  n_subjets  : AUC = {st['nsubjets']['auc']:.4f}   sep = {st['nsubjets']['sep']:+.3f}")
        log(f"  -> wrote {pdf_name}")
        panels.append((f"$\\eta \\in [{lo:.1f}, {hi:.1f})$", sub))
        rows.append({
            "bin": f"[{lo:+.1f},{hi:+.1f})",
            "n_qq": st["n_qq"], "n_gg": st["n_gg"],
            "auc_nsd": st["nsd"]["auc"], "sep_nsd": st["nsd"]["sep"],
            "auc_nsj": st["nsubjets"]["auc"], "sep_nsj": st["nsubjets"]["sep"],
        })

    # --- combined grid -----------------------------------------------------
    if panels:
        plot_grid(panels, "roc_all_bins.pdf", sample_name)
        log("\n-> wrote roc_all_bins.pdf (grid view)")

    # --- summary table -----------------------------------------------------
    log("\n" + "=" * 78)
    log("SUMMARY")
    log("=" * 78)
    log(f"{'bin':<16}{'N_QQ':>10}{'N_GG':>10}"
        f"{'AUC_nSD':>12}{'sep_nSD':>10}"
        f"{'AUC_nSJ':>12}{'sep_nSJ':>10}")
    log("-" * 80)
    for r in rows:
        if np.isnan(r["auc_nsd"]):
            log(f"{r['bin']:<16}{r['n_qq']:>10}{r['n_gg']:>10}"
                f"{'  —':>12}{'  —':>10}{'  —':>12}{'  —':>10}")
        else:
            log(f"{r['bin']:<16}{r['n_qq']:>10}{r['n_gg']:>10}"
                f"{r['auc_nsd']:>12.4f}{r['sep_nsd']:>+10.3f}"
                f"{r['auc_nsj']:>12.4f}{r['sep_nsj']:>+10.3f}")

    log_path.write_text("\n".join(log_lines) + "\n")
    print(f"\nLog written to {log_path.resolve()}")


if __name__ == "__main__":
    main()
