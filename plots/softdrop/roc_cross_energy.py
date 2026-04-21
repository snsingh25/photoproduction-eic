#!/usr/bin/env python3
"""One master PDF overlaying ROC curves across energies.

For each eta bin (plus the inclusive panel), plots ROC curves for every
sample on a single axis, one color per energy. Produces two PDFs at
`data-jets/`:

    roc_cross_energy_nSD.pdf          n_SD curves,       5 panels
    roc_cross_energy_nSubjets.pdf     n_subjets curves,  5 panels
    roc_cross_energy_master.pdf       both observables,  2x5 grid

Usage (from repo root):
    python plots/softdrop/roc_cross_energy.py
    python plots/softdrop/roc_cross_energy.py --bins=-1,0,1,2,3
    python plots/softdrop/roc_cross_energy.py --out-dir=data-jets
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt


DEFAULT_SAMPLES = [
    ("eic64_pT7",   "EIC 64 GeV",   64,  "C0"),
    ("eic105_pT7",  "EIC 105 GeV",  105, "C2"),
    ("eic141_pT7",  "EIC 141 GeV",  141, "C3"),
    ("hera300_pT7", "HERA 300 GeV", 300, "C4"),
]
# Gets overwritten by --samples if the user provides one.
SAMPLES = list(DEFAULT_SAMPLES)

CATEGORIES = ["QQ_Events", "GG_Events"]


# -----------------------------------------------------------------------------
def find_root(data_jets_dir, sample):
    # Match both alljets_*.root (default) and dijets_*.root (DIJET_ONLY mode).
    for pat in ("alljets_*.root", "dijets_*.root", "*.root"):
        hits = sorted(glob.glob(f"{data_jets_dir}/{sample}/{pat}"))
        if hits:
            return hits[0]
    return None


def load(root_path):
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


def slice_eta(data, lo, hi):
    out = {}
    for cat, d in data.items():
        m = (d["eta"] >= lo) & (d["eta"] < hi)
        if m.sum() == 0:
            continue
        out[cat] = {k: v[m] for k, v in d.items()}
    return out


def roc_integer(sig, bkg):
    vmin = int(min(sig.min(), bkg.min()))
    vmax = int(max(sig.max(), bkg.max()))
    grid = np.arange(vmin, vmax + 2)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    auc = np.trapezoid(eff_sig[order], eff_bkg[order])
    return eff_sig[order], eff_bkg[order], auc


# -----------------------------------------------------------------------------
def collect_curves(sample_data, bins, min_jets=20):
    """For each (sample, bin, observable), return {sample_key -> {bin_label -> (es, eb, auc, n_qq, n_gg)}}."""
    def _fmt(x):
        return f"{x:.0f}" if float(x).is_integer() else f"{x:.1f}"
    labels = ["inclusive"] + [
        f"$\\eta \\in [{_fmt(lo)}, {_fmt(hi)})$" for lo, hi in bins
    ]

    results = {}   # sample_key -> {bin_label -> {'nsd': (...), 'nsubjets': (...)}}
    for sname, label, _sqrts, _col in SAMPLES:
        d = sample_data.get(sname)
        if d is None:
            continue
        per_sample = {}
        per_sample[labels[0]] = make_roc_dict(d, min_jets)
        for (lo, hi), blabel in zip(bins, labels[1:]):
            sub = slice_eta(d, lo, hi)
            per_sample[blabel] = make_roc_dict(sub, min_jets)
        results[sname] = per_sample
    return results, labels


def make_roc_dict(sub, min_jets):
    """For a per-bin slice, return {'nsd': (es, eb, auc, n_qq, n_gg), 'nsubjets': ...}. None if stats too thin."""
    if "QQ_Events" not in sub or "GG_Events" not in sub:
        return None
    qq = sub["QQ_Events"]
    gg = sub["GG_Events"]
    if qq["nsd"].size < min_jets or gg["nsd"].size < min_jets:
        return None
    out = {}
    for obs in ("nsd", "nsubjets"):
        es, eb, auc = roc_integer(gg[obs], qq[obs])
        out[obs] = {"es": es, "eb": eb, "auc": auc,
                    "n_qq": qq[obs].size, "n_gg": gg[obs].size}
    return out


# -----------------------------------------------------------------------------
def plot_row(ax_list, bin_labels, obs_key, obs_display, results):
    """Plot one row: one panel per bin, each with 4-energy overlay for `obs_key`."""
    for ax, blabel in zip(ax_list, bin_labels):
        for sname, slabel, _sqrts, col in SAMPLES:
            r = results.get(sname, {}).get(blabel)
            if r is None or obs_key not in r:
                continue
            curve = r[obs_key]
            ax.plot(curve["eb"], curve["es"], "o-", color=col,
                    linewidth=1.8, markersize=4,
                    label=f"{slabel}  AUC={curve['auc']:.3f}")
        ax.plot([0, 1], [0, 1], "--", color="grey", alpha=0.4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
        ax.set_xlabel(r"$\varepsilon$(QQ)")
        ax.set_ylabel(r"$\varepsilon$(GG)")
        ax.set_title(f"{obs_display}  —  {blabel}", fontsize=11)
        ax.legend(frameon=False, loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)


# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bins", default="-1,0,1,2,3",
                    help="comma-separated eta bin edges (default: -1,0,1,2,3)")
    ap.add_argument("--data-jets", default="data-jets",
                    help="root of per-sample output dirs")
    ap.add_argument("--out-dir", default=None,
                    help="where to write output PDFs (default: same as --data-jets)")
    ap.add_argument("--min-jets", type=int, default=20)
    ap.add_argument("--samples",
                    help="semicolon-separated sample specs, e.g. "
                         "'eic64_antikt_dijets,EIC 64,64,C0;hera300_kt_dijets,HERA 300,300,C4'. "
                         "Each spec is <dirname>,<label>,<sqrts>,<color>. "
                         "If omitted, the default 4-sample list is used.")
    args = ap.parse_args()

    global SAMPLES
    if args.samples:
        SAMPLES = []
        for spec in args.samples.split(";"):
            parts = [p.strip() for p in spec.split(",")]
            if len(parts) != 4:
                raise SystemExit(f"Bad --samples entry '{spec}'; expected 4 fields")
            SAMPLES.append((parts[0], parts[1], int(parts[2]), parts[3]))

    bin_edges = [float(x) for x in args.bins.split(",")]
    bins = list(zip(bin_edges[:-1], bin_edges[1:]))

    data_jets = Path(args.data_jets).resolve()
    out_dir = Path(args.out_dir or args.data_jets).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scanning {data_jets} for samples: {[s[0] for s in SAMPLES]}")
    sample_data = {}
    for sname, label, _sqrts, _col in SAMPLES:
        root = find_root(data_jets, sname)
        if root is None:
            print(f"  [skip] {sname}: no alljets_*.root")
            continue
        print(f"  [load] {sname}: {root}")
        sample_data[sname] = load(root)

    if not sample_data:
        print("No samples loaded — nothing to plot.")
        sys.exit(1)

    results, labels = collect_curves(sample_data, bins, args.min_jets)

    n_panels = len(labels)

    # --- one-row PDFs (n_SD alone, n_subjets alone) ------------------------
    for obs_key, obs_display, outname in [
        ("nsd",      r"$n_{\mathrm{SD}}$",      "roc_cross_energy_nSD.pdf"),
        ("nsubjets", r"$n_{\mathrm{subjets}}$", "roc_cross_energy_nSubjets.pdf"),
    ]:
        fig, axes = plt.subplots(1, n_panels, figsize=(4.5 * n_panels, 4.5))
        axes = np.atleast_1d(axes)
        plot_row(axes, labels, obs_key, obs_display, results)
        fig.suptitle(f"ROC curves by $\\eta$ bin and centre-of-mass energy  —  {obs_display}",
                     fontsize=13, y=1.02)
        fig.tight_layout()
        fig.savefig(out_dir / outname)
        plt.close(fig)
        print(f"Wrote {out_dir / outname}")

    # --- master 2x5 grid --------------------------------------------------
    fig, axes = plt.subplots(2, n_panels, figsize=(4.5 * n_panels, 9))
    plot_row(axes[0], labels, "nsd",      r"$n_{\mathrm{SD}}$",      results)
    plot_row(axes[1], labels, "nsubjets", r"$n_{\mathrm{subjets}}$", results)
    fig.suptitle("ROC curves for $n_{SD}$ (top) and $n_{subjets}$ (bottom) across energies and $\\eta$ bins",
                 fontsize=13, y=1.00)
    fig.tight_layout()
    master_path = out_dir / "roc_cross_energy_master.pdf"
    fig.savefig(master_path)
    plt.close(fig)
    print(f"Wrote {master_path}")

    # --- compact text summary --------------------------------------------
    log_lines = ["Cross-energy ROC summary", "=" * 78,
                 f"{'Sample':<14}{'bin':<20}{'N_QQ':>8}{'N_GG':>8}"
                 f"{'AUC_nSD':>11}{'AUC_nSJ':>11}"]
    for sname, slabel, _sqrts, _col in SAMPLES:
        per_sample = results.get(sname, {})
        for blabel in labels:
            r = per_sample.get(blabel)
            if r is None:
                log_lines.append(f"{slabel:<14}{blabel:<20}{'—':>8}{'—':>8}{'—':>11}{'—':>11}")
                continue
            log_lines.append(f"{slabel:<14}{blabel:<20}"
                             f"{r['nsd']['n_qq']:>8}{r['nsd']['n_gg']:>8}"
                             f"{r['nsd']['auc']:>11.4f}{r['nsubjets']['auc']:>11.4f}")
        log_lines.append("-" * 78)
    log_path = out_dir / "roc_cross_energy_summary.log"
    log_path.write_text("\n".join(log_lines) + "\n")
    print(f"Wrote {log_path}")
    print("\n".join(log_lines))


if __name__ == "__main__":
    main()
