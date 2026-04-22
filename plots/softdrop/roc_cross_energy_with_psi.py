#!/usr/bin/env python3
"""3-way ROC master plot (Psi(r=0.3) + n_SD + n_subjets) across 4 samples.

Same convention as roc_cross_energy.py but extended to include the
paper's primary integrated-jet-shape discriminant Psi(r=0.3). Psi is
quark-larger so the ROC is evaluated on (-Psi) to keep a uniform
"signal-larger" convention.

Default samples are the paper-config dijets:
    HERA 300 : kT     + dijets + ET>17
    EIC 141  : anti-kT+ dijets + ET>10
    EIC 105  : anti-kT+ dijets + ET>10
    EIC 64   : anti-kT+ dijets + ET>10

Usage:
    python plots/softdrop/roc_cross_energy_with_psi.py
    (override default --samples / --out-dir similar to roc_cross_energy.py)
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt


DEFAULT_SAMPLES = [
    ("eic64_antikt_dijets",   "EIC 64 (anti-kT dijets)",   64,  "C0"),
    ("eic105_antikt_dijets",  "EIC 105 (anti-kT dijets)",  105, "C2"),
    ("eic141_antikt_dijets",  "EIC 141 (anti-kT dijets)",  141, "C3"),
    ("hera300_kt_dijets",     "HERA 300 (kT dijets)",      300, "C4"),
]

OBS_META = [
    ("jet_psi03",    r"$\Psi(r=0.3)$",          "C2", "continuous", True),
    ("jet_nsd",      r"$n_{\mathrm{SD}}$",      "C1", "integer",    False),
    ("jet_nsubjets", r"$n_{\mathrm{subjets}}$", "C4", "integer",    False),
]   # last field: signal_smaller? (True for Psi, since Psi is quark-larger)


def find_root(data_jets_dir, sample):
    for pat in ("dijets_*.root", "alljets_*.root", "*.root"):
        hits = sorted(glob.glob(f"{data_jets_dir}/{sample}/{pat}"))
        if hits:
            return hits[0]
    return None


def load(root_path):
    data = {}
    with uproot.open(root_path) as f:
        for cat in ("QQ_Events", "GG_Events"):
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            t = f[key]
            out = {}
            for br in ("jet_psi03", "jet_nsd", "jet_nsubjets", "jet_eta"):
                if br not in t.keys():
                    return None
                a = t[br].array(library="np")
                out[br] = np.concatenate([np.asarray(x) for x in a if len(x) > 0])
            data[cat] = out
    return data


def slice_eta(d, lo, hi):
    out = {}
    for cat, dd in d.items():
        m = (dd["jet_eta"] >= lo) & (dd["jet_eta"] < hi)
        if m.sum() == 0:
            continue
        out[cat] = {k: v[m] for k, v in dd.items()}
    return out


def roc(sig_arr, bkg_arr, kind, signal_smaller):
    """Return (eff_sig, eff_bkg, auc) for one observable."""
    if signal_smaller:
        sig_arr = -sig_arr
        bkg_arr = -bkg_arr
    if kind == "integer":
        vmin = int(min(sig_arr.min(), bkg_arr.min()))
        vmax = int(max(sig_arr.max(), bkg_arr.max()))
        grid = np.arange(vmin, vmax + 2)
    else:
        lo = float(min(sig_arr.min(), bkg_arr.min()))
        hi = float(max(sig_arr.max(), bkg_arr.max()))
        grid = np.linspace(lo - 1e-6, hi + 1e-6, 400)
    eff_sig = np.array([(sig_arr >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg_arr >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    auc = np.trapezoid(eff_sig[order], eff_bkg[order])
    return eff_sig[order], eff_bkg[order], auc


def collect(sample_data, bins, min_jets=20):
    """For each sample, compute per-bin ROC curves for each observable."""
    labels = ["inclusive"] + [f"$\\eta \\in [{lo:g}, {hi:g})$" for lo, hi in bins]
    results = {}
    for sname, d in sample_data.items():
        per_sample = {}
        per_sample[labels[0]] = _obs_block(d, min_jets)
        for (lo, hi), blabel in zip(bins, labels[1:]):
            sub = slice_eta(d, lo, hi)
            per_sample[blabel] = _obs_block(sub, min_jets)
        results[sname] = per_sample
    return results, labels


def _obs_block(sub, min_jets):
    if "QQ_Events" not in sub or "GG_Events" not in sub:
        return None
    qq, gg = sub["QQ_Events"], sub["GG_Events"]
    if qq["jet_psi03"].size < min_jets or gg["jet_psi03"].size < min_jets:
        return None
    out = {}
    for br, _lbl, _c, kind, sigsmall in OBS_META:
        eff_s, eff_b, auc = roc(gg[br], qq[br], kind, sigsmall)
        out[br] = {"eff_sig": eff_s, "eff_bkg": eff_b, "auc": auc,
                   "n_qq": qq[br].size, "n_gg": gg[br].size}
    return out


def plot_row(ax_list, bin_labels, obs_key, obs_display, obs_kind, color_per_sample, results):
    """One row: one panel per bin, all samples overlaid for a single observable."""
    for ax, blabel in zip(ax_list, bin_labels):
        for sname, slabel, col in color_per_sample:
            r = results.get(sname, {}).get(blabel)
            if r is None or obs_key not in r:
                continue
            ax.plot(r[obs_key]["eff_bkg"], r[obs_key]["eff_sig"],
                    "o-" if obs_kind == "integer" else "-",
                    color=col, linewidth=1.8,
                    markersize=4 if obs_kind == "integer" else 0,
                    label=f"{slabel}  AUC={r[obs_key]['auc']:.3f}")
        ax.plot([0, 1], [0, 1], "--", color="grey", alpha=0.4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
        ax.set_xlabel(r"$\varepsilon$(QQ)")
        ax.set_ylabel(r"$\varepsilon$(GG)")
        ax.set_title(f"{obs_display}  —  {blabel}", fontsize=11)
        ax.legend(frameon=False, loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bins", default="-1,0,1,1.5,2",
                    help="comma-separated eta bin edges (default paper bins)")
    ap.add_argument("--data-jets", default="data-jets",
                    help="root of per-sample output dirs")
    ap.add_argument("--out-dir", default="data-jets/cross_energy_paperconfig",
                    help="where to write output PDFs")
    ap.add_argument("--samples",
                    help=("semicolon-separated sample specs "
                          "'dir,label,sqrts,color;...'. If omitted, uses "
                          "the 4 paper-config dijet samples."))
    ap.add_argument("--min-jets", type=int, default=20)
    args = ap.parse_args()

    edges = [float(x) for x in args.bins.split(",")]
    bins = list(zip(edges[:-1], edges[1:]))

    samples = list(DEFAULT_SAMPLES)
    if args.samples:
        samples = []
        for spec in args.samples.split(";"):
            parts = [p.strip() for p in spec.split(",")]
            samples.append((parts[0], parts[1], int(parts[2]), parts[3]))

    data_jets = Path(args.data_jets).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scanning {data_jets}")
    sample_data = {}
    for sname, label, _s, _col in samples:
        rf = find_root(data_jets, sname)
        if rf is None:
            print(f"  [skip] {sname}: no ROOT"); continue
        d = load(rf)
        if d is None:
            print(f"  [skip] {sname}: missing branches"); continue
        sample_data[sname] = d
        print(f"  [load] {sname}: {rf}")

    if not sample_data:
        print("no data"); sys.exit(1)

    results, labels = collect(sample_data, bins, args.min_jets)
    color_per_sample = [(s[0], s[1], s[3]) for s in samples if s[0] in sample_data]

    n_panels = len(labels)

    # --- master grid: 3 rows (Psi, n_SD, n_subjets) x n_panels cols ---
    fig, axes = plt.subplots(3, n_panels, figsize=(4.5 * n_panels, 13))
    for row, (obs_key, obs_display, _col, kind, _sigsmall) in enumerate(OBS_META):
        plot_row(axes[row], labels, obs_key, obs_display, kind,
                 color_per_sample, results)
    fig.suptitle(
        r"ROC curves across $\eta$ bins and centre-of-mass energies  "
        r"(rows: $\Psi(r=0.3)$ [top] / $n_{SD}$ / $n_{subjets}$ [bottom])",
        fontsize=13, y=1.00)
    fig.tight_layout()
    out_path = out_dir / "roc_cross_energy_psi_nsd_nsubjets.pdf"
    fig.savefig(out_path); plt.close(fig)
    print(f"\nWrote {out_path}")

    # --- AUC summary table ---
    log_lines = ["AUC summary — paper-config dijet samples, paper eta bins",
                 "=" * 92,
                 f"{'Sample':<30}{'bin':<16}{'N_QQ':>8}{'N_GG':>8}"
                 f"{'Psi_AUC':>10}{'nSD_AUC':>10}{'nSJ_AUC':>10}"]
    for sname, slabel, _s, _col in samples:
        if sname not in results:
            continue
        for blabel in labels:
            r = results[sname].get(blabel)
            if r is None:
                log_lines.append(f"{slabel:<30}{blabel:<16}{'—':>8}{'—':>8}"
                                  f"{'—':>10}{'—':>10}{'—':>10}")
                continue
            log_lines.append(
                f"{slabel:<30}{blabel:<16}"
                f"{r['jet_psi03']['n_qq']:>8}{r['jet_psi03']['n_gg']:>8}"
                f"{r['jet_psi03']['auc']:>10.3f}"
                f"{r['jet_nsd']['auc']:>10.3f}"
                f"{r['jet_nsubjets']['auc']:>10.3f}")
        log_lines.append("-" * 92)
    log_path = out_dir / "roc_cross_energy_psi_nsd_nsubjets.log"
    log_path.write_text("\n".join(log_lines) + "\n")
    print("\n".join(log_lines))


if __name__ == "__main__":
    main()
