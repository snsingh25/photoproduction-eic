#!/usr/bin/env python3
"""Paper figures for the EIC_Observables ROC comparison.

Produces the two figures the paper needs, both on identical jets
(anti-kT, R=1.0, exclusive 2-jet, ET > 10 GeV) at all four energies:

  fig_roc_panels.pdf   2x2 square panels, one per sqrt(s), each overlaying
                       the three discriminants. Full text width.
  fig_auc_vs_energy.pdf  AUC vs sqrt(s), one line per discriminant.
                       Single column, square.

The ROC definition is imported from roc_cross_energy_with_psi so the paper
figures and the scan tables cannot drift apart.

Usage:
    python plots/softdrop/paper_roc_figures.py
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "style"))
import paper_style as ps  # noqa: E402


def _load_roc_module():
    """Import the ROC definitions from the scan script (single source of truth)."""
    spec = importlib.util.spec_from_file_location(
        "roc_cross_energy_with_psi", HERE / "roc_cross_energy_with_psi.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rc = _load_roc_module()

# anti-kT, EtMin=10 for every sample -> "identical jets" is literally true
SAMPLES = [
    (64,  "eic64_antikt_dijets"),
    (105, "eic105_antikt_dijets"),
    (141, "eic141_antikt_dijets"),
    (300, "hera300_antikt_dijets"),
]

# Psi is quark-larger, so its ROC is evaluated on -Psi to keep one
# "signal-larger" convention across all three observables.
KIND = {"jet_psi03": ("continuous", True),
        "jet_nsd": ("integer", False),
        "jet_nsubjets": ("integer", False)}


def compute(data_jets):
    """Return {sqrts: {obs: {eff_sig, eff_bkg, auc, n_qq, n_gg}}}."""
    out = {}
    for sqrts, sample in SAMPLES:
        rf = rc.find_root(data_jets, sample)
        if rf is None:
            print(f"  [skip] {sample}: no ROOT file")
            continue
        d = rc.load(rf)
        if d is None or "QQ_Events" not in d or "GG_Events" not in d:
            print(f"  [skip] {sample}: missing branches")
            continue
        print(f"  [load] {sqrts:>3} GeV  {Path(rf).name}")
        qq, gg = d["QQ_Events"], d["GG_Events"]
        per_obs = {}
        for obs in ps.OBS_ORDER:
            kind, sigsmall = KIND[obs]
            eff_s, eff_b, auc = rc.roc(gg[obs], qq[obs], kind, sigsmall)
            per_obs[obs] = {"eff_sig": eff_s, "eff_bkg": eff_b, "auc": auc,
                            "n_qq": qq[obs].size, "n_gg": gg[obs].size}
        out[sqrts] = per_obs
    return out


def fig_roc_panels(results, out_path):
    """2x2 square panels, one per energy, three discriminants each."""
    fig, axes = plt.subplots(2, 2, figsize=(ps.FULL_W, ps.FULL_W * 0.99))

    for idx, (sqrts, _sample) in enumerate(SAMPLES):
        ax = axes[idx // 2][idx % 2]
        ps.square(ax)
        ps.chance_line(ax)
        res = results.get(sqrts)
        if res is None:
            ax.set_visible(False)
            continue

        for obs in ps.OBS_ORDER:
            r = res[obs]
            dashes = ps.OBS_DASH[obs]
            line, = ax.plot(r["eff_bkg"], r["eff_sig"],
                            color=ps.OBS_COLOR[obs], linewidth=1.4,
                            solid_capstyle="round",
                            label=ps.auc_legend_label(obs, r["auc"]))
            if dashes[0] is not None:
                line.set_dashes(dashes)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(0, 1.01, 0.2))
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.set_title(ps.energy_label(sqrts), pad=6)

        # only outer panels carry axis labels
        if idx // 2 == 1:
            ax.set_xlabel(r"$\varepsilon_{\mathrm{QQ}}$")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel(r"$\varepsilon_{\mathrm{GG}}$")
        else:
            ax.set_yticklabels([])

        leg = ax.legend(loc="lower right",
                        title=r"\makebox[1.6cm][l]{}AUC",
                        title_fontsize=7, borderpad=0.4)
        leg.get_title().set_color(ps.MUTED)

    fig.subplots_adjust(wspace=0.06, hspace=0.14)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_auc_vs_energy(results, out_path):
    """AUC vs sqrt(s), one line per discriminant. Single column, square."""
    fig, ax = plt.subplots(figsize=(ps.COL_W, ps.COL_W))
    ps.square(ax)

    energies = [s for s, _ in SAMPLES if s in results]

    ax.axhline(0.5, color=ps.MUTED, linestyle=(0, (2.5, 2.5)),
               linewidth=0.7, zorder=0)
    ax.text(315, 0.506, "no discrimination", ha="right", va="bottom",
            fontsize=6.5, color=ps.MUTED)

    for obs in ps.OBS_ORDER:
        aucs = [results[s][obs]["auc"] for s in energies]
        dashes = ps.OBS_DASH[obs]
        line, = ax.plot(energies, aucs, color=ps.OBS_COLOR[obs],
                        linewidth=1.4, zorder=2, label=ps.OBS_LABEL[obs])
        if dashes[0] is not None:
            line.set_dashes(dashes)
        # markers drawn separately so the dash pattern never eats them
        ax.plot(energies, aucs, linestyle="none", marker="o",
                markersize=3.2, color=ps.OBS_COLOR[obs], zorder=3)

    ax.set_xlabel(r"$\sqrt{s}$ [GeV]")
    ax.set_ylabel(r"ROC AUC, gluon vs quark")
    ax.set_xlim(44, 320)
    ax.set_ylim(0.48, 0.82)
    ax.set_xticks(energies)
    ax.set_yticks(np.arange(0.50, 0.81, 0.05))
    ax.legend(loc="lower left", bbox_to_anchor=(0.015, 0.10))

    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


ETA_BINS = [(-1.0, 0.0), (0.0, 1.0), (1.0, 1.5), (1.5, 2.0)]


def compute_eta(data_jets, min_jets=50):
    """Per-eta-bin AUC, plus the inclusive value, for every sample."""
    out = {}
    for sqrts, sample in SAMPLES:
        rf = rc.find_root(data_jets, sample)
        if rf is None:
            continue
        d = rc.load(rf)
        if d is None:
            continue
        per = {"incl": {}, "bins": {}}
        qq, gg = d["QQ_Events"], d["GG_Events"]
        for obs in ps.OBS_ORDER:
            kind, sigsmall = KIND[obs]
            per["incl"][obs] = rc.roc(gg[obs], qq[obs], kind, sigsmall)[2]
        for lo, hi in ETA_BINS:
            sub = rc.slice_eta(d, lo, hi)
            if "QQ_Events" not in sub or "GG_Events" not in sub:
                continue
            q, g = sub["QQ_Events"], sub["GG_Events"]
            if q["jet_eta"].size < min_jets or g["jet_eta"].size < min_jets:
                continue
            per["bins"][(lo, hi)] = {
                obs: rc.roc(g[obs], q[obs], *KIND[obs])[2] for obs in ps.OBS_ORDER}
        out[sqrts] = per
    return out


def fig_auc_vs_eta(eta_results, out_path):
    """2x2 square panels: per-bin AUC vs eta, inclusive shown as a flat line."""
    fig, axes = plt.subplots(2, 2, figsize=(ps.FULL_W, ps.FULL_W * 0.99))
    centres = [0.5 * (lo + hi) for lo, hi in ETA_BINS]

    for idx, (sqrts, _s) in enumerate(SAMPLES):
        ax = axes[idx // 2][idx % 2]
        ps.square(ax)
        res = eta_results.get(sqrts)
        if res is None:
            ax.set_visible(False)
            continue

        for obs in ps.OBS_ORDER:
            col = ps.OBS_COLOR[obs]
            xs = [c for c, b in zip(centres, ETA_BINS) if b in res["bins"]]
            ys = [res["bins"][b][obs] for b in ETA_BINS if b in res["bins"]]
            # inclusive value as a flat reference line
            ax.axhline(res["incl"][obs], color=col, linewidth=0.8,
                       linestyle=(0, (1.0, 2.0)), alpha=0.85, zorder=1)
            line, = ax.plot(xs, ys, color=col, linewidth=1.4, zorder=2,
                            label=ps.OBS_LABEL[obs])
            dashes = ps.OBS_DASH[obs]
            if dashes[0] is not None:
                line.set_dashes(dashes)
            ax.plot(xs, ys, linestyle="none", marker="o", markersize=3.2,
                    color=col, zorder=3)

        ax.set_xlim(-1.25, 2.25)
        ax.set_ylim(0.50, 0.82)
        ax.set_xticks(centres)
        ax.set_xticklabels([rf"${lo:g}$--${hi:g}$" for lo, hi in ETA_BINS])
        ax.set_yticks(np.arange(0.50, 0.81, 0.05))
        ax.set_title(ps.energy_label(sqrts), pad=6)
        if idx // 2 == 1:
            ax.set_xlabel(r"jet $\eta$")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel(r"ROC AUC, gluon vs quark")
        else:
            ax.set_yticklabels([])
        if idx == 0:
            ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02))

    fig.subplots_adjust(wspace=0.06, hspace=0.14)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


# Binning-width scan. Uniform bins over one fixed range, so the 1-bin case is
# the genuine pooled limit of the same jets and the comparison is like-for-like.
ETA_SCAN_RANGE = (-1.0, 2.0)
ETA_SCAN_NBINS = [1, 2, 3, 6, 12]


def compute_eta_scan(data_jets, min_jets=100):
    """Jet-weighted mean per-bin AUC as a function of eta bin width."""
    lo0, hi0 = ETA_SCAN_RANGE
    out = {}
    for sqrts, sample in SAMPLES:
        rf = rc.find_root(data_jets, sample)
        if rf is None:
            continue
        d = rc.load(rf)
        if d is None:
            continue
        per_scheme = {}
        for nb in ETA_SCAN_NBINS:
            edges = np.linspace(lo0, hi0, nb + 1)
            width = (hi0 - lo0) / nb
            acc = {obs: [] for obs in ps.OBS_ORDER}
            weights, kept, dropped = [], 0, 0
            for lo, hi in zip(edges[:-1], edges[1:]):
                sub = rc.slice_eta(d, lo, hi)
                if "QQ_Events" not in sub or "GG_Events" not in sub:
                    dropped += 1
                    continue
                q, g = sub["QQ_Events"], sub["GG_Events"]
                nq, ng = q["jet_eta"].size, g["jet_eta"].size
                if nq < min_jets or ng < min_jets:
                    dropped += 1
                    continue
                for obs in ps.OBS_ORDER:
                    acc[obs].append(rc.roc(g[obs], q[obs], *KIND[obs])[2])
                weights.append(nq + ng)
                kept += 1
            if not weights:
                continue
            w = np.array(weights, dtype=float)
            per_scheme[nb] = {
                "width": width,
                "auc": {obs: float(np.average(acc[obs], weights=w))
                        for obs in ps.OBS_ORDER},
                "n_bins_kept": kept,
                "n_bins_dropped": dropped,
                "n_jets": int(w.sum()),
            }
            if dropped:
                print(f"  [note] sqrt(s)={sqrts}, {nb} bins: dropped {dropped} "
                      f"bin(s) below {min_jets} jets/class")
        out[sqrts] = per_scheme
    return out


def fig_auc_vs_binwidth(scan, out_path):
    """AUC against eta bin width. Flat = pooling-immune, rising = contaminated."""
    fig, axes = plt.subplots(2, 2, figsize=(ps.FULL_W, ps.FULL_W * 0.99))
    for idx, (sqrts, _s) in enumerate(SAMPLES):
        ax = axes[idx // 2][idx % 2]
        ps.square(ax)
        res = scan.get(sqrts)
        if not res:
            ax.set_visible(False)
            continue
        widths = [res[nb]["width"] for nb in sorted(res, reverse=True)]
        for obs in ps.OBS_ORDER:
            ys = [res[nb]["auc"][obs] for nb in sorted(res, reverse=True)]
            col = ps.OBS_COLOR[obs]
            line, = ax.plot(widths, ys, color=col, linewidth=1.4, zorder=2,
                            label=ps.OBS_LABEL[obs])
            dashes = ps.OBS_DASH[obs]
            if dashes[0] is not None:
                line.set_dashes(dashes)
            ax.plot(widths, ys, linestyle="none", marker="o", markersize=3.2,
                    color=col, zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(0.21, 3.5)
        ax.set_ylim(0.61, 0.80)
        ax.set_xticks([0.25, 0.5, 1.0, 1.5, 3.0])
        ax.set_xticklabels(["0.25", "0.5", "1", "1.5", "3"])
        ax.set_yticks(np.arange(0.62, 0.80, 0.04))
        ax.set_title(ps.energy_label(sqrts), pad=6)
        if idx // 2 == 1:
            ax.set_xlabel(r"$\eta$ bin width")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel(r"ROC AUC, gluon vs quark")
        else:
            ax.set_yticklabels([])
        if idx == 0:
            ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02))
    fig.subplots_adjust(wspace=0.06, hspace=0.14)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


def write_scan_table(scan, out_path):
    lo0, hi0 = ETA_SCAN_RANGE
    lines = [f"Jet-weighted mean AUC vs eta bin width, uniform bins over "
             f"[{lo0:g},{hi0:g}]",
             "anti-kT, R=1.0, ET>10 GeV. Width 3 = fully pooled.",
             "=" * 78,
             f"{'sqrt(s)':>8}{'bins':>6}{'width':>8}{'jets':>10}"
             f"{'Psi(0.3)':>10}{'n_SD':>9}{'n_subjets':>11}"]
    for sqrts, _ in SAMPLES:
        res = scan.get(sqrts)
        if not res:
            continue
        for nb in sorted(res):
            v = res[nb]
            lines.append(f"{sqrts:>8}{nb:>6}{v['width']:>8.2f}{v['n_jets']:>10}"
                         f"{v['auc']['jet_psi03']:>10.3f}"
                         f"{v['auc']['jet_nsd']:>9.3f}"
                         f"{v['auc']['jet_nsubjets']:>11.3f}")
        lines.append("-" * 78)
    Path(out_path).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def write_eta_table(eta_results, out_path):
    lines = ["Per-eta-bin AUC (dotted = inclusive). anti-kT, R=1.0, ET>10 GeV",
             "=" * 78,
             f"{'sqrt(s)':>8}  {'bin':<12}{'Psi(0.3)':>10}{'n_SD':>9}"
             f"{'n_subjets':>11}   {'nSJ-nSD':>9}"]
    for sqrts, _ in SAMPLES:
        r = eta_results.get(sqrts)
        if r is None:
            continue
        lines.append(f"{sqrts:>8}  {'inclusive':<12}"
                     f"{r['incl']['jet_psi03']:>10.3f}{r['incl']['jet_nsd']:>9.3f}"
                     f"{r['incl']['jet_nsubjets']:>11.3f}   "
                     f"{r['incl']['jet_nsubjets']-r['incl']['jet_nsd']:>+9.3f}")
        for b in ETA_BINS:
            if b not in r["bins"]:
                continue
            v = r["bins"][b]
            lines.append(f"{'':>8}  {f'[{b[0]:g},{b[1]:g})':<12}"
                         f"{v['jet_psi03']:>10.3f}{v['jet_nsd']:>9.3f}"
                         f"{v['jet_nsubjets']:>11.3f}   "
                         f"{v['jet_nsubjets']-v['jet_nsd']:>+9.3f}")
        lines.append("-" * 78)
    Path(out_path).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def write_table(results, out_path):
    lines = ["AUC — identical jets: anti-kT, R=1.0, exclusive 2-jet, ET>10 GeV",
             "=" * 74,
             f"{'sqrt(s)':>8}{'N_QQ':>10}{'N_GG':>10}"
             f"{'Psi(0.3)':>12}{'n_SD':>10}{'n_subjets':>12}"]
    for sqrts, _ in SAMPLES:
        if sqrts not in results:
            continue
        r = results[sqrts]
        lines.append(
            f"{sqrts:>8}{r['jet_psi03']['n_qq']:>10}{r['jet_psi03']['n_gg']:>10}"
            f"{r['jet_psi03']['auc']:>12.3f}{r['jet_nsd']['auc']:>10.3f}"
            f"{r['jet_nsubjets']['auc']:>12.3f}")
    lines.append("=" * 74)
    Path(out_path).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-jets", default="data-jets")
    ap.add_argument("--out-dir",
                    default=str(HERE / "output/paper_roc_antikt_etmin10"))
    args = ap.parse_args()

    data_jets = Path(args.data_jets).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ps.use()
    print(f"Scanning {data_jets}")
    results = compute(data_jets)
    if not results:
        print("no data")
        sys.exit(1)

    fig_roc_panels(results, out_dir / "fig_roc_panels.pdf")
    fig_auc_vs_energy(results, out_dir / "fig_auc_vs_energy.pdf")
    write_table(results, out_dir / "auc_table.log")

    eta_results = compute_eta(data_jets)
    fig_auc_vs_eta(eta_results, out_dir / "fig_auc_vs_eta.pdf")
    write_eta_table(eta_results, out_dir / "auc_eta_table.log")

    scan = compute_eta_scan(data_jets)
    fig_auc_vs_binwidth(scan, out_dir / "fig_auc_vs_binwidth.pdf")
    write_scan_table(scan, out_dir / "auc_binwidth_table.log")


if __name__ == "__main__":
    main()
