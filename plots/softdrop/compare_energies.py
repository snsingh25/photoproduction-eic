#!/usr/bin/env python3
"""Cross-energy comparison of n_SD vs n_subjets as q/g discriminants.

Reads ROOT outputs of jetreco_softdrop from hera/, eic141/, eic105/, eic64/
and shows how separation/AUC scale with sqrt(s).

Usage (run from softdrop/ dir):
    python compare_energies.py
"""

import os
import sys
from pathlib import Path
import numpy as np
import uproot
import matplotlib.pyplot as plt


_OUT_DIR = Path(__file__).resolve().parent / "output" / "cross_energy_paperconfig"
_OUT_DIR.mkdir(parents=True, exist_ok=True)


SAMPLES = [
    ("eic64",   "EIC 64 GeV",   64),
    ("eic105",  "EIC 105 GeV",  105),
    ("eic141",  "EIC 141 GeV",  141),
    ("hera",    "HERA 300 GeV", 300),
]

CATEGORIES = ["QQ_Events", "GG_Events"]


def load_sample(folder):
    """Find the alljets_*.root in folder and return dict cat -> {nsd, nsubjets, eta}."""
    import glob
    root_files = glob.glob(f"{folder}/alljets_*.root")
    if not root_files:
        return None
    d = {}
    with uproot.open(root_files[0]) as f:
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
            d[cat] = {
                "nsd":      np.concatenate([np.asarray(a) for a in nsd if len(a) > 0]),
                "nsubjets": np.concatenate([np.asarray(a) for a in nsj if len(a) > 0]),
                "eta":      np.concatenate([np.asarray(a) for a in eta if len(a) > 0]),
            }
    return d


def separation(a, b):
    return (b.mean() - a.mean()) / np.sqrt(a.var(ddof=1) + b.var(ddof=1))


def auc_roc(sig, bkg):
    grid = np.arange(min(sig.min(), bkg.min()), max(sig.max(), bkg.max()) + 2)
    eff_sig = np.array([(sig >= c).mean() for c in grid])
    eff_bkg = np.array([(bkg >= c).mean() for c in grid])
    order = np.argsort(eff_bkg)
    return np.trapezoid(eff_sig[order], eff_bkg[order])


# -----------------------------------------------------------------------------
rows = []
for folder, label, sqrts in SAMPLES:
    data = load_sample(folder)
    if data is None or "QQ_Events" not in data or "GG_Events" not in data:
        print(f"[skip] {folder} — missing ROOT or QQ/GG trees")
        continue
    qq, gg = data["QQ_Events"], data["GG_Events"]
    rows.append({
        "folder":     folder,
        "label":      label,
        "sqrts":      sqrts,
        "n_qq":       qq["nsd"].size,
        "n_gg":       gg["nsd"].size,
        "sep_nsd":    separation(qq["nsd"],      gg["nsd"]),
        "sep_nsj":    separation(qq["nsubjets"], gg["nsubjets"]),
        "auc_nsd":    auc_roc(gg["nsd"],      qq["nsd"]),
        "auc_nsj":    auc_roc(gg["nsubjets"], qq["nsubjets"]),
        "mean_nsd_qq": qq["nsd"].mean(),
        "mean_nsd_gg": gg["nsd"].mean(),
        "mean_nsj_qq": qq["nsubjets"].mean(),
        "mean_nsj_gg": gg["nsubjets"].mean(),
    })

if not rows:
    print("No samples loaded — did you run jetreco_softdrop in each subfolder?")
    sys.exit(1)

# --- summary table -----------------------------------------------------------
out_lines = []
def p(s):
    print(s); out_lines.append(s)

p("=" * 88)
p(f"{'Sample':<16}{'sqrt(s)':>10}{'N_QQ':>8}{'N_GG':>8}"
  f"{'sep_nSD':>10}{'sep_nSJ':>10}{'AUC_nSD':>10}{'AUC_nSJ':>10}")
p("-" * 88)
for r in rows:
    p(f"{r['label']:<16}{r['sqrts']:>10d}{r['n_qq']:>8d}{r['n_gg']:>8d}"
      f"{r['sep_nsd']:>10.3f}{r['sep_nsj']:>10.3f}"
      f"{r['auc_nsd']:>10.4f}{r['auc_nsj']:>10.4f}")

p("")
p("Means:")
p(f"{'Sample':<16}{'<n_SD>_QQ':>12}{'<n_SD>_GG':>12}"
  f"{'<n_subj>_QQ':>14}{'<n_subj>_GG':>14}")
for r in rows:
    p(f"{r['label']:<16}{r['mean_nsd_qq']:>12.2f}{r['mean_nsd_gg']:>12.2f}"
      f"{r['mean_nsj_qq']:>14.2f}{r['mean_nsj_gg']:>14.2f}")

_summary_path = _OUT_DIR / "compare_energies_summary.txt"
with open(_summary_path, "w") as fh:
    fh.write("\n".join(out_lines) + "\n")
p(f"\nWrote {_summary_path}")


# --- figure: separation and AUC vs sqrt(s) ----------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

sqrts = [r["sqrts"] for r in rows]
labels = [r["label"] for r in rows]
ax = axes[0]
ax.plot(sqrts, [r["sep_nsd"] for r in rows], "o-", color="C1", linewidth=2,
        markersize=9, label=r"$n_{\mathrm{SD}}$")
ax.plot(sqrts, [r["sep_nsj"] for r in rows], "s-", color="C4", linewidth=2,
        markersize=9, label=r"$n_{\mathrm{subjets}}$")
ax.set_xlabel(r"$\sqrt{s}$ [GeV]")
ax.set_ylabel(r"mean separation  $(\langle X\rangle_{GG}-\langle X\rangle_{QQ})/\sigma_{\mathrm{pooled}}$")
ax.set_title("Q/G separation vs. centre-of-mass energy")
ax.set_xticks(sqrts)
ax.set_xticklabels(labels, rotation=15)
ax.legend(frameon=False)
ax.grid(alpha=0.3)
ax.set_ylim(0, max([r["sep_nsd"] for r in rows] + [r["sep_nsj"] for r in rows]) * 1.3)

ax = axes[1]
ax.plot(sqrts, [r["auc_nsd"] for r in rows], "o-", color="C1", linewidth=2,
        markersize=9, label=r"$n_{\mathrm{SD}}$")
ax.plot(sqrts, [r["auc_nsj"] for r in rows], "s-", color="C4", linewidth=2,
        markersize=9, label=r"$n_{\mathrm{subjets}}$")
ax.axhline(0.5, linestyle="--", color="grey", alpha=0.7, label="random")
ax.set_xlabel(r"$\sqrt{s}$ [GeV]")
ax.set_ylabel("ROC AUC (GG vs QQ)")
ax.set_title("Q/G tagging AUC vs. centre-of-mass energy")
ax.set_xticks(sqrts)
ax.set_xticklabels(labels, rotation=15)
ax.legend(frameon=False)
ax.grid(alpha=0.3)
ax.set_ylim(0.45, 0.85)

fig.tight_layout()
_pdf_path = _OUT_DIR / "compare_energies.pdf"
fig.savefig(_pdf_path)
print(f"\nWrote {_pdf_path}")
