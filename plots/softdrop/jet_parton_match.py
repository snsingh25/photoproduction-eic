#!/usr/bin/env python3
"""Jet-to-parton Delta-R matching sanity check.

Reads the jetreco_softdrop output of a *validation* sample (one where
the input event tree had status-23 outgoing-parton branches, so that
`jet_parton_dR` and `jet_parton_pdgId` are populated per jet).

For each labelled category (QQ, GG, GQ) and separately for leading and
subleading jets, plots Delta-R to the nearest hard parton and reports:

  - Mean and median Delta-R.
  - Fraction of jets matched within Delta-R < 0.5 (= R/2 for anti-kT R=1).
  - Fraction unmatched or badly matched (Delta-R > 1.0).
  - Whether the nearest parton's PDG id is consistent with the label
    (e.g. for GG events, the nearest parton of BOTH leading and
    subleading should be a gluon).

Writes:
    jet_parton_match.pdf
    jet_parton_match_summary.log

Outputs into the current working directory.

Usage:
    python jet_parton_match.py <root_file>
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

MATCH_TIGHT = 0.5      # "good" match
MATCH_LOOSE = 1.0      # anything beyond this is effectively unmatched


def load(root_path):
    """Return dict: cat -> {'lead_dR','lead_pdg','sub_dR','sub_pdg',
    'lead_eta','sub_eta'} flat arrays."""
    data = {}
    with uproot.open(root_path) as f:
        for cat in CATEGORIES:
            key = f"{cat}/jets_{cat}"
            if key not in f:
                continue
            tree = f[key]
            if "jet_parton_dR" not in tree.keys():
                print(f"  [warn] {cat}: no jet_parton_dR branch — need a validation ROOT")
                continue
            jdr = tree["jet_parton_dR"].array(library="np")
            jpdg = tree["jet_parton_pdgId"].array(library="np")
            jeta = tree["jet_eta"].array(library="np")
            n_jets = tree["n_jets"].array(library="np")

            lead_dR, lead_pdg, lead_eta = [], [], []
            sub_dR,  sub_pdg,  sub_eta  = [], [], []
            for i in range(len(jdr)):
                n = n_jets[i]
                if n < 2 or len(jdr[i]) < 2:
                    continue
                lead_dR.append(jdr[i][0]);  lead_pdg.append(jpdg[i][0]);  lead_eta.append(jeta[i][0])
                sub_dR.append(jdr[i][1]);   sub_pdg.append(jpdg[i][1]);   sub_eta.append(jeta[i][1])
            if not lead_dR:
                continue
            data[cat] = {
                "lead_dR":  np.asarray(lead_dR, dtype=float),
                "lead_pdg": np.asarray(lead_pdg, dtype=int),
                "lead_eta": np.asarray(lead_eta, dtype=float),
                "sub_dR":   np.asarray(sub_dR,  dtype=float),
                "sub_pdg":  np.asarray(sub_pdg, dtype=int),
                "sub_eta":  np.asarray(sub_eta, dtype=float),
            }
    return data


def pct(n, tot):
    return f"{100*n/tot:5.2f}%" if tot else "  —   "


def pdg_category(pdg_array, label):
    """Fraction of matched partons consistent with the label's expectation."""
    abs_pdg = np.abs(pdg_array)
    is_q = (abs_pdg >= 1) & (abs_pdg <= 6)
    is_g = abs_pdg == 21
    if label == "QQ_Events":
        return is_q.sum(), len(pdg_array)
    if label == "GG_Events":
        return is_g.sum(), len(pdg_array)
    if label == "GQ_Events":
        # GQ events have one quark + one gluon outgoing. We can't tell
        # which jet matches which at this stage (that's a separate test)
        # — just check the parton is a valid hard parton.
        return (is_q | is_g).sum(), len(pdg_array)
    return 0, len(pdg_array)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file")
    args = ap.parse_args()

    data = load(args.root_file)
    if not data:
        print(f"No usable jet_parton_dR data in {args.root_file}")
        sys.exit(1)

    sample_name = Path(args.root_file).stem
    log_lines = []
    def log(s=""):
        print(s); log_lines.append(s)

    log("=" * 94)
    log(f"Jet-parton Delta-R matching check  —  {sample_name}")
    log("=" * 94)
    log(f"Tight match: Delta-R < {MATCH_TIGHT}    Loose match: Delta-R < {MATCH_LOOSE}")
    log("")

    log(f"{'Category':<12}{'Jet':<10}{'N':>8}{'<dR>':>8}{'median':>9}"
        f"{'<0.5 (tight)':>14}{'<1.0 (loose)':>14}{'flavour OK':>13}")
    log("-" * 94)

    for cat in CATEGORIES:
        if cat not in data:
            continue
        d = data[cat]
        for jet_kind, drs, pdgs in [("leading", d["lead_dR"], d["lead_pdg"]),
                                    ("subleading", d["sub_dR"], d["sub_pdg"])]:
            n = drs.size
            good_tight = (drs < MATCH_TIGHT).sum()
            good_loose = (drs < MATCH_LOOSE).sum()
            # flavour-OK uses the tight-matched subset (only then does
            # 'nearest parton == event label' carry signal)
            tight_mask = drs < MATCH_TIGHT
            flav_ok, flav_tot = pdg_category(pdgs[tight_mask], cat)
            log(f"{LABELS[cat]:<12}{jet_kind:<10}{n:>8}"
                f"{drs.mean():>8.3f}{np.median(drs):>9.3f}"
                f"{good_tight:>6} {pct(good_tight, n):>8}"
                f"{good_loose:>6} {pct(good_loose, n):>8}"
                f"{flav_ok:>5} {pct(flav_ok, flav_tot):>8}")

    log("")
    log("Interpretation:")
    log("  - For anti-kT R=1 jets, tight-match fraction >~95% is expected:")
    log("    hard partons should be close to the jets they seeded.")
    log("  - >5% with Delta-R > 1.0 would indicate the LO-parton -> jet")
    log("    identification is unreliable (e.g. remnant contamination).")
    log("  - 'flavour OK' shows fraction of tightly-matched jets whose")
    log("    nearest parton is consistent with the event label.")

    # --- plots ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    for ax, jet_kind, key in [(axes[0], "Leading jet",    "lead_dR"),
                              (axes[1], "Subleading jet", "sub_dR")]:
        for cat in CATEGORIES:
            if cat not in data:
                continue
            drs = data[cat][key]
            ax.hist(drs, bins=80, range=(0, 2.0),
                    density=True, histtype="step", linewidth=2,
                    color=COLORS[cat], label=LABELS[cat])
        ax.axvline(MATCH_TIGHT, linestyle="--", color="grey", alpha=0.6,
                   label=f"tight cut ($\\Delta R$<{MATCH_TIGHT})")
        ax.axvline(MATCH_LOOSE, linestyle=":", color="grey", alpha=0.6,
                   label=f"loose cut ($\\Delta R$<{MATCH_LOOSE})")
        ax.set_xlabel(r"$\Delta R$(jet, nearest hard parton)")
        ax.set_ylabel("fraction of dijet events")
        ax.set_title(f"{sample_name} — {jet_kind}")
        ax.set_xlim(0, 2.0)
        ax.legend(frameon=False, loc="upper right")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    out_pdf = Path("jet_parton_match.pdf")
    fig.savefig(out_pdf); plt.close(fig)
    log(f"\nWrote {out_pdf.resolve()}")

    out_log = Path("jet_parton_match_summary.log")
    out_log.write_text("\n".join(log_lines) + "\n")
    print(f"Wrote {out_log.resolve()}")


if __name__ == "__main__":
    main()
