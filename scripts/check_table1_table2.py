#!/usr/bin/env python3
"""Cross-check Tables I and II of the paper against regenerated data.

Table I — total event counts per sample:
    EIC (10,100)   -> sqrt(s)=63.2   998,511 events
    EIC (10,275)   -> sqrt(s)=104.9  999,444 events
    EIC (18,275)   -> sqrt(s)=141    999,633 events
    HERA (27.5,820)-> sqrt(s)=300    999,881 events

Table II — subprocess contribution fractions. The `Contribution`
column is relative to the ALL-event count (Resolved + Direct = 100%).

Reads each sample's Combined_Events tree and tabulates:
    - total events
    - Resolved / Direct event counts and fractions
    - Per-process-code fractions (within total)

Usage:
    python scripts/check_table1_table2.py [<root1> <root2> ...]
    (default: data/<sample>/*.root for hera300, eic64, eic105, eic141)
"""

import argparse
import glob
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import uproot


SAMPLES = [
    ("eic64",   "EIC 64 GeV (10, 100)",    63.2,  998_511),
    ("eic105",  "EIC 105 GeV (10, 275)",   104.9, 999_444),
    ("eic141",  "EIC 141 GeV (18, 275)",   141.0, 999_633),
    ("hera300", "HERA 300 GeV (27.5, 820)", 300.0, 999_881),
]

# Paper Table II reference values: percentages of total events, broken
# down by subprocess. Keys are (sample_key, process_code) -> percent.
# Only the subprocesses explicitly listed in the paper are included;
# unspecified codes contribute to the "Minor (<1% each)" residual.
PAPER_TABLE2 = {
    # EIC 64 GeV
    ("eic64",   "Resolved"): 50.2,
    ("eic64",   113): 27.09,
    ("eic64",   114): 16.47,
    ("eic64",   111): 5.42,
    ("eic64",   "Direct"):   49.8,
    ("eic64",   274): 26.73,
    ("eic64",   271): 14.29,
    # EIC 105 GeV
    ("eic105",  "Resolved"): 62.6,
    ("eic105",  113): 35.26,
    ("eic105",  114): 15.81,
    ("eic105",  111): 10.33,
    ("eic105",  "Direct"):   37.4,
    ("eic105",  274): 16.24,
    ("eic105",  271): 12.91,
    # EIC 141 GeV
    ("eic141",  "Resolved"): 61.4,
    ("eic141",  113): 34.66,
    ("eic141",  114): 16.62,
    ("eic141",  111): 8.92,
    ("eic141",  "Direct"):   38.6,
    ("eic141",  274): 16.80,
    ("eic141",  271): 12.74,
    # HERA 300 GeV
    ("hera300", "Resolved"): 74.4,
    ("hera300", 113): 42.36,
    ("hera300", 111): 17.37,
    ("hera300", 114): 13.48,
    ("hera300", "Direct"):   25.6,
    ("hera300", 271): 10.15,
    ("hera300", 274): 8.10,
}

TOLERANCE_PP = 0.20   # flag if |data - paper| > this many percentage points

# Process codes known to be in the "dominant" rows of Table II per sample.
# Used just for display ordering; we print all non-zero codes anyway.
PROCESS_NAMES = {
    111: "g g -> g g",
    112: "g g -> q qbar",
    113: "q g -> q g",
    114: "q q' -> q q'",
    115: "q qbar -> g g",
    116: "q qbar -> q' qbar'",
    121: "g g -> c cbar",
    122: "q qbar -> c cbar",
    123: "g g -> b bbar",
    124: "q qbar -> b bbar",
    271: "g gamma -> q qbar (uds)",
    272: "g gamma -> c cbar",
    273: "g gamma -> b bbar",
    274: "q gamma -> q g",
    281: "gamma g -> q qbar (uds)",
    282: "gamma g -> c cbar",
    283: "gamma g -> b bbar",
    284: "gamma q -> q g",
}

RESOLVED_CODES = set(range(111, 125))
DIRECT_CODES   = set(range(271, 285))


def find_root(sample_dir):
    hits = sorted(glob.glob(f"{sample_dir}/*.root"))
    return hits[0] if hits else None


def audit_sample(root_path, display_name, sqrts, paper_N, sample_key=None):
    print("=" * 86)
    print(f"Sample: {display_name}   (sqrt(s) = {sqrts} GeV)")
    print(f"  file: {root_path}")
    print("=" * 86)

    with uproot.open(root_path) as f:
        if "Combined_Events/Combined_Events" not in f:
            print("  Combined_Events tree missing — skipping")
            return 0
        t = f["Combined_Events/Combined_Events"]
        codes = t["processCode"].array(library="np")

    N = len(codes)
    resolved_N = int(((codes >= 111) & (codes <= 124)).sum())
    direct_N   = int(((codes >= 271) & (codes <= 284)).sum())
    other_N    = N - resolved_N - direct_N

    fails = 0   # incremented per |delta| > TOLERANCE_PP

    # --- Table I check ---
    print(f"\n[Table I] Combined events:")
    dN = N - paper_N
    pct = 100 * dN / paper_N if paper_N else 0
    print(f"  Paper:    {paper_N:>10,}")
    print(f"  Data:     {N:>10,}")
    print(f"  Diff:     {dN:>+10,}  ({pct:+.3f}%)")
    status = "PASS" if abs(dN) < 0.01 * paper_N else "CHECK"
    print(f"  Status:   {status}  (threshold: |diff| < 1% of paper)")

    # --- Table II check ---
    print(f"\n[Table II] Photon-interaction decomposition (fractions of total):")

    def vs_paper(cat_or_code, data_pct):
        """Return a ' (paper X.XX, Δ +0.YY pp) [FLAG]' string if a paper
        value is known for this sample_key, else ''."""
        if sample_key is None:
            return ""
        paper_pct = PAPER_TABLE2.get((sample_key, cat_or_code))
        if paper_pct is None:
            return ""
        delta = data_pct - paper_pct
        marker = "" if abs(delta) <= TOLERANCE_PP else "  ** FLAG **"
        return f"  (paper {paper_pct:5.2f}%, Δ {delta:+.2f} pp{marker})"

    res_pct = 100 * resolved_N / N
    dir_pct = 100 * direct_N   / N
    pr = vs_paper("Resolved", res_pct)
    pd = vs_paper("Direct",   dir_pct)
    if pr and "FLAG" in pr: fails += 1
    if pd and "FLAG" in pd: fails += 1
    print(f"  Resolved  : {resolved_N:>8,}  ({res_pct:6.2f}%){pr}")
    print(f"  Direct    : {direct_N:>8,}  ({dir_pct:6.2f}%){pd}")
    if other_N:
        print(f"  Uncategorised: {other_N:>8,}  ({100*other_N/N:6.2f}%)")

    # Per-process-code table, sorted by count within each category
    counts = Counter(codes.tolist())
    print(f"\n  {'Category':<10}{'Code':<6}{'Process':<30}{'N':>10}{'% of total':>14}  vs paper")
    print("  " + "-" * 88)
    for label, code_set in [("Resolved", sorted(RESOLVED_CODES)),
                            ("Direct",   sorted(DIRECT_CODES))]:
        sub = [(code, counts.get(code, 0)) for code in code_set if counts.get(code, 0) > 0]
        sub.sort(key=lambda x: -x[1])
        for code, n in sub:
            pct = 100 * n / N
            pp = vs_paper(code, pct)
            if pp and "FLAG" in pp: fails += 1
            print(f"  {label:<10}{code:<6}{PROCESS_NAMES.get(code, '?'):<30}{n:>10,}{pct:>13.2f}%{pp}")
    print()
    return fails


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_files", nargs="*",
                    help="If empty, uses data/<sample>/*.root for each known sample")
    args = ap.parse_args()

    total_fails = 0
    if args.root_files:
        for rf in args.root_files:
            name = Path(rf).stem
            matched_key = None
            for key, disp, s, n in SAMPLES:
                if key in name.lower():
                    matched_key = key
                    total_fails += audit_sample(rf, disp, s, n, sample_key=key) or 0
                    break
            if matched_key is None:
                total_fails += audit_sample(rf, name, 0.0, 0) or 0
    else:
        for sample, display, sqrts, paper_N in SAMPLES:
            sdir = f"data/{sample}"
            rf = find_root(sdir)
            if rf is None:
                print(f"[skip] {sample}: no ROOT under {sdir}")
                continue
            total_fails += audit_sample(rf, display, sqrts, paper_N, sample_key=sample) or 0

    print("=" * 86)
    if total_fails == 0:
        print(f"OVERALL: all paper Table I / II values match within {TOLERANCE_PP} pp tolerance.")
    else:
        print(f"OVERALL: {total_fails} value(s) outside {TOLERANCE_PP} pp tolerance — see ** FLAG ** above.")
    print("=" * 86)


if __name__ == "__main__":
    main()
