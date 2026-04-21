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


def audit_sample(root_path, display_name, sqrts, paper_N):
    print("=" * 86)
    print(f"Sample: {display_name}   (sqrt(s) = {sqrts} GeV)")
    print(f"  file: {root_path}")
    print("=" * 86)

    with uproot.open(root_path) as f:
        if "Combined_Events/Combined_Events" not in f:
            print("  Combined_Events tree missing — skipping")
            return
        t = f["Combined_Events/Combined_Events"]
        codes = t["processCode"].array(library="np")

    N = len(codes)
    resolved_N = int(((codes >= 111) & (codes <= 124)).sum())
    direct_N   = int(((codes >= 271) & (codes <= 284)).sum())
    other_N    = N - resolved_N - direct_N

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
    print(f"  Resolved  : {resolved_N:>8,}  ({100*resolved_N/N:6.2f}%)")
    print(f"  Direct    : {direct_N:>8,}  ({100*direct_N/N:6.2f}%)")
    if other_N:
        print(f"  Uncategorised: {other_N:>8,}  ({100*other_N/N:6.2f}%)")

    # Per-process-code table, sorted by count within each category
    counts = Counter(codes.tolist())
    print(f"\n  {'Category':<10}{'Code':<6}{'Process':<30}{'N':>10}{'% of total':>14}")
    print("  " + "-" * 68)
    for label, code_set in [("Resolved", sorted(RESOLVED_CODES)),
                            ("Direct",   sorted(DIRECT_CODES))]:
        # sort by count within category, descending
        sub = [(code, counts.get(code, 0)) for code in code_set if counts.get(code, 0) > 0]
        sub.sort(key=lambda x: -x[1])
        for code, n in sub:
            print(f"  {label:<10}{code:<6}{PROCESS_NAMES.get(code, '?'):<30}{n:>10,}{100*n/N:>13.2f}%")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_files", nargs="*",
                    help="If empty, uses data/<sample>/*.root for each known sample")
    args = ap.parse_args()

    if args.root_files:
        for rf in args.root_files:
            name = Path(rf).stem
            sqrts = 0.0
            paper_N = 0
            # Try to infer sample from filename
            for key, disp, s, n in SAMPLES:
                if key in name.lower():
                    sqrts, paper_N = s, n
                    audit_sample(rf, disp, sqrts, paper_N)
                    break
            else:
                audit_sample(rf, name, sqrts, paper_N)
    else:
        for sample, display, sqrts, paper_N in SAMPLES:
            sdir = f"data/{sample}"
            rf = find_root(sdir)
            if rf is None:
                print(f"[skip] {sample}: no ROOT under {sdir}")
                continue
            audit_sample(rf, display, sqrts, paper_N)


if __name__ == "__main__":
    main()
