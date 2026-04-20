#!/usr/bin/env python3
"""Audit the QQ/GG/GQ event-labelling convention (process-code based).

Background: the ROOT files store only final-state hadrons (status 63/84/91
from PYTHIA's hadronisation), so there is no direct per-event
"outgoing-parton PDG" to check. Instead we use the per-event
`processCode` branch, which is the PYTHIA `info.code()` that survived
from the hard 2->2. Each code has a canonical outgoing flavour:

    111: g g -> g g            (GG out)
    112: g g -> q qbar          (QQ out)
    113: q g -> q g             (QG out)
    114: q q' -> q q'           (QQ out)
    115: q qbar -> g g          (GG out)
    116: q qbar -> q' qbar'     (QQ out)
    121: g g -> c cbar          (QQ out)
    122: q qbar -> c cbar       (QQ out)
    123: g g -> b bbar          (QQ out)
    124: q qbar -> b bbar       (QQ out)
    271: g gamma -> q qbar      (QQ out)
    272: g gamma -> c cbar      (QQ out)
    273: g gamma -> b bbar      (QQ out)
    274: q gamma -> q g         (QG out)
    281: gamma g -> q qbar      (QQ out)
    282: gamma g -> c cbar      (QQ out)
    283: gamma g -> b bbar      (QQ out)
    284: gamma q -> q g         (QG out)

The script then answers, for each labelled sample: "what fraction
actually has the outgoing configuration that the label implies?"
  - "QQ_Events" should be ~100% QQ-outgoing
  - "GG_Events" should be ~100% GG-outgoing
  - "GQ_Events" should be ~100% QG-outgoing

If any tree contains codes with a different canonical outgoing config,
that's contamination from whichever classification convention was used
at generation time.

Usage:
  python scripts/audit_truth_labels.py <root_file> [<root_file2> ...]
  (default: all ROOTs under data/allevents_pt7GeV/)
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import uproot


CATEGORIES = ["QQ_Events", "GG_Events", "GQ_Events"]

# code -> (human-readable name, canonical outgoing config)
PROCESS = {
    111: ("g g -> g g",              "GG"),
    112: ("g g -> q qbar",           "QQ"),
    113: ("q g -> q g",              "QG"),
    114: ("q q' -> q q'",            "QQ"),
    115: ("q qbar -> g g",           "GG"),
    116: ("q qbar -> q' qbar'",      "QQ"),
    121: ("g g -> c cbar",           "QQ"),
    122: ("q qbar -> c cbar",        "QQ"),
    123: ("g g -> b bbar",           "QQ"),
    124: ("q qbar -> b bbar",        "QQ"),
    271: ("g gamma -> q qbar (uds)", "QQ"),
    272: ("g gamma -> c cbar",       "QQ"),
    273: ("g gamma -> b bbar",       "QQ"),
    274: ("q gamma -> q g",          "QG"),
    281: ("gamma g -> q qbar (uds)", "QQ"),
    282: ("gamma g -> c cbar",       "QQ"),
    283: ("gamma g -> b bbar",       "QQ"),
    284: ("gamma q -> q g",          "QG"),
}

# What each label claims about the outgoing configuration.
LABEL_EXPECT = {"QQ_Events": "QQ", "GG_Events": "GG", "GQ_Events": "QG"}


def audit_tree_codes(tree):
    """Return Counter(processCode)."""
    codes = tree["processCode"].array(library="np")
    return Counter(codes.tolist())


def pct(n, total):
    return f"{100*n/total:5.2f}%" if total else "  —   "


def audit_sample(sample_name, trees_by_cat):
    print("=" * 92)
    print(f"Sample: {sample_name}")
    print("=" * 92)

    # Aggregate: for each label, counts of outgoing config
    label_outgoing = {cat: Counter() for cat in CATEGORIES}
    label_code_breakdown = {cat: Counter() for cat in CATEGORIES}

    for cat, tree in trees_by_cat.items():
        code_counts = audit_tree_codes(tree)
        label_code_breakdown[cat] = code_counts
        for code, n in code_counts.items():
            out_cfg = PROCESS.get(int(code), (None, "UNKNOWN"))[1]
            label_outgoing[cat][out_cfg] += n

    # --- headline audit table ---
    print(f"\n{'Label':<12}{'N events':>12}"
          f"{'QQ outgoing':>15}{'GG outgoing':>15}{'QG outgoing':>15}{'UNKNOWN':>11}"
          f"  {'match?':>8}")
    print("-" * 92)
    for cat in CATEGORIES:
        if cat not in trees_by_cat:
            continue
        outg = label_outgoing[cat]
        N = sum(outg.values())
        expected = LABEL_EXPECT[cat]
        match_frac = outg.get(expected, 0) / N * 100 if N else 0
        status = "✓" if match_frac >= 95 else ("⚠" if match_frac >= 80 else "✗")
        print(f"{cat:<12}{N:>12}"
              f"{outg['QQ']:>7} {pct(outg['QQ'], N)}"
              f"{outg['GG']:>7} {pct(outg['GG'], N)}"
              f"{outg['QG']:>7} {pct(outg['QG'], N)}"
              f"{outg['UNKNOWN']:>4} {pct(outg['UNKNOWN'], N)}"
              f"  {match_frac:>5.2f}% {status}")

    # --- per-label, per-code breakdown ---
    for cat in CATEGORIES:
        if cat not in trees_by_cat:
            continue
        print(f"\n  {cat}   (expected outgoing: {LABEL_EXPECT[cat]})")
        print(f"  {'code':<6}{'process':<30}{'outgoing':>10}{'N':>10}"
              f"{'% of label':>12}{'OK?':>5}")
        total = sum(label_code_breakdown[cat].values())
        for code in sorted(label_code_breakdown[cat]):
            n = label_code_breakdown[cat][code]
            name, out_cfg = PROCESS.get(int(code), (f"unknown({code})", "UNKNOWN"))
            ok = "✓" if out_cfg == LABEL_EXPECT[cat] else "✗"
            print(f"  {code:<6}{name:<30}{out_cfg:>10}{n:>10}"
                  f"{100*n/total:>11.2f}%  {ok:>4}")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_files", nargs="*")
    args = ap.parse_args()

    if not args.root_files:
        base = Path("data/allevents_pt7GeV")
        args.root_files = sorted(str(p) for p in base.rglob("*.root"))

    for rf in args.root_files:
        with uproot.open(rf) as f:
            trees_by_cat = {}
            for cat in CATEGORIES:
                key = f"{cat}/{cat}"
                if key in f:
                    trees_by_cat[cat] = f[key]
            if not trees_by_cat:
                print(f"[skip] {rf}: no QQ/GG/GQ trees")
                continue
            audit_sample(Path(rf).stem, trees_by_cat)


if __name__ == "__main__":
    main()
