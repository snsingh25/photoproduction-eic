#!/usr/bin/env python3
"""Study 1 sanity check for softdrop multiplicity output from jetreco_softdrop.cc.

Flags anything anomalous so we can catch bugs before running on full EIC samples.
"""

import sys
import numpy as np
import uproot

if len(sys.argv) < 2:
    print("Usage: python check_nsd.py <root_file>")
    sys.exit(1)

root_path = sys.argv[1]
categories = ["QQ_Events", "GG_Events", "GQ_Events",
              "Combined_Events", "Resolved_Events", "Direct_Events"]

print(f"\nReading {root_path}\n")

# --- main summary table ------------------------------------------------------
print(f"{'Category':<20}{'N_jets':>10}{'<n_sd>':>10}"
      f"{'median':>10}{'min':>6}{'max':>6}{'std':>8}")
print("-" * 70)

results = {}  # category -> flat numpy array of n_sd
eta_results = {}  # category -> flat numpy array of jet eta

with uproot.open(root_path) as f:
    for cat in categories:
        key = f"{cat}/jets_{cat}"
        if key not in f:
            print(f"{cat:<20}  (not found)")
            continue

        tree = f[key]
        nsd_jagged = tree["jet_nsd"].array(library="np")
        eta_jagged = tree["jet_eta"].array(library="np")

        if len(nsd_jagged) == 0:
            print(f"{cat:<20}  (empty tree)")
            continue

        # Flatten over events
        nsd_flat = np.concatenate(
            [np.asarray(a) for a in nsd_jagged if len(a) > 0]
        ) if len(nsd_jagged) else np.array([])
        eta_flat = np.concatenate(
            [np.asarray(a) for a in eta_jagged if len(a) > 0]
        ) if len(eta_jagged) else np.array([])

        if nsd_flat.size == 0:
            print(f"{cat:<20}  (no jets)")
            continue

        results[cat] = nsd_flat
        eta_results[cat] = eta_flat

        print(f"{cat:<20}{nsd_flat.size:>10d}"
              f"{nsd_flat.mean():>10.2f}"
              f"{np.median(nsd_flat):>10.1f}"
              f"{int(nsd_flat.min()):>6d}"
              f"{int(nsd_flat.max()):>6d}"
              f"{nsd_flat.std(ddof=1):>8.2f}")

# --- automatic sanity checks -------------------------------------------------
print("\n" + "=" * 70)
print("AUTOMATIC SANITY CHECKS")
print("=" * 70)

issues = []

for cat, nsd in results.items():
    if nsd.min() < 0:
        issues.append(f"[FAIL] {cat}: min(n_sd) = {nsd.min()} is negative!")
    if nsd.max() > 15:
        issues.append(f"[WARN] {cat}: max(n_sd) = {nsd.max()} looks high "
                      "(>15 is rare for pT ~ 17 GeV jets)")
    if nsd.mean() < 0.5:
        issues.append(f"[WARN] {cat}: <n_sd> = {nsd.mean():.2f} is very low "
                      "(cut too tight? tree walk broken?)")
    if nsd.mean() > 8:
        issues.append(f"[WARN] {cat}: <n_sd> = {nsd.mean():.2f} is very high "
                      "(cut too loose? counting soft junk?)")

# QQ/GG separation check
if "QQ_Events" in results and "GG_Events" in results:
    sep = results["GG_Events"].mean() - results["QQ_Events"].mean()
    print(f"\n⟨n_sd⟩(GG) − ⟨n_sd⟩(QQ) = {sep:+.3f}")
    if sep <= 0:
        issues.append(f"[FAIL] GG does not have higher <n_sd> than QQ! "
                      f"Separation = {sep:.3f}")
    elif sep < 0.3:
        issues.append(f"[WARN] GG/QQ separation is small ({sep:.3f}); "
                      "may indicate cut or tree-walk issue")
    else:
        print(f"  → Separation is in expected direction and magnitude.")

# GQ should sit between QQ and GG
if all(c in results for c in ["QQ_Events", "GG_Events", "GQ_Events"]):
    qq_mean = results["QQ_Events"].mean()
    gg_mean = results["GG_Events"].mean()
    gq_mean = results["GQ_Events"].mean()
    if not (qq_mean - 0.2 < gq_mean < gg_mean + 0.2):
        issues.append(f"[WARN] GQ mean ({gq_mean:.2f}) not between "
                      f"QQ ({qq_mean:.2f}) and GG ({gg_mean:.2f}) as expected")
    else:
        print(f"  ⟨n_sd⟩: QQ={qq_mean:.2f} < GQ={gq_mean:.2f} < GG={gg_mean:.2f}")
        print(f"  → Ordering is correct (expected from QCD color factors)")

# Resolved vs Direct check
if "Resolved_Events" in results and "Direct_Events" in results:
    res_mean = results["Resolved_Events"].mean()
    dir_mean = results["Direct_Events"].mean()
    print(f"\n⟨n_sd⟩ Resolved = {res_mean:.2f}, Direct = {dir_mean:.2f}")
    if res_mean > dir_mean:
        print("  → Resolved > Direct, consistent with more gluons in resolved processes")

print("\n" + "=" * 70)
if issues:
    print("ISSUES FOUND:")
    for issue in issues:
        print(f"  {issue}")
else:
    print("All automatic checks PASSED.")
print("=" * 70)

# --- histogram of n_sd values (ASCII, quick visual) --------------------------
print("\nn_sd distribution (fraction of jets) — ASCII histogram:")
print(f"{'n_sd':>5}", end="")
for cat in ["QQ_Events", "GG_Events"]:
    if cat in results:
        print(f"{cat:>18}", end="")
print()

for n in range(0, 11):
    print(f"{n:>5}", end="")
    for cat in ["QQ_Events", "GG_Events"]:
        if cat in results:
            frac = (results[cat] == n).sum() / results[cat].size
            bar = "█" * int(frac * 100)
            print(f"  {frac*100:5.1f}% {bar:<10}", end="")
    print()

print("""
Expected:
  QQ peaks around n_sd = 2-3
  GG peaks around n_sd = 3-4 (shifted right of QQ)
  Both distributions should smoothly fall off toward high n_sd
""")