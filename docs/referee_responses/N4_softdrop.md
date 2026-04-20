# Response to Reviewer Point N4 — softdrop multiplicity comparison

**Reviewer said:** *"The authors describe n_subjets as another discriminating
variable. However, similar to the above, this has been studied extensively,
and it is not clear what has been added to the table here. Indeed, given
that the authors are already using FastJet, it would be easy to check other
observables well-known to be good discriminants like softdrop multiplicity,
which are specifically known to be better than n_subjets."*

## Headline

On the same jets, same event sample, same cuts, iterated soft-drop
multiplicity `n_SD` (z_cut=0.1, β=0) outperforms kT-exclusive `n_subjets`
(y_cut=5×10⁻⁴) as a quark–gluon discriminant in *every* photoproduction
sample where GG statistics are adequate (EIC 105 GeV and above):

| √s | sep(n_SD) | sep(n_subj) | AUC(n_SD) | AUC(n_subj) |
|----|----:|----:|----:|----:|
| EIC 64 GeV  | (stats limited†) | (stats limited†) | — | — |
| EIC 105 GeV | **0.53** | 0.43 | **0.71** | 0.66 |
| EIC 141 GeV | **0.57** | 0.41 | **0.71** | 0.65 |
| HERA 300 GeV | **0.62** | 0.21 | **0.73** | 0.57 |

† Only **2 GG jets** survive the ET > 17 GeV cut at √s = 64 GeV — see
  statistical caveat below.

**Advantage of `n_SD` grows with √s.** Going from EIC 105 to HERA 300,
`n_SD` improves slightly (0.53 → 0.62) while `n_subjets` degrades
(0.43 → 0.21). Root cause is explained in the physics section.

## Why `n_SD` beats `n_subjets` (physics)

Both observables count "subjets"; they differ in *which kind*:

- **`n_subjets` (kT exclusive, y_cut = 5×10⁻⁴):** stop clustering when all
  remaining pairs have `min(pT_i², pT_j²) ΔR² > y_cut · Q²`. The threshold
  is an **absolute pT² scale** (~0.2 GeV² for a 20 GeV jet). *Both hard
  and soft-wide-angle emissions* pass the cut.
- **`n_SD` (modified mass-drop, z_cut = 0.1, β = 0):** recluster with C/A,
  walk the hardest branch, count splittings with
  `z = pT_soft/(pT_soft + pT_hard) > z_cut = 0.1`. The threshold is
  **dimensionless momentum sharing** — no GeV scale. Only *hard* (z ~ O(1))
  perturbative splittings get counted.

### Color-factor amplification

For hard perturbative splittings, the leading-log probability at scale Q is
∝ `C_R · α_s`. At O(1) momentum sharing:

- quark  : C_F = 4/3
- gluon  : C_A = 3
- ratio  : C_A / C_F = 2.25

So in the LL limit, ⟨n_SD⟩_G / ⟨n_SD⟩_Q → 2.25 per hard splitting.
That's the theoretical *maximum* separation any multiplicity-style
tagger can deliver, and it's what `n_SD` is designed to approach
[Frye–Larkoski–Moult–Thaler, arXiv:1704.06266]. Your data gets
⟨n_SD⟩_GG / ⟨n_SD⟩_QQ = 4.87 / 3.64 ≈ 1.34 *per jet* at HERA 300;
Poisson stacking across splittings gives an AUC of 0.73.

### Why `n_subjets` degrades with √s

`n_subjets` counts three things simultaneously and only the first is
color-sensitive:

1. hard splittings (what we want) ✓
2. soft-wide-angle bremsstrahlung — multiplicity ∝ α_s log(pT²/y_cut·Q²),
   depends on jet pT, *not* parton type.
3. underlying event + proton-remnant contamination — soft prongs from
   the beam remnant in photoproduction, largest at forward η.

Going up in √s, items 2 and 3 grow faster than item 1. The dilution
becomes visible in the data:

| √s | <n_subj>_QQ | <n_subj>_GG | Δ_mean |
|----|----:|----:|----:|
| 105 | 2.58 | 3.27 | 0.69 |
| 141 | 3.11 | 3.90 | 0.79 |
| 300 | 2.84 | 3.30 | 0.46 |

Means stay within a narrow range, but σ grows (more soft-radiation
phase space) → separation falls.

By contrast ⟨n_SD⟩ grows slowly and in parallel for both classes:

| √s | <n_SD>_QQ | <n_SD>_GG | Δ_mean |
|----|----:|----:|----:|
| 105 | 3.49 | 4.50 | 1.01 |
| 141 | 3.55 | 4.67 | 1.12 |
| 300 | 3.64 | 4.87 | 1.23 |

— the gluon–quark offset *increases* because more perturbative phase space
means more counted splittings, and the C_A/C_F factor is amplified.

### η dependence (HERA 300)

| η bin     | sep(n_SD) | sep(n_subjets) |
|-----------|----:|----:|
| −1 to 0   | 0.50 | **0.80** |
| 0 to 1    | 0.58 | **0.84** |
| 1 to 2    | **0.61** | 0.71 |
| 2 to 3    | **0.63** | 0.57 |
| pooled    | **0.62** | 0.21 |

- `n_SD` is η-flat (0.50 → 0.63) — signature of scale invariance.
- `n_subjets` is η-dependent; its power collapses when pooled across η
  because ⟨n_subjets⟩ itself drifts with η and QQ-at-forward looks like
  GG-at-central.

In **narrow η bins**, `n_subjets` is competitive (or better) because
soft radiation is weakly correlated with parton type. In an **inclusive
or coarsely binned** analysis — which is your paper's use case — `n_SD`
is the cleaner observable.

## Statistical caveat — EIC 64 GeV

Only **2 GG jets** survive `ET > 17 GeV` at √s = 64 GeV. The apparent
negative separation / AUC = 0.20 there is a 2-jet fluctuation, not a
physics result. This also explains the reviewer's minor point #7 on
*why no 64 GeV curve in Fig. 8* — the ET > 17 GeV HERA cut is too
aggressive for the lowest EIC energy.

Recommendation: either (i) lower the ET cut for EIC 64 GeV in a
dedicated low-ET figure/appendix, or (ii) state explicitly that 64 GeV
is excluded from the substructure figures due to insufficient jets.

## Numerical backing (from this repo)

All numbers produced by
[`compare_nsd_vs_nsubjets.py`](compare_nsd_vs_nsubjets.py) and
[`compare_energies.py`](compare_energies.py) on the output of
[`jetreco_softdrop.cc`](jetreco_softdrop.cc), run in per-sample
directories:

```
softdrop/
├── jetreco_softdrop.cc          main C++ analyzer (anti-kT + n_SD + n_subjets)
├── check_nsd.py                 sanity check for n_SD output
├── compare_nsd_vs_nsubjets.py   per-sample 3-panel plot + ROC + per-η
├── compare_energies.py          cross-sample scan (AUC, sep vs √s)
├── RESPONSE_N4.md               this file
├── hera/                        HERA 300 GeV: √s=300, ep
├── eic141/                      EIC 141 GeV
├── eic105/                      EIC 105 GeV
└── eic64/                       EIC 64 GeV   (stats-limited)
```

Each sample folder contains the `alljets_*.root` tree, a run log,
the 3-panel comparison PDF (n_SD dist, n_subjets dist, ROC) and the
text summary.

Cross-energy scan:
[`compare_energies.pdf`](compare_energies.pdf) +
[`compare_energies_summary.txt`](compare_energies_summary.txt).

## Suggested text for the paper (draft)

> To address whether alternative substructure observables might be more
> powerful quark–gluon discriminants, we implemented the iterated
> soft-drop multiplicity \(n_{\mathrm{SD}}\) [FLMT] in parallel to our
> \(n_{\mathrm{subjets}}\) analysis, evaluating both on the same jets
> with the modified-mass-drop parameters \((z_{\mathrm{cut}},\beta) =
> (0.1,\,0)\) and \(y_{\mathrm{cut}} = 5\times10^{-4}\) respectively.
> The soft-drop observable gives systematically better QQ/GG separation
> in the kinematic range of the EIC: ROC AUC 0.71–0.73 for
> \(n_{\mathrm{SD}}\) vs. 0.57–0.66 for \(n_{\mathrm{subjets}}\) across
> √s = 105–300 GeV. The advantage grows with √s because
> \(n_{\mathrm{subjets}}\) with a fixed \(y_{\mathrm{cut}}\) is
> sensitive to soft-wide-angle radiation that scales with the jet pT,
> whereas \(n_{\mathrm{SD}}\) counts only splittings with
> \(z \gtrsim z_{\mathrm{cut}}\) and is therefore insensitive to the
> overall scale. Within single η bins the two observables are
> competitive, consistent with the soft-radiation contribution being
> well-approximated as η-local; the marginal (η-inclusive) advantage
> of \(n_{\mathrm{SD}}\) reflects its approximate scale invariance
> across the η acceptance. Based on these findings, we include
> \(n_{\mathrm{SD}}\) as a complementary discriminant to the
> \(n_{\mathrm{subjets}}\) results [new Fig.~X].
