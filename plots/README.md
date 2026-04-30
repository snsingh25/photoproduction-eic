# `plots/` — plotting suites

Each subdirectory below holds the scripts (and, where applicable, the
generated PDFs) for one analysis topic. Generated artifacts always live in
`<topic>/output/` so the source files stay clean. The active suite for the
PRD revision is **`softdrop/`**; the others are kept for reference and have
not all been re-run against the regenerated `data-jets/` ROOT trees.

## softdrop/ &nbsp;— ACTIVE: paper-figure reproducers + diagnostics

Live-data plotters that read the per-jet observables (`jet_psi_curve_flat`,
`jet_psi03`, `jet_nsd*`, `jet_nsubjets`, `jet_parton_dR`, …) written by
`bin/jetreco_softdrop`. Outputs default to
`softdrop/output/<sample>/` (per-sample) or
`softdrop/output/cross_energy_paperconfig/` (cross-sample summaries).

| File | What it produces |
|---|---|
| `paper_fig5_psi_r.py` | Paper Fig 5 — 4 PDFs of Ψ(r) vs r, one per η bin, all 4 samples overlaid |
| `paper_fig8_fig9_regression.py` | Paper Fig 8/9 reproducers (⟨n_subjets⟩, ⟨Ψ(0.3)⟩ vs η) |
| `roc_eta_scan.py` | Per-η-bin ROC for n_SD vs n_subjets (one sample at a time) |
| `roc_psi_vs_others.py` | 3-way ROC: Ψ(0.3) vs n_SD vs n_subjets on a single sample |
| `roc_cross_energy.py` / `_with_psi.py` | Cross-energy ROC summary (4 samples × 4 η bins) |
| `dijet_balance_check.py` | ⟨pT_2/pT_1⟩, Δφ for the dijet selection (Tier-1 sanity) |
| `jet_parton_match.py` | ΔR(jet, nearest hard parton) (Tier-1 sanity) |
| `compare_nsd_vs_nsubjets.py` | Per-sample shape comparison of n_SD and n_subjets |
| `sanity_eta_means.py` | ⟨observable⟩ vs η for QQ/GG (sanity) |
| `softdrop_sensitivity.py` | n_SD scan over (z_cut, β) variants |
| `algo_compare_kt_vs_antikt.py` | Side-by-side anti-kT vs kT outputs |
| `compare_energies.py` | Cross-sample summary (separation + AUC vs √s) |
| `check_nsd.py` | Quick stdout diagnostic on n_SD branches |

## jet_kinematics/

2-D and 1-D jet kinematic plotters (E_T, η, E_T-vs-η). The Python script
`plot_combined_et_vs_eta.py` is currently the only one wired to the
regenerated repo data (auto-discovers the ROOT files via
`<repo>/data-jets/`); it writes to `jet_kinematics/output/`. The `.cc`
versions and the other Python files still need to be re-pointed at the new
ROOT layout.

## jet_shapes/

Reproductions of the paper Ψ(r) and ρ(r) panels. The Python scripts here
use **frozen paper numbers** (no ROOT input) — they regenerate the original
PDFs verbatim. Use `softdrop/paper_fig5_psi_r.py` if you want Ψ(r) computed
from live data instead. PDFs land in `jet_shapes/output/`.

## subjets/

n-subjets multiplicity plotters. Two scripts here:

- `plotsubjetmultiplicity.py` — single-bin histogram (QQ vs GG); reads a
  `.txt` produced by `bin/nsubjets <root> <etaMin> <etaMax>`. Output PDFs
  land in `subjets/output/` next to the data txts in `subjets/output/data/`.
- `plot_mean_nsubjet_vs_eta.py` — ⟨n_subjets⟩ vs η, single-sample or
  multi-sample overlay. Inputs are the aggregated per-sample txts in
  `subjets/output/data_mean_vs_eta/`.

## efficiency_purity/

Subprocess purity / efficiency analyses (HERA-300 paper context). Mostly
ROOT-C++ plotters; not currently re-run against the regenerated samples.

## subprocess_fractions/

Stacked QQ / GG / GQ fractions and thin-vs-thick jet fractions. Both `.cc`
plotters that need their input paths updated before they run on this repo's
data.

## range_studies/

η-range studies for the integrated shape Ψ(r). Two Python plotters; not
currently re-run against the regenerated samples.

## shape_cuts/

Cut-study plotters for the integrated jet shape (sensitivity to the η or
E_T window). One Python script; not currently re-run.

---

## Conventions

- **Outputs always under `<topic>/output/`** so the script directories stay
  clean. The active suite (`softdrop/`) and `subjets/`, `jet_shapes/`,
  `jet_kinematics/` follow this convention. The legacy suites still write
  to CWD; treat their PDFs as reference material rather than pipeline
  artifacts.
- **Self-contained input discovery**: live-data plotters in `softdrop/`,
  `jet_kinematics/`, and `subjets/` derive their input paths from
  `<repo>/data-jets/` (and `<repo>/data/` for event-level scripts). No
  hard-coded absolute paths.
- **`*.root` is gitignored** under `data/` and `data-jets/`; PDFs and `.log`
  files are tracked, so the plot output committed to the repo is the
  on-disk truth.
