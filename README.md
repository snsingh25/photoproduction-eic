# photoproduction-eic

## Abstract

Jets are studied in photoproduction ep collisions at the EIC design energies,
$\sqrt{s} = 64$–$141$ GeV, alongside a HERA $\sqrt{s} = 300$ GeV reference.
The contribution of the photoproduction sub-processes — direct and resolved —
is studied separately. Events are generated with the **PYTHIA8** event
generator, jets are reconstructed with the longitudinally-invariant
$k_T$ / anti-$k_T$ algorithm from **FastJet3**, and the substructure of
gluon- vs quark-initiated jets is analyzed via integrated jet shapes
$\Psi(r)$, iterated soft-drop multiplicity $n_{\mathrm{SD}}$, and exclusive
$k_T$ subjet multiplicity.

## Papers

Two papers come out of this repository.

**I — Study of quark and gluon jet identification in photoproduction at EIC.**
[arXiv:2603.09545](https://arxiv.org/abs/2603.09545), under revision for PRD.
Establishes the dijet photoproduction samples and uses the integrated jet
shape $\Psi(r)$ together with subjet multiplicity to produce quark- and
gluon-enriched jet samples. The manuscript source lives outside this repo;
the samples, figure code and cross-checks behind it are here, with the
referee-response material under `docs/referee_responses/`.

**II — Comparative study of jet-substructure quark/gluon discriminants in
photoproduction at the EIC.** Source in `EIC_Observables/`, in preparation.
Asks which substructure observable actually separates quark from gluon jets
best at EIC energies, by computing $\Psi(r=0.3)$, $n_{\mathrm{subjets}}$ and
$n_{\mathrm{SD}}$ on *identical* jets — same clustering, same cuts, same
event loop — at all four centre-of-mass energies.

Headline results of II:

- The integrated jet shape gives the largest separation at every energy:
  ROC AUC $0.75$–$0.77$, against $0.66$–$0.69$ for $n_{\mathrm{SD}}$ and
  $0.53$–$0.65$ for $n_{\mathrm{subjets}}$.
- The ordering of the two *counting* observables is not fixed. Pooled over
  $\eta$, $n_{\mathrm{SD}}$ leads by up to $0.16$ in AUC; binned in $\eta$,
  $n_{\mathrm{subjets}}$ leads in all fourteen bins measured. A bin-width
  scan traces this to the absolute momentum threshold in
  $n_{\mathrm{subjets}}$, which leaves it partly sensitive to jet
  kinematics rather than colour flow.

## Workflow

```text
  src/evtgen/        generate events           PYTHIA8 → ROOT
       │                                       data/<sample>/*.root
       ▼
  src/jetreco/       cluster jets + observables FastJet → ROOT
       │                                       data-jets/<sample>/*.root
       │                                       per-jet: psi(r), n_SD, n_subjets,
       │                                                jet_parton_dR, ...
       ▼
  plots/softdrop/    Python plotters           uproot + matplotlib
       │
       ▼
  plots/softdrop/    PDFs + logs               (tracked in repo)
  output/<sample>/
```

## Repository layout

```text
photoproduction-eic/
├── Makefile, build.mk         single top-level build: `make <name>`
├── src/
│   ├── evtgen/                evtgen.cc — PYTHIA8 photoproduction generator
│   ├── jetreco/               jetreco_softdrop.cc — main jet reco + observables
│   │                          jetreco.cc — earlier reco (kept for reference)
│   ├── jetshapes/             integrated/ + differential/ shape analyses
│   ├── subjets/               n-subjet multiplicity
│   ├── thinthick/             thick vs thin jet classification
│   ├── dijets/                older dijet-specific analyses
│   └── analysis/              jet counts, efficiency/purity, all-shapes
├── EIC_Observables/           PAPER II manuscript (REVTeX 4.2)
│   ├── apssamp.tex            main file (template name, real content)
│   ├── EIC_Observables.bib    bibliography
│   ├── figures/               figure PDFs, copied from plots/*/output/
│   └── iith_abstract/         conference abstract
├── plots/
│   ├── style/
│   │   └── paper_style.py     shared publication style — fonts, palettes,
│   │                          REVTeX widths, square panels
│   ├── softdrop/              ACTIVE plotting suite (Python; uproot/matplotlib)
│   │   ├── paper_roc_figures.py           PAPER II Figs 1-4 + AUC tables
│   │   ├── roc_cross_energy_with_psi.py   ROC definition — single source of
│   │   │                                  truth, imported by the above
│   │   ├── roc_cross_energy.py, roc_eta_scan.py, roc_psi_vs_others.py
│   │   ├── sanity_eta_means.py            eta variance decomposition
│   │   ├── paper_fig5_psi_r.py            paper I Fig 5 — Psi(r) per eta bin
│   │   ├── paper_fig8_fig9_regression.py  paper I Fig 8/9 reproducers
│   │   ├── dijet_balance_check.py, jet_parton_match.py
│   │   ├── softdrop_sensitivity.py, algo_compare_kt_vs_antikt.py
│   │   ├── compare_energies.py, compare_nsd_vs_nsubjets.py, check_nsd.py
│   │   └── output/            generated PDFs + summary logs (tracked)
│   │       ├── <sample>/                      per-sample artifacts
│   │       ├── paper_roc_antikt_etmin10/      PAPER II figures + tables
│   │       ├── cross_energy_antikt_etmin10/   all-anti-kT cross-energy scan
│   │       ├── cross_energy_paperconfig/      earlier mixed-algorithm scan
│   │       └── cross_energy_legacy/           pre-paper-config summaries
│   ├── jet_kinematics/
│   │   └── paper_pt_vs_eta.py  jet pT-vs-eta density, one panel per energy
│   ├── jet_shapes/            legacy integrated/differential shape plotters
│   ├── shape_cuts/, range_studies/, subjets/, efficiency_purity/, subprocess_fractions/
│   └─                         legacy plotters (kept for reference)
├── scripts/
│   ├── audit_truth_labels.py            QQ/QG/GG label purity audit (Tier 1)
│   ├── check_table1_table2.py           validate samples vs paper Table I/II
│   ├── run_dijet_balance.sh             driver — sweep dijet_balance_check
│   ├── run_roc_scan.sh                  driver — sweep roc_eta_scan
│   └── output/                          script stdout captures (tracked)
├── docs/
│   ├── methodology_softdrop_and_roc.md  formulas/conventions for the referee
│   ├── tier1_checks.md                  sanity checks log
│   ├── tier2_checks.md                  internal-consistency checks log
│   └── referee_responses/               per-point response drafts
├── data/                      PYTHIA event ROOT (*.root gitignored, logs tracked)
│   └── <sample>/                <sample> = hera300, eic64, eic105, eic141, ...
└── data-jets/                 jetreco_softdrop output ROOT
    └── <sample>/              (*.root gitignored; jetreco run logs tracked)
```

## Building

One Makefile, no per-directory Makefiles. It auto-discovers every `.cc` under
`src/` and `plots/`, so `make <basename>` works from the repo root:

```bash
make evtgen                   # bin/evtgen   (links PYTHIA8)
make jetreco_softdrop         # bin/jetreco_softdrop  (links ROOT + FastJet)
make all                      # builds every .cc in the repo
make list                     # prints every target and its source path
make clean                    # removes bin/
```

Binaries land in `bin/` (git-ignored). Edit the paths in `build.mk` once per
machine (PYTHIA prefix, FastJet include/lib, ROOT via `root-config`).

## Running the analysis

The two-stage pipeline:

```bash
# 1) generate events (preset = HERA_300 | EIC_64 | EIC_105 | EIC_141)
nohup bin/evtgen HERA_300 1000000 paperconfig &

# 2) reconstruct jets (algo = antikt|kt; selection = alljets|dijets)
bin/jetreco_softdrop data/hera300/<file>.root kt dijets 17
```

The Python plotters then read the per-sample ROOT under `data-jets/<sample>/`
and write PDFs/logs to `plots/softdrop/output/<sample>/` automatically:

```bash
python plots/softdrop/roc_eta_scan.py data-jets/hera300_kt_dijets/dijets_*.root
python plots/softdrop/paper_fig5_psi_r.py        # cross-sample
python plots/softdrop/roc_cross_energy_with_psi.py
```

Each script accepts `--out-dir` to override the default output location.

### Regenerating the Paper II figures

Both scripts discover their own inputs and need no arguments:

```bash
python plots/softdrop/paper_roc_figures.py       # Figs 1-4 + AUC tables
python plots/jet_kinematics/paper_pt_vs_eta.py   # jet pT vs eta density
```

They write to `plots/*/output/paper_*/`. Getting a refreshed figure into the
manuscript is a **manual copy** into `EIC_Observables/figures/` — nothing
links the two automatically.

Python deps: `uproot`, `awkward`, `numpy`, `matplotlib`, `scipy` (and a working
TeX install for the paper-figure reproducers — pass `--no-tex` to fall back to
matplotlib's mathtext).

## Data files

Only the bulky ROOT files are gitignored:

- `data/**/*.root`        — PYTHIA event trees
- `data-jets/**/*.root`   — jet trees from `jetreco_softdrop`

Everything else (PDFs, summary logs, per-run `.log` text, run records) is
tracked so plots and provenance stay with the repo. That includes the Paper II
manuscript and its figure PDFs under `EIC_Observables/`; only the LaTeX build
artefacts there are ignored, via a nested `.gitignore`.

## Conventions

- `src/<topic>/` holds the analysis binary source (one `.cc` per executable).
- `plots/<topic>/` holds plotters; `plots/softdrop/` is the active suite. C++
  plotters link ROOT; Python plotters use `uproot`+`matplotlib`.
- Plotter outputs default to `plots/softdrop/output/<sample>/`, derived from
  the input ROOT file's parent directory. Pass `--out-dir` to redirect.
- Per-sample ROOT files in `data-jets/<sample>/` always sit next to a matching
  `dijets_*.log` or `alljets_*.log` recording the `jetreco_softdrop` run.
- Sample directory names encode the cuts (`hera300_antikt_dijets`,
  `hera300_kt_dijets_etmin10`, `hera300_pTHat3`, ...). Variants get their own
  directory rather than overwriting an existing one, so every configuration
  that was ever run stays on disk and reproducible.
- Paper II figures are built on **identical jets** across all four energies:
  anti-$k_T$, $R = 1.0$, exactly two jets, $E_T > 10$ / $7$ GeV. Mixing jet
  algorithms or $E_T$ cuts between samples invalidates the comparison, which
  is the whole point of that paper.
