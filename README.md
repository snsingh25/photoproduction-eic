# photoproduction-eic

## Abstract

Jets are studied in photoproduction ep collisions at the proposed EIC energies,
$\sqrt{s} = 30$–$140$ GeV. The contribution of the photoproduction sub-processes —
direct and resolved — is studied separately. Events are generated with the
**PYTHIA8** event generator, jets are reconstructed with the longitudinally-invariant
$k_T$ / anti-$k_T$ algorithm from **FastJet3**, and the substructure of
gluon- vs quark-initiated jets is analyzed (subjet multiplicities, thick/thin
jet shapes, differential and integrated jet-shape functions).

## Workflow

```
  src/evtgen/        generate events   (PYTHIA8 -> ROOT)
       |
       v
  src/jetreco/       cluster jets      (FastJet)
       |
       +--> src/jetshapes/      integrated & differential jet shapes
       +--> src/subjets/        subjet multiplicity
       +--> src/thinthick/      thick/thin jet classification
       +--> src/dijets/         dijet selections
       +--> src/analysis/       counts, efficiency/purity
       |
       v
  plots/             comparison plots (ROOT / matplotlib / uproot)
```

## Repository layout

```
photoproduction-eic/
├── Makefile                 single top-level build: `make <name>`
├── build.mk                 ROOT/FastJet/PYTHIA paths (edit for your machine)
├── src/
│   ├── evtgen/              event generation
│   ├── jetreco/             jet reconstruction (anti-kT)
│   ├── jetshapes/
│   │   ├── integrated/      psi(r) integrated shape
│   │   └── differential/    rho(r) differential shape
│   ├── subjets/             n-subjet multiplicity
│   ├── thinthick/           thick vs thin jet classification
│   ├── dijets/              dijet-specific analyses
│   └── analysis/            jet counts, efficiency/purity, all-shapes
└── plots/
    ├── jet_kinematics/      ET, eta, 2D ET-vs-eta
    ├── jet_shapes/          integrated & differential shape plotters
    ├── shape_cuts/          cut-study plotters
    ├── subjets/             subjet multiplicity plotters
    ├── range_studies/       eta-range plotters
    ├── efficiency_purity/   purity/efficiency plotters
    └── subprocess_fractions/ QQ/GG/GQ and thin/thick fraction plotters
```

## Building

One Makefile, no per-directory Makefiles. It auto-discovers every `.cc` under
`src/` and `plots/`, so `make <basename>` works from the repo root:

```bash
make evtgen                 # builds bin/evtgen   (links PYTHIA8)
make jetreco                # builds bin/jetreco  (links ROOT + FastJet)
make plot_combined_et_vs_eta
make all                    # builds every .cc in the repo
make list                   # prints every target and its source path
make clean                  # removes bin/
```

Binaries land in `bin/` (git-ignored). Edit the paths in `build.mk` once per
machine (PYTHIA prefix, FastJet include/lib, ROOT via `root-config`).

Python plotters run directly:

```bash
python plots/jet_kinematics/plot_combined_et_vs_eta.py
```

Dependencies for Python plots: `uproot`, `awkward`, `numpy`, `matplotlib`
(some use `ROOT` instead — see the top of each file).

## Data files

ROOT files, logs, and generated PDFs/PNGs are **not tracked** (see
[.gitignore](.gitignore)). Event samples live outside the repo — each script
has its input path near the top; update it to point at your local data.

## Conventions

- `src/<topic>/` holds the analysis binary source (one `.cc` per executable).
- `plots/<topic>/` holds plotters. C++ plotters link ROOT; Python plotters use
  either `uproot`+`matplotlib` or `ROOT`.
- When two scripts do different things but happened to share a filename, they
  were renamed to reflect what they actually plot (e.g. `plot_subprocess_fractions.cc`
  vs `plot_thinthick_fractions.cc`).
