# Methodology: event labels, substructure observables, ROC, and purity/efficiency

This document pins down, in enough detail to reproduce, how events are
classified into QQ/GG/GQ samples, how the two jet-substructure discriminants
($n_{\mathrm{subjets}}$ and $n_{\mathrm{SD}}$) are computed, how the ROC
curves are built from them, and how the ROC relates to the fixed-cut
purity/efficiency analysis used in the paper. It is written so the
relevant sections can be quoted verbatim in the referee response.

Referee comments addressed here:

- **m2** — explicit definitions of QQ / QG / GG labels, PYTHIA configuration.
- **m4** — ensemble-level vs jet-level purity.
- **N4** — soft-drop multiplicity calculation and comparison to $n_{\mathrm{subjets}}$.

---

## 1. Event generation and truth labels

Event generation source: [`src/evtgen/evtgen.cc`](../src/evtgen/evtgen.cc).

### 1.1 PYTHIA 8 configuration

```text
PYTHIA 8.312, frameType = 2  (asymmetric e-p beams)

Beams:idA = 2212                        # proton
Beams:idB = 11                          # electron
Beams:eA  = config.protonEnergy
Beams:eB  = config.electronEnergy
PDF:beamB2gamma     = on                # Weizsaecker-Williams photon flux
Photon:Q2max        = config.Q2max
Photon:ProcessType  = 0                 # automatic direct + resolved mix
PhaseSpace:pTHatMin = config.pTHatMin   # 7 GeV for the pT_hat > 7 samples

HardQCD:all         = on                # resolved-photon 2->2 QCD
PhotonParton:all    = on                # direct-photon 2->2 QCD

MultipartonInteractions:pT0Ref = config.pT0Ref     # MPI tune
StringZ:aLund                  = config.StringZ_aLund
StringZ:bLund                  = config.StringZ_bLund
StringPT:sigma                 = config.StringPT_sigma
```

`Photon:ProcessType = 0` lets PYTHIA mix direct and resolved
photoproduction automatically. `HardQCD:all` enables the Born 2->2 QCD
matrix elements for the resolved case (the photon resolves into a
hadronic state whose partons then scatter off the proton partons);
`PhotonParton:all` enables the direct case (pointlike photon scattering
on a single proton parton).

### 1.2 Subprocess classification (QQ / GG / GQ)

Each event is labelled by its PYTHIA process code `pythia.info.code()`
using the function `classifySubprocess` in `evtgen.cc`. The convention
is **outgoing hard-process parton flavour**:

| Label | Outgoing partons | Codes | Processes |
|---|---|---|---|
| **QQ** | quark + quark | 112, 114, 116, 121, 122, 123, 124 (resolved); 271, 272, 273, 281, 282, 283 (direct BGF) | $gg\to q\bar q$; $qq'\to qq'$; $q\bar q\to q'\bar q'$; $gg\to c\bar c$; $q\bar q\to c\bar c$; $gg\to b\bar b$; $q\bar q\to b\bar b$; $g\gamma\to q\bar q$ (uds,c,b); $\gamma g\to q\bar q$ (uds,c,b) |
| **GG** | gluon + gluon | 111, 115 | $gg\to gg$; $q\bar q\to gg$ |
| **GQ** | quark + gluon | 113 (resolved); 274, 284 (direct QCDC) | $qg\to qg$; $q\gamma\to qg$; $\gamma q\to qg$ |

Two orthogonal flags tag the photon interaction mode:

- `isResolved = True`  iff  $111 \le$ code $\le 124$
- `isDirect   = True`  iff  $271 \le$ code $\le 284$

**Audit.** Running
[`scripts/audit_truth_labels.py`](../scripts/audit_truth_labels.py)
on every event in the four samples confirms the labelling is
**consistent with the outgoing-flavour convention at 99.99% or
better**: the only contamination is 14–28 events per sample (<0.007%)
of process 284 ($\gamma q \to qg$) in the QQ_Events tree instead of
the GQ_Events tree, which is small enough to be neglected. Per-sample
match rates (the remaining fraction being this 284 leakage):

| Sample | QQ_Events | GG_Events | GQ_Events |
|---|---:|---:|---:|
| HERA 300 GeV | 100.00% | 100.00% | 100.00% |
| EIC 141 GeV  | 99.99%  | 100.00% | 100.00% |
| EIC 105 GeV  | 99.99%  | 100.00% | 100.00% |
| EIC 64 GeV   | 99.99%  | 100.00% | 100.00% |

The physics language in the paper ("QQ ≈ quark–quark dijets", etc.)
is therefore fully justified.

---

## 2. Jet reconstruction

All analyses use the same jet definition (see
[`src/jetreco/jetreco_softdrop.cc`](../src/jetreco/jetreco_softdrop.cc)):

- Algorithm: **anti-$k_T$** with $R = 1.0$ (FastJet).
- Particle-level cuts: $p_T > 0.1$ GeV, $|\eta| < 5$.
- Jet-level cuts: $E_T > 17$ GeV, $-4 < \eta < 4$.

Both substructure observables ($n_{\mathrm{subjets}}$ and $n_{\mathrm{SD}}$)
are computed on the same jets inside the same event loop, so every jet
yields a matched pair (`jet_nsubjets`, `jet_nsd`) written to the output
tree ([`jetreco_softdrop.cc:216-266`](../src/jetreco/jetreco_softdrop.cc#L216-L266)).

---

## 3. $n_{\mathrm{subjets}}$ — kT exclusive multiplicity

For each anti-$k_T$ jet with constituents $\{p_i\}$:

1. Recluster with the **$k_T$** algorithm at the same $R = 1$.
2. Use the $k_T$ pair distance
   $d_{ij} = \min(p_{T,i}^2, p_{T,j}^2)\,\Delta R_{ij}^2 / R^2$.
3. Stop clustering at the point where all remaining $d_{ij} > y_{\mathrm{cut}}\,Q^2$
   (where $Q$ is the jet scale). In FastJet, this is
   `cs_kt.exclusive_jets_ycut(y_cut).size()`.
4. $n_{\mathrm{subjets}}$ is the number of proto-jets at the stopping point.

Parameter: $y_{\mathrm{cut}} = 5\times 10^{-4}$, chosen to match the
ZEUS/H1 convention and ensure the result can be compared to archival
HERA measurements.

The $k_T$ pairwise threshold is an **absolute $p_T^2$ scale**
($\sim 0.2~\mathrm{GeV}^2$ for a 20 GeV jet at $y_{\mathrm{cut}}=5\times 10^{-4}$).
Both hard perturbative splittings and soft-wide-angle bremsstrahlung
survive it, so $n_{\mathrm{subjets}}$ is sensitive to both color-structure
information and overall jet kinematics (pT, underlying event, remnants).

---

## 4. $n_{\mathrm{SD}}$ — iterated soft-drop multiplicity

**Reference**: Frye, Larkoski, Moult, Thaler, JHEP 07 (2017) 064
[arXiv:1704.06266] (hereafter **FLMT**).

### 4.1 Formula

$$
n_{\mathrm{SD}}(z_{\mathrm{cut}},\beta)
\;=\;
\sum_{\text{primary-branch splittings } i}
\Theta\!\left[\,z_i \;>\; z_{\mathrm{cut}}\!\left(\frac{\Delta R_i}{R_0}\right)^{\!\beta}\,\right]
$$

where for each 1$\to$2 splitting into prongs $(p_1,p_2)$ with
$p_{T,1}\ge p_{T,2}$:

$$
z_i \equiv \frac{p_{T,2}}{p_{T,1}+p_{T,2}} \in (0, \tfrac{1}{2}],
\qquad
\Delta R_i \equiv \sqrt{(\Delta\eta)^2+(\Delta\phi)^2}.
$$

$\Theta$ is the Heaviside step. The sum runs over splittings along the
**primary branch**: at every split we follow the harder prong $p_1$
inward and never enter $p_2$'s history. Termination is when the harder
prong is irreducible (a single constituent).

### 4.2 Parameter choice in this paper

| symbol | value | role |
|---|---|---|
| $z_{\mathrm{cut}}$ | 0.1 | minimum momentum fraction of the softer prong |
| $\beta$ | 0 | angular exponent; $\beta=0$ gives an angle-independent cut |
| $R_0$ | $R = 1$ | angular normalisation; drops out for $\beta=0$ |

For $\beta=0$ the Heaviside condition reduces to $z > z_{\mathrm{cut}} = 0.1$
— the softer prong must carry at least 10% of the pair's $p_T$. This
is exactly the **modified Mass-Drop Tagger (mMDT)** condition of
Dasgupta–Fregoso–Marzani–Salam [arXiv:1307.0007], applied iteratively
(we keep counting past the first pass, unlike the Soft-Drop groomer
which stops at the first pass).

### 4.3 Algorithm

For each anti-$k_T$ jet $J$ ([`jetreco_softdrop.cc:224-254`](../src/jetreco/jetreco_softdrop.cc#L224-L254)):

1. Extract the constituents $\{p_k\}$ of $J$. If fewer than 2, return
   $n_{\mathrm{SD}}=0$.
2. Recluster with **Cambridge/Aachen** at radius $R_0 = R$. C/A's
   clustering metric $d_{ij} = \Delta R_{ij}^2$ is angular-ordered, which
   matches the QCD coherence/angular-ordering property of parton
   showers. C/A is the canonical choice for declustering;
   $k_T$ would interleave emissions at different angles, and anti-$k_T$
   yields no meaningful substructure tree.
3. Take the leading C/A jet (highest $p_T$) as the root of the tree.
4. Initialise `current` $\leftarrow$ root, `n_SD` $\leftarrow 0$.
5. While `current.has_parents(p1, p2)`:
   1. If $p_{T,1} < p_{T,2}$ swap so $p_1$ is the harder prong.
   2. Compute $\Delta R = p_1.\text{delta\_R}(p_2)$ and
      $z = p_{T,2}/(p_{T,1}+p_{T,2})$.
   3. If $z > z_{\mathrm{cut}}\,(\Delta R/R_0)^\beta$, increment `n_SD`.
   4. Set `current` $\leftarrow p_1$ (follow the harder branch
      regardless of whether the split passed).
6. Return `n_SD`.

Step 5.3 is **counting**, not grooming — we never discard $p_2$. Step
5.4 is what makes this "iterated": we keep walking past a passing
split, counting every qualifying emission off the primary branch.

### 4.4 Why C/A, and why the primary branch

- **Cambridge/Aachen**: pure angular clustering produces an
  angular-ordered declustering tree. Each declustering step
  corresponds to one emission in angle — exactly the ordering generated
  by coherent QCD showers. Using $k_T$ instead would give a
  $p_T$-ordered history and mix up emissions at different angles.
- **Primary (hardest) branch**: at every splitting the harder subjet
  is the continuation of the original hard parton; the softer subjet
  is the emission. Walking the hard branch traces the initial parton
  through its emission history — the "primary Lund plane" trajectory.

### 4.5 Color-factor argument — why it tags q vs g

At leading log, the primary-branch splitting density in the Lund plane is

$$
\frac{dP}{d\log(1/z)\,d\log(R_0/\Delta R)}
\;\approx\;
\frac{\alpha_s(p_T \Delta R)}{\pi}\, C_R\, \bar P_i(z)
$$

with $\bar P_i(z)$ the integrated DGLAP splitting kernel and $C_R$ the
Casimir of the initiating parton:

- Quark: $C_R = C_F = 4/3$
- Gluon: $C_R = C_A = 3$

Integrating over the allowed Lund region (splittings with
$z > z_{\mathrm{cut}}(\Delta R/R_0)^\beta$) gives
$\langle n_{\mathrm{SD}}\rangle \propto C_R$, so the theoretical ratio

$$
\frac{\langle n_{\mathrm{SD}}\rangle_{\mathrm{gluon}}}{\langle n_{\mathrm{SD}}\rangle_{\mathrm{quark}}}
\;\longrightarrow\; \frac{C_A}{C_F} = 2.25
$$

in the LL limit. This is the theoretical maximum any
counting-style tagger can reach from perturbative QCD. In the paper's
HERA-300 sample we measure $\langle n_{\mathrm{SD}}\rangle_{GG}/\langle n_{\mathrm{SD}}\rangle_{QQ} = 4.87/3.64\approx 1.34$ per jet; Poissonian stacking
across multiple splittings delivers the observed ROC AUC of 0.73.

### 4.6 Infrared and collinear safety

- **IR-safe for all $\beta$**: $z\to 0$ emissions fail the cut and
  are not counted.
- **Collinear-safe for $\beta > 0$**: $(\Delta R/R_0)^\beta\to 0$ tightens
  the cut in the collinear limit, removing the collinear divergence.
- **$\beta = 0$** (our choice) is collinear-*unsafe* at fixed order —
  the rectangular Lund region extends to arbitrarily small angles —
  but it is **Sudakov-safe**: the divergence exponentiates into the
  Sudakov form factor and gives a finite all-orders answer. FLMT
  compute $n_{\mathrm{SD}}$ at modified-LL / NLL accuracy with modified
  factorisation. In practice, a natural IR cutoff (hadronisation scale
  / track resolution / constituent mass) renders the count finite per
  event and well-defined for Monte-Carlo comparison.

### 4.7 Lund-plane picture

Each primary-branch splitting lives at a point
$(\log(1/\Delta R),\,\log(1/z))$ on the primary Lund plane. $n_{\mathrm{SD}}$
counts points below the line

$$
\log\!\frac{1}{z} \;<\; \log\!\frac{1}{z_{\mathrm{cut}}} + \beta\,\log\!\frac{R_0}{\Delta R}.
$$

For $\beta=0$ this is a horizontal cut at $\log(1/z) < \log 10 \approx 2.3$:
all splittings with $z > 0.1$, at any angle, are counted.

---

## 5. ROC curve computation

Implementation: [`plots/softdrop/roc_eta_scan.py:roc_integer`](../plots/softdrop/roc_eta_scan.py),
[`plots/softdrop/roc_cross_energy.py`](../plots/softdrop/roc_cross_energy.py).

### 5.1 Definitions

Let $X$ denote a discriminant (e.g. $n_{\mathrm{SD}}$ or $n_{\mathrm{subjets}}$).
With the convention that gluon jets have *larger* $X$:

- **Signal** $=$ jets from GG events
- **Background** $=$ jets from QQ events

For every threshold $c$:

$$
\varepsilon_S(c) \;=\; \Pr(X \ge c \mid \mathrm{GG}),
\qquad
\varepsilon_B(c) \;=\; \Pr(X \ge c \mid \mathrm{QQ}).
$$

The ROC curve is the locus $\bigl(\varepsilon_B(c),\varepsilon_S(c)\bigr)$ as
$c$ sweeps its support.

Both $n_{\mathrm{SD}}$ and $n_{\mathrm{subjets}}$ are integer-valued, so $c$ is
swept over integers from $X_{\min}$ to $X_{\max}+1$. This produces a
piecewise-linear ROC with one node per achievable cut; interpolation
across nodes is not implied (and not needed for AUC).

### 5.2 AUC

$$
\mathrm{AUC} \;=\; \int_0^1 \varepsilon_S(\varepsilon_B)\,d\varepsilon_B
$$

computed by trapezoidal integration on the $(\varepsilon_B,\varepsilon_S)$
nodes sorted by ascending $\varepsilon_B$. AUC = 0.5 is random-guess; AUC = 1
is perfect separation; AUC > 0.5 means the signal class has systematically
larger $X$.

### 5.3 Mean-separation metric

A one-number summary used alongside AUC:

$$
\mathcal{S} \;=\; \frac{\langle X\rangle_{\mathrm{GG}} - \langle X\rangle_{\mathrm{QQ}}}
                      {\sqrt{\sigma^2_{\mathrm{QQ}} + \sigma^2_{\mathrm{GG}}}}
$$

Positive when GG has the higher mean. $\mathcal{S}$ is sensitive to mean
shift only; AUC additionally captures shape differences.

### 5.4 $\eta$-binning and Simpson's-paradox caveat

For per-$\eta$-bin ROCs we slice both populations on $\eta \in [\eta_{\min},
\eta_{\max})$ *before* computing the ROC. The resulting per-bin AUCs can
be substantially larger than the pooled (inclusive) AUC. This is not a
bug: the pooled variance decomposes as

$$
\sigma^2_{\mathrm{pool}}
\;=\;
\underbrace{\langle\sigma^2_{\mathrm{bin}}\rangle_{\mathrm{bins}}}_{\text{within}}
\;+\;
\underbrace{\mathrm{Var}_{\mathrm{bins}}(\langle X\rangle_{\mathrm{bin}})}_{\text{between}}
$$

so any observable whose mean varies with $\eta$ carries extra pooled
width (the "between" term), widening the pooled QQ and GG distributions
and hurting inclusive AUC. For the HERA-300 sample, the between-bin
component is ~46% of $\sigma^2_{\mathrm{pool}}$ for $n_{\mathrm{subjets}}$
(QQ) and ~62% for $n_{\mathrm{subjets}}$ (GG); for $n_{\mathrm{SD}}$ it is
< 1%. This quantitatively explains the inclusive-vs-per-bin AUC
inversion. Diagnostic script:
[`plots/softdrop/sanity_eta_means.py`](../plots/softdrop/sanity_eta_means.py).

---

## 6. Relationship to purity and efficiency

### 6.1 What the paper does

The paper's
[`src/analysis/efficiency_purity_analysis.cc`](../src/analysis/efficiency_purity_analysis.cc)
picks a **fixed classification cut** on the integrated jet shape $\Psi(r=0.3)$:

- *thin* jet $\equiv \Psi(0.3) > 0.8$  (quark-like)
- *thick* jet $\equiv \Psi(0.3) < 0.6$  (gluon-like)

and measures:

- **Efficiency** $\varepsilon_X = \Pr(\text{classified as expected} \mid \text{truth label } X)$
  (per-jet or per-event).
- **Purity** (sample composition) $\Pr(\text{truth label } X \mid \text{classified as } X)$,
  quoted as "fraction of thin jets from QQ events" etc.

### 6.2 ROC as the continuous generalisation

A ROC curve is the locus of $(\varepsilon_B,\varepsilon_S)$ as the classification
cut is swept. The paper's fixed-cut figure is therefore one point on a
ROC we've now traced fully. Given any operating point
$(\varepsilon_B,\varepsilon_S)$ and pre-selection populations $N_{\mathrm{QQ}}$,
$N_{\mathrm{GG}}$, $N_{\mathrm{GQ}}$, the purity of the passing sample is

$$
\mathrm{Purity}_{\mathrm{GG}}
\;=\;
\frac{N_{\mathrm{GG}}\,\varepsilon_S}
     {N_{\mathrm{QQ}}\,\varepsilon_B + N_{\mathrm{GG}}\,\varepsilon_S + N_{\mathrm{GQ}}\,\varepsilon_{\mathrm{GQ}}}.
$$

So reading purity off the ROC is a one-line calculation. The ROC is
strictly more informative: it answers "at any required gluon purity,
what's the achievable efficiency and how does it depend on $\eta$ and
$\sqrt{s}$?" — whereas the fixed-cut analysis fixes the trade-off.

### 6.3 Ensemble-level vs jet-level (referee m4)

Both the paper's purity numbers and the AUCs quoted here are
**ensemble-level**: they refer to fractions within populations of
classified jets (e.g. "70% of passing jets come from QQ events"), not
to per-jet probabilities. The ROC's $\varepsilon_S$ and $\varepsilon_B$
are expectation values over the QQ and GG populations; a jet-level
interpretation of the discriminant would require a per-jet probability
map (e.g. from a calibrated classifier output), which is not what is
done in this paper.

---

## 7. Quick cross-reference for the referee response

- **m2**: §1 gives the explicit PYTHIA configuration and the exact
  process-code lists for QQ/GG/GQ, plus the proton-side-initial-parton
  convention caveat.
- **m4**: §6.3 states the purity is ensemble-level.
- **N4**: §4 provides the formula and algorithm for $n_{\mathrm{SD}}$;
  §5 explains the ROC construction; detailed headline numbers,
  per-$\eta$ breakdown, and cross-energy summary live in
  [`docs/referee_responses/N4_softdrop.md`](referee_responses/N4_softdrop.md).
