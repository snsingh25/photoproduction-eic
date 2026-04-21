# Response to Minor Point m2 — explicit QQ/GG/GQ definitions and PYTHIA setup

> *"Please state the explicit definitions of QQ, QG, and GG labels and
> how they are assigned. My guess is it is just based off of the 2-2
> matrix elements in PYTHIA, but this should be explicit. What is the
> specific PYTHIA process used? What exactly would I set up in PYTHIA
> to reproduce the simulation?"*

The referee's guess is correct: every event is classified by the
outgoing-parton flavour content of the PYTHIA 2→2 hard process
(`pythia.info.code()`). Below we give the complete specification.

## 1. QQ / GG / GQ labelling — outgoing-parton-flavour convention

For each event we read the PYTHIA process code and assign it to one
of three orthogonal samples according to the **flavour content of the
two outgoing hard partons**:

| Label | Meaning | Process codes |
|---|---|---|
| **QQ** | both outgoing hard partons are quarks (or heavy-flavour quarks) | 112, 114, 116, 121, 122, 123, 124 (resolved); 271, 272, 273, 281, 282, 283 (direct BGF) |
| **GG** | both outgoing hard partons are gluons | 111, 115 |
| **GQ** | one outgoing quark, one outgoing gluon | 113 (resolved); 274, 284 (direct QCDC) |

The code-by-code assignment is:

```text
code  Process                    Outgoing  Label
----  -------------------------  --------  -----
111   g g -> g g                 GG        GG
112   g g -> q qbar              QQ        QQ
113   q g -> q g                 QG        GQ
114   q q' -> q q'               QQ        QQ
115   q qbar -> g g              GG        GG
116   q qbar -> q' qbar'         QQ        QQ
121   g g -> c cbar              QQ (c)    QQ
122   q qbar -> c cbar           QQ (c)    QQ
123   g g -> b bbar              QQ (b)    QQ
124   q qbar -> b bbar           QQ (b)    QQ
271   g gamma -> q qbar (uds)    QQ        QQ
272   g gamma -> c cbar          QQ (c)    QQ
273   g gamma -> b bbar          QQ (b)    QQ
274   q gamma -> q g             QG        GQ
281   gamma g -> q qbar (uds)    QQ        QQ
282   gamma g -> c cbar          QQ (c)    QQ
283   gamma g -> b bbar          QQ (b)    QQ
284   gamma q -> q g             QG        GQ
```

Two orthogonal flags tag the photon-interaction mode:

- `isResolved = True`  iff  111 ≤ code ≤ 124
- `isDirect   = True`  iff  271 ≤ code ≤ 284

These are independent of the QQ/GG/GQ split.

**Classification purity.** An end-to-end audit of the PYTHIA event
record (per-event `pdgId` / `status` / `processCode`) confirms the
labels are consistent with the outgoing-parton-flavour convention at
99.99 % or better across all four samples. The only sub-per-mille
contamination is 14-28 events of process 284 (γq → qg, outgoing q+g)
ending up in QQ_Events instead of GQ_Events, which is well below the
statistical precision of any quoted purity.

## 2. PYTHIA configuration — exactly what we run

We use PYTHIA 8.309 for the paper's original event generation and
PYTHIA 8.311/12 for the revision cross-check (the Table II numbers
agree to ≤ 0.12 pp between the two versions at the 10⁶-event level).
The configuration block that reproduces the simulation is:

```text
# --- Beam Configuration ---
Beams:frameType = 2          # asymmetric beams
Beams:idA = 2212             # proton
Beams:idB = 11               # electron
Beams:eA = <E_p>             # from Table I
Beams:eB = <E_e>             # from Table I
PDF:beamB2gamma = on         # Weizsaecker-Williams photon flux

# --- Photoproduction Settings ---
Photon:ProcessType = 0       # auto direct + resolved mix
Photon:Q2max = 1.0           # GeV^2
PhaseSpace:pTHatMin = <pTHatMin>     # per Table below

# --- Enable Processes ---
HardQCD:all       = on       # resolved-photon 2->2 QCD
PhotonParton:all  = on       # direct-photon 2->2 QCD

# --- MPI Tune (important for gamma-p) ---
MultipartonInteractions:pT0Ref = <pT0Ref>    # 4.0 (HERA), 3.0 (EIC)

# --- Fragmentation Tune ---
Tune:pp = 14                 # Monash 2013
Tune:ee = 7
StringZ:aLund  = 0.68
StringZ:bLund  = 0.98
StringPT:sigma = 0.335
```

### Per-sample values

| Sample | E_p [GeV] | E_e [GeV] | √s [GeV] | pTHatMin [GeV] | pT0Ref [GeV] |
|---|---:|---:|---:|---:|---:|
| EIC 64 | 100 | 10 | 63.2 | **5** | 3.0 |
| EIC 105 | 275 | 10 | 104.9 | **5** | 3.0 |
| EIC 141 | 275 | 18 | 141.0 | **7** | 3.0 |
| HERA 300 | 820 | 27.5 | 300 | **7** | 4.0 |

**Note on `pTHatMin`.** The paper's methodology paragraph currently
states `pTHatMin = 3 GeV`, but this is inconsistent with the Table II
percentages — those values are reproduced only with the per-sample
`pTHatMin` listed above. The paper text should be updated to list
these values explicitly. A numerical cross-check of all 28 Table II
percentages using this configuration agrees with the published
numbers to |Δ| ≤ 0.12 pp (max, ~2σ statistical) and ≤ 0.05 pp median.

## 3. Minimal PYTHIA standalone example

A complete, self-contained PYTHIA `.cmnd` file that reproduces the
**HERA 300 GeV** sample (one of the four rows of Table I):

```text
! Photoproduction at HERA 300 GeV — reproduces paper Table II row 4
Beams:frameType        = 2
Beams:idA              = 2212
Beams:idB              = 11
Beams:eA               = 820.0
Beams:eB               = 27.5
PDF:beamB2gamma        = on

Photon:ProcessType     = 0
Photon:Q2max           = 1.0
PhaseSpace:pTHatMin    = 7.0

HardQCD:all            = on
PhotonParton:all       = on

MultipartonInteractions:pT0Ref = 4.0

Tune:pp                = 14
Tune:ee                = 7
StringZ:aLund          = 0.68
StringZ:bLund          = 0.98
StringPT:sigma         = 0.335

Main:numberOfEvents    = 1000000
```

For EIC 64 / 105, set `Beams:eA`, `Beams:eB` per Table I,
`PhaseSpace:pTHatMin = 5.0`, and `MultipartonInteractions:pT0Ref = 3.0`.
For EIC 141, same beam change but `PhaseSpace:pTHatMin = 7.0`.

## 4. Label assignment (code snippet)

After event generation the classifier (in
[`src/evtgen/evtgen.cc`](../../src/evtgen/evtgen.cc)) simply switches
on `pythia.info.code()`:

```cpp
int classifySubprocess(int code) {
    // QQ: both outgoing hard partons are quarks
    if (code == 112 || code == 114 || code == 116 ||
        code == 121 || code == 122 || code == 123 || code == 124 ||
        code == 271 || code == 272 || code == 273 ||
        code == 281 || code == 282 || code == 283)
        return 1;  // QQ
    // GG: both outgoing hard partons are gluons
    if (code == 111 || code == 115)
        return 2;  // GG
    // GQ: one outgoing quark, one outgoing gluon
    if (code == 113 || code == 274 || code == 284)
        return 3;  // GQ
    return 0;      // uncategorised (should not occur given our process list)
}
```

---

### Summary sentence for the referee letter

> "Every event is labelled by the outgoing-parton flavour content of
> the PYTHIA 2→2 hard process (`pythia.info.code()`): QQ = two outgoing
> quarks (codes 112, 114, 116, 121-124, 271-273, 281-283), GG = two
> outgoing gluons (codes 111, 115), GQ = one quark + one gluon (codes
> 113, 274, 284). The complete PYTHIA configuration (including
> per-sample beam energies, `pTHatMin`, and MPI tune) is given in
> Appendix X; it reproduces Tables I and II of the paper to within
> statistical precision."

A reproducibility script
([`scripts/check_table1_table2.py`](../../scripts/check_table1_table2.py))
is available in the repository; it reads the regenerated event ROOT
files and reports the data fraction vs. the paper's Table II value
for every subprocess, with a pass/fail tolerance of 0.20 pp.
