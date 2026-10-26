# Ti — titanium (d²)

cc-pVDZ. Reference **ROHF −848.406489 Eh** (3d²4s², gated good-basin).
M_S = 1. Ground term ³F. Provenance: `ladder.json` (Methods compared — regenerated in-place),
the lattice logs (window lattice), `report.md/json`, `scan.md`. NOTES.md: the gated reference, escape/swap,
JT, ⟨L²⟩ = 0 MO windows.

## Term table — the combinatorial contract

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| CAS(2,5) / {3d} | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 |
| CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 |

## Multiplet energies
CASCI, CAS(2,5), averaged

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³F | 7 | −1.4687 | — | −3…3 |
| ¹D | 5 | −1.4237 | +0.0449 | −2…2 |
| ³P | 3 | −1.4140 | +0.0547 | −1…1 |
| ¹G | 9 | −1.3982 | +0.0704 | −4…4 |
| ¹S | 1 | −1.2951 | +0.1736 | 0 |

## Window lattice
CASCI, different CAS, averaged; ROHF reference −848.406 Eh; the d-atom window already correlates.
Constructions (the start axis): shell = raw AOs; PySCF MO = PySCF's
canonical pick; NO = natural orbitals; custom = our constructed window
(core-orthogonalized / radial products). — = not yet run.

| window | CAS | dets | shell, Eh | ⟨L²⟩ | PySCF MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|---|---|
| 3d | (2,5) | 50 | −848.3959 | 12.000 | — | — | — | — | — | — |
| 3d4s | (4,6) | 120 | −848.3950 | 12.000 | −848.4050 | 0.000 | −848.3990 | 12.000 | — | — |
| 3d4d | (4,10) | 1,200 | −848.2320 | 6.000 | — | — | — | — | — | — |
| 3d4s4p | (4,9) | 756 | −848.4290 | 12.000 | — | — | −848.4320 | 12.000 | — | — |
| 3d4s4d | (4,11) | 1,815 | −848.4090 | 12.000 | −848.4140 | 0.010 | — | — | — | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.4200 | 0.201 | **−848.438** | 12.000 | — | — |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.4390 | 12.000 | — | — | — | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — | — | — | — |

## Methods compared at CAS(4,6)

The full factor space: objective × start × mask × substrate. Pending
rows = planned runs (filled after enable).

| method | start | mask | substrate | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[Ĥ,L̂z]‖ | fold shift | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | −848.4065 | 0 | — | — | 0 | — |
| CASCI | shell {3d,4s} | — | raw fold | −848.3959 | +10.5 | 12.0000 | 4.9e-02 | ~0 | baseline (detector) |
| CASCI | shell {3d,4s} | — | averaged fold | −848.3959 | +10.5 | 12.0000 | 4.2e-09 | 2.8 | exact — the ruler |
| CASCI | shell {3d,4s} | — | averaged fold, pur | −848.3946 | +11.9 | 12.0000 | 8.9e-14 | 0.0 | exact (purified core) |
| CASCI | shell {3d,4s} | — | fractional fold | −848.3585 | +48.0 | 12.0000 | 1.1e-09 | ~0 | different functional |
| CASSCF | shell | free | polarized | −848.4377 | −31.3 | 0.0503 | 4.3e-01 | 2.8 | escape — broken |
| CASSCF | PySCF MOs | free | polarized | −848.4377 | −31.3 | 0.0503 | 4.3e-01 | 2.8 | escape — broken (same basin) |
| CASSCF | virtual-mix | free | polarized | −848.4377 | −31.2 | 0.0496 | 4.3e-01 | 2.8 | escape — broken (diagnostic) |
| CASSCF | customized MOs | free | polarized | pending | | | | | |
| CASSCF | NOs | free | polarized | pending | | | | | |
| CASSCF | shell | free | purified | pending | | | | | |
| CASSCF | PySCF MOs | free | purified | pending | | | | | |
| CASSCF | customized MOs | free | purified | pending | | | | | |
| CASSCF | NOs | free | purified | pending | | | | | |
| CASSCF | shell | radial restriction | polarized | −848.4112 | −4.7 | 6.9323 | 3.2e+00 | 2.7 | gauge — do not pin |
| CASSCF | PySCF MOs | radial restriction | polarized | pending | | | | | |
| CASSCF | customized MOs | radial restriction | polarized | pending | | | | | |
| CASSCF | NOs | radial restriction | polarized | pending | | | | | |
| CASSCF | shell | radial restriction | purified | −848.4100 | −3.5 | 6.9638 | 3.2e+00 | 0.0 | gauge — do not pin |
| CASSCF | PySCF MOs | radial restriction | purified | pending | | | | | |
| CASSCF | customized MOs | radial restriction | purified | pending | | | | | |
| CASSCF | NOs | radial restriction | purified | pending | | | | | |
| CASSCF | shell | frozen window | polarized | −848.3961 | +10.4 | 12.0000 | 4.2e-09 | 2.7 | baseline |
| CASSCF | PySCF MOs | frozen window | polarized | pending | | | | | |
| CASSCF | customized MOs | frozen window | polarized | pending | | | | | |
| CASSCF | NOs | frozen window | polarized | pending | | | | | |
| CASSCF | shell | frozen window | purified | −848.3961 | +10.4 | 12.0000 | 2.8e-14 | 0.0 | machine-exact |
| CASSCF | PySCF MOs | frozen window | purified | pending | | | | | |
| CASSCF | customized MOs | frozen window | purified | pending | | | | | |
| CASSCF | NOs | frozen window | purified | pending | | | | | |
| SA-CASSCF | shell | free | polarized | −848.4166 | −10.1 | 12.0029 | 1.1e+00 | partition changed | swap — broken |
| SA-CASSCF | PySCF MOs | free | polarized | −848.4316 | −25.1 | 0.0255 | 3.0e-01 | 2.8 | escape — broken |
| SA-CASSCF | customized MOs | free | polarized | pending | | | | | |
| SA-CASSCF | NOs | free | polarized | pending | | | | | |
| SA-CASSCF | shell | free | purified | −848.4166 | −10.1 | 12.0029 | 1.0e+00 | partition changed | swap survives purification |
| SA-CASSCF | PySCF MOs | free | purified | pending | | | | | |
| SA-CASSCF | customized MOs | free | purified | pending | | | | | |
| SA-CASSCF | NOs | free | purified | pending | | | | | |
| SA-CASSCF | shell | radial restriction | polarized | −848.4069 | −0.4 | 12.0000 | 1.3e-03 | 2.8 | flat-gauge floor |
| SA-CASSCF | PySCF MOs | radial restriction | polarized | pending | | | | | |
| SA-CASSCF | customized MOs | radial restriction | polarized | pending | | | | | |
| SA-CASSCF | NOs | radial restriction | polarized | pending | | | | | |
| SA-CASSCF | shell | radial restriction | purified | −848.4055 | +1.0 | 12.0000 | 1.3e-04 | 0.0 | flat-gauge floor |
| SA-CASSCF | PySCF MOs | radial restriction | purified | pending | | | | | |
| SA-CASSCF | customized MOs | radial restriction | purified | pending | | | | | |
| SA-CASSCF | NOs | radial restriction | purified | pending | | | | | |
| SA-CASSCF | shell | frozen window | polarized | pending | | | | | |
| SA-CASSCF | PySCF MOs | frozen window | polarized | pending | | | | | |
| SA-CASSCF | customized MOs | frozen window | polarized | pending | | | | | |
| SA-CASSCF | NOs | frozen window | polarized | pending | | | | | |
| SA-CASSCF | shell | frozen window | purified | −848.3948 | +11.7 | 12.0000 | 1.0e-13 | 0.0 | exact — certified |
| SA-CASSCF | PySCF MOs | frozen window | purified | pending | | | | | |
| SA-CASSCF | customized MOs | frozen window | purified | pending | | | | | |
| SA-CASSCF | NOs | frozen window | purified | pending | | | | | |

E = the as-run variational total; ΔE against the ROHF anchor. Fold
shift = |E(own fold) − E(averaged fold)|, max over roots; the two swap
rows read −849.27 in our fold — the partition changed (a core orbital
rotated into the window), so totals only. Starts: shell = raw {3d,4s}
AOs; PySCF MOs = PySCF's own pick; customized MOs = our constructed
(core-orthogonalized, orthonormalized) window; NOs = occupation groups
of the term-averaged CI; virtual-mix = the deliberately m-broken
registry diagnostic. The single-root radial rows are chaotic gauges
(fact-class: never pin them). Source: `ladder.json` (regenerated gated
run, in-place).

### Methods compared at CAS(4,10) / {3d,4d} — planned

No 4s in the window: the 4s² electrons are forced into 3d/4d — a
constrained partition, annotated here before the runs. Rows are filled
by study cells (window × casci / relax / sa-casscf), the pattern of the
Studies insert below.

| treatment | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | verdict |
|---|---|---|---|---|---|
| ROHF (anchor) | −848.4065 | 0 | — | — | — |
| CASCI, averaged | pending | | | | |
| CASSCF, free | pending | | | | |
| CASSCF, core-relax | pending | | | | |
| SA-CASSCF, frozen | pending | | | | |

### Methods compared at CAS(4,14) / {3d,4s,4p,4d} — planned

The physically complete valence window; frozen-CASCI reference values
already in the window lattice above. Same study-cell filling.

| treatment | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | verdict |
|---|---|---|---|---|---|
| ROHF (anchor) | −848.4065 | 0 | — | — | — |
| CASCI, averaged | −848.4380 | −31.5 | 12.0000 | pending | pending |
| CASSCF, free | pending | | | | |
| CASSCF, core-relax | pending | | | | |
| SA-CASSCF, frozen | pending | | | | |

## Frozen-scan verdicts (Motivation; curves in NOTES.md)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| rises (+9.8 mEh @ θ = 0.05); the ~11 mEh gain lives on the orthogonalized s slot | descends (−163 µEh @ δ = 0.1; pair patterns vary) — JT-like |

Ticket-24 amendment: the ceilings are window statements (5d −4.36,
4f −11.46 mEh cross below the escape).

## Open questions

⟨L²⟩ = 0 MO-window phenomenon (the U ⁵S° class; mechanism open);
semicore (deprioritized); spin–orbit.

## Studies (generated between the markers; ticket 28)

<!-- study:begin --> (generated by atomterm.study)

| cell | status | E, Eh | <L2> | |[P,Lz]| | |[P,Lp]| |
|---|---|---|---|---|---|
| study_ti_no_no2_casci | ok | -848.398849 | 10.6709 | 9.4e-15 | 2.6e-14 |
| study_ti_no_no2_relax | ok | -848.400255 | 10.6939 | 9.4e-15 | 2.6e-14 |
| study_ti_no_no2_sa-casscf | ok | -848.416563 | 11.0461 | 9.4e-15 | 2.6e-14 |
| study_ti_no_no3_casci | ok | -848.432430 | 8.5168 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no3_relax | ok | -848.433822 | 8.5669 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no3_sa-casscf | ok | -848.436561 | 9.2758 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no4_casci | ok | -848.438477 | 8.4219 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no4_relax | ok | -848.439848 | 8.4347 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no4_sa-casscf | ok | -848.442374 | 8.0006 | 2.2e-13 | 4.4e-13 |

<!-- study:end -->
