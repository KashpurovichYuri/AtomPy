# Ti — titanium

cc-pVDZ. Reference **ROHF −848.406489 Eh** (3d²4s², gated good-basin).
M_S = 1. Ground term ³F. Provenance: `rungs.json` (Methods compared — regenerated in-place),
the lattice logs (window lattice), `report.md/json`, `scan.md`. Records: the gated reference —
AGENTS fact 21a; escape/swap/JT — the manifest verdicts + README ·
Motivation; ⟨L²⟩ = 0 MO windows — open, elements/U/.

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

### One run, three fields — CAS(4,6) Ms=1 per root (the split demonstration)

The same CAS(4,6) Ms=1 CI read through three core fields, M_L resolved
in the complex (|l,m⟩) basis (`scripts/multiplet_split.py`; transform
validated by the ladder test — C CAS(2,3): 0,1,1; Ti CAS(2,5):
0,0,1,1,1,1,4,4,9,9 — exact). Contrast in the raw fold: carbon's core
anisotropy is tiny ([Ĥ,L̂z] ~ 1e-6) so M_L survives there (its ³P
splits axially, M_L = 0 vs the ±1 pair); titanium's percent-class
field breaks it to D₂.

**Isotropic fields** (averaged ρ̄_c and fractional AOC): M_L is a good
quantum number; term degeneracies exact.

| M_L | ⟨L̂z²⟩ | ³F/⁵F | E averaged, Eh | E fractional, Eh |
|---|---|---|---|---|
| −3 | 9 | ³F | −3.343644 | −3.462172 |
| −2 | 4 | ³F | −3.343644 | −3.462172 |
| −1 | 1 | ³F | −3.343644 | −3.462172 |
| 0 | 0 | ³F | −3.343644 | −3.462172 |
| +1 | 1 | ³F | −3.343644 | −3.462172 |
| +2 | 4 | ³F | −3.343644 | −3.462172 |
| +3 | 9 | ³F | −3.343644 | −3.462172 |
| −3 | 9 | ⁵F | −3.296365 | −3.457635 |
| −2 | 4 | ⁵F | −3.296365 | −3.457635 |

Both isotropic fields give the exact M_L ladder; the term *gaps* are
field-dependent (the ³F–⁵F gap: 47.3 mEh averaged vs 4.5 mEh
fractional) — degeneracies validate, gaps belong to the field.

**Raw (polarized) field**: M_L is not conserved (the real-AO SCF field
breaks axial symmetry to D₂ — this is what the old 1/3/3/2 root-group
pattern actually was); the septet splits and ⟨L̂z²⟩ measures each
root's M_L² content. The septet trace Σ⟨L̂z²⟩ = 28 is conserved — the
split redistributes it.

| raw root | E, Eh | ⟨L̂z²⟩ |
|---|---|---|
| 1 | −3.346418 | 4.0000 |
| 2 | −3.344106 | 3.9036 |
| 3 | −3.344106 | 4.1928 |
| 4 | −3.344106 | 3.9036 |
| 5 | −3.342274 | 5.9020 |
| 6 | −3.342274 | 0.1927 |
| 7 | −3.342274 | 5.9020 |

## Window lattice
CASCI, different CAS, averaged; ROHF reference −848.406 Eh.

| window | CAS | dets | AOs, Eh | ⟨L²⟩ | PySCF MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|---|---|
| 3d | (2,5) | 50 | −848.3959 | 12.000 | — | — | — | — | — | — |
| 3d4s | (4,6) | 120 | −848.3950 | 12.000 | −848.4050 | 0.000 | −848.3990 | 12.000 | — | — |
| 3d4d | (4,10) | 1,200 | −848.2320 | 6.000 | — | — | — | — | — | — |
| 3d4s4p | (4,9) | 756 | −848.4290 | 12.000 | — | — | −848.4320 | 12.000 | — | — |
| 3d4s4d | (4,11) | 1,815 | −848.4090 | 12.000 | −848.4140 | 0.010 | — | — | — | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.4200 | 0.201 | **−848.438** | 12.000 | — | — |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.4390 | 12.000 | — | — | — | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — | — | — | — |

## CAS(4,6) / {3d,4s}

| method | start | mask | substrate | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | fold shift | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | −848.4065 | 0 | −23086.3167 | +0.0000 | — | — | — | — | 0 | — |
| CASCI | AOs {3d,4s} | — | raw fold | −848.3959 | +10.5 | −23086.0285 | +0.2857 | 12.0000 | 4.0000 | 4.9e-02 | 2.6e-01 | ~0 | baseline (detector) |
|  |  | — | averaged fold | −848.3959 | +10.5 | −23086.0285 | +0.2857 | 12.0000 | 4.8804* | 4.2e-09 | 2.3e-08 | 2.8 | exact |
|  |  | — | averaged fold, pur | −848.3946 | +11.9 | −23085.9931 | +0.3238 | 12.0000 | 2.5883* | 8.9e-14 | 5.8e-13 | 0.0 | exact (pur) |
|  |  | — | fractional fold | −848.3585 | +48.0 | −23085.0108 | +1.3061 | 12.0000 | 5.5492* | 1.1e-09 | 7.6e-09 | ~0 | different functional |
| CASSCF | AOs {3d,4s} | free | polarized | −848.4377 | −31.3 | −23087.1660 | −0.8517 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | 2.8 | escape |
|  | PySCF MOs |  |  | −848.4377 | −31.3 | −23087.1660 | −0.8517 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | 2.8 | escape (same minimum) |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | virtual-mix |  |  | −848.4377 | −31.2 | −23087.1660 | −0.8490 | 0.0496 | 0.0165 | 4.3e-01 | 1.1e+00 | 2.8 | escape (diagnostic) |
|  | AOs {3d,4s} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} | radial | polarized | −848.4112 | −4.7 | −23086.4449 | −0.1279 | 6.9323 | 2.3102 | 3.2e+00 | 1.2e+01 | 2.7 | gauge |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} |  | purified | −848.4100 | −3.5 | −23086.4122 | −0.0952 | 6.9638 | 2.3214 | 3.2e+00 | 1.2e+01 | 0.0 | gauge |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} | frozen window | polarized | −848.3961 | +10.4 | −23086.0340 | +0.2830 | 12.0000 | 4.8838* | 4.2e-09 | 2.3e-08 | 2.7 | baseline |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} |  | purified | −848.3961 | +10.4 | −23086.0340 | +0.2830 | 12.0000 | 8.9049* | 2.8e-14 | 3.7e-13 | 0.0 | exact |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
| SA-CASSCF | AOs {3d,4s} | free | polarized | −848.4166 | −10.1 | −23086.5918 | −0.2748 | 12.0029 | 1.8438 | 1.1e+00 | 6.0e+00 | partition changed | swap |
|  | PySCF MOs |  |  | −848.4316 | −25.1 | −23087.0000 | −0.6830 | 0.0255 | 0.0085 | 3.0e-01 | 7.4e-01 | 2.8 | escape |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} |  | purified | −848.4166 | −10.1 | −23086.5918 | −0.2748 | 12.0029 | 3.1567 | 1.0e+00 | 6.0e+00 | partition changed | swap (pur) |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} | radial | polarized | −848.4069 | −0.4 | −23086.3278 | −0.0109 | 12.0000 | 4.0000* | 1.3e-03 | 6.9e-03 | 2.8 | floor |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} |  | purified | −848.4055 | +1.0 | −23086.2898 | +0.0272 | 12.0000 | 3.9981* | 1.3e-04 | 5.4e-04 | 0.0 | floor |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s} |  | purified | −848.3948 | +11.7 | −23085.9986 | +0.3184 | 12.0000 | 1.0687* | 1.0e-13 | 4.0e-13 | 0.0 | exact |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |

The full factor space: objective × start × mask × substrate — the
reading (axes, inheritance, fold shift, ⟨L̂z²⟩ gauge) is README ·
Methods. Source: `rungs.json`.

### CAS(4,10) / {3d,4d} — planned (constrained)

No 4s: the 4s² electrons are forced into 3d/4d. Same factor
format as CAS(4,6); reading — README · Methods.

| method | start | mask | substrate | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | fold shift | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | −848.4065 | 0 | −23086.3167 | +0.0000 | — | — | — | — | 0 | — |
| CASCI | AOs {3d,4d} |  | raw fold | pending |  |  |  |  |  |  |  |  |  |
|  |  |  | averaged fold | pending |  |  |  |  |  |  |  |  |  |
|  |  |  | averaged fold, pur | pending |  |  |  |  |  |  |  |  |  |
|  |  |  | fractional fold | pending |  |  |  |  |  |  |  |  |  |
| CASSCF |  | free | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | virtual-mix |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} | radial | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
| SA-CASSCF | AOs {3d,4d} | free | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} | radial | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |

### CAS(4,14) / {3d,4s,4p,4d} — planned

The physically complete valence window; the averaged-fold values
already in the window lattice above. Reading — README · Methods.

| method | start | mask | substrate | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | fold shift | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | −848.4065 | 0 | −23086.3167 | +0.0000 | — | — | — | — | 0 | — |
| CASCI | AOs {3d,4s,4p,4d} |  | raw fold | pending |  |  |  |  |  |  |  |  |  |
|  |  |  | averaged fold | −848.4380 | -31.5 | −23087.1741 | -0.8572 | 12.0000 |  |  |  |  |  |
|  |  |  | averaged fold, pur | pending |  |  |  |  |  |  |  |  |  |
|  |  |  | fractional fold | pending |  |  |  |  |  |  |  |  |  |
| CASSCF |  | free | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | virtual-mix |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} | radial | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
| SA-CASSCF | AOs {3d,4s,4p,4d} | free | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} | radial | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  |
|  | PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  |
|  | NOs |  |  | pending |  |  |  |  |  |  |  |  |  |

## Frozen-scan verdicts (Motivation; curves in `scan.md`)

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
