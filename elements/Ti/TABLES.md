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

The same CAS(4,6) Ms=1 CI through three core fields
(`scripts/multiplet_split.py`; the transform validated on carbon and
titanium CAS(2,5) — ladder exact). Carbon's raw field keeps axial
symmetry (README · Motivation); titanium's does not — stronger
open-shell mixing. ⟨L̂z²⟩ after each energy column: it reads M_L² for
the isotropic fields, and shows the breaking for raw (septet trace
Σ = 28 conserved; the raw column pairs its roots by energy order —
M_L is not conserved there).

| M_L | E averaged, Eh | ⟨L̂z²⟩ | E fractional, Eh | ⟨L̂z²⟩ | E raw, Eh | ⟨L̂z²⟩ |
|---|---|---|---|---|---|---|
| −3 | −3.343644 | 9 | −3.462172 | 9 | −3.346418 | 4.0000 |
| −2 | −3.343644 | 4 | −3.462172 | 4 | −3.344106 | 3.9036 |
| −1 | −3.343644 | 1 | −3.462172 | 1 | −3.344106 | 4.1928 |
| 0 | −3.343644 | 0 | −3.462172 | 0 | −3.344106 | 3.9036 |
| +1 | −3.343644 | 1 | −3.462172 | 1 | −3.342274 | 5.9020 |
| +2 | −3.343644 | 4 | −3.462172 | 4 | −3.342274 | 0.1927 |
| +3 | −3.343644 | 9 | −3.462172 | 9 | −3.342274 | 5.9020 |
| ⁵F −3 | −3.296365 | 9 | −3.457635 | 9 | −3.297772 | 5.9023 |
| ⁵F −2 | −3.296365 | 4 | −3.457635 | 4 | −3.297772 | 0.1927 |

## Window lattice — d² series
CASCI, averaged; ROHF reference −848.406 Eh.

| window | CAS | dets | AOs, Eh | ⟨L²⟩ | PySCF MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|---|---|
| 3d | (2,5) | 50 | −848.3959 | 12.000 | — | — | — | — | — | — |
| 3d4p | (2,8) | pending | pending | | | | | | | |
| 3d4p4d | (2,13) | pending | pending | | | | | | | |

## Window lattice — 3d²4s² series

| window | CAS | dets | AOs, Eh | ⟨L²⟩ | PySCF MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|---|---|
| 3d4s | (4,6) | 120 | −848.3950 | 12.000 | −848.4050 | 0.000 | −848.3990 | 12.000 | — | — |
| 3d4s4p | (4,9) | 756 | −848.4290 | 12.000 | — | — | −848.4320 | 12.000 | — | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.4200 | 0.201 | **−848.438** | 12.000 | — | — |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.4390 | 12.000 | — | — | — | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — | — | — | — |
| 3d4d | (4,10) constrained | 1,200 | −848.2320 | 6.000 | — | — | — | — | — | — |
| 3d4s4d | (4,11) constrained | 1,815 | −848.4090 | 12.000 | −848.4140 | 0.010 | — | — | — | — |

## CAS(4,6) / {3d,4s}

C = CASSCF (single root); S = SA-CASSCF (the M_L average). A blank
factor cell inherits from the row above; reading — README · Methods.

**The folds (CASCI — no start, no mask). ROHF anchor: −848.4065 Eh
(−23086.3167 eV):**

| start | mask | substrate | fold | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | fold shift | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CASCI | AOs {3d,4s} | — | raw fold | −848.3959 | +10.5 | −23086.0285 | +0.2857 | 12.0000 | 4.0000 | 4.9e-02 | 2.6e-01 | ~0 | baseline (detector) |
|  |  | — | averaged fold | −848.3959 | +10.5 | −23086.0285 | +0.2857 | 12.0000 | 4.8804* | 4.2e-09 | 2.3e-08 | 2.8 | exact |
|  |  | — | averaged fold, pur | −848.3946 | +11.9 | −23085.9931 | +0.3238 | 12.0000 | 2.5883* | 8.9e-14 | 5.8e-13 | 0.0 | exact (pur) |
|  |  | — | fractional fold | −848.3585 | +48.0 | −23085.0108 | +1.3061 | 12.0000 | 5.5492* | 1.1e-09 | 7.6e-09 | ~0 | different functional |

**The optimizer grid (paired by objective):**

| start | mask | substrate | C E, Eh | C ΔE (mEh) | C E, eV | C ΔE (eV) | C ⟨L²⟩ | C ⟨L̂z²⟩ | C ‖[Ĥ,L̂z]‖ | C ‖[Ĥ,L²]‖ | C fold shift | C verdict | S E, Eh | S ΔE (mEh) | S E, eV | S ΔE (eV) | S ⟨L²⟩ | S ⟨L̂z²⟩ | S ‖[Ĥ,L̂z]‖ | S ‖[Ĥ,L²]‖ | S fold shift | S verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AOs {3d,4s} | free | polarized | −848.4377 | -31.3 | −23087.1660 | -0.8517 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | 2.8 | escape | −848.4166 | -10.1 | −23086.5918 | -0.2748 | 12.0029 | 1.8438 | 1.1e+00 | 6.0e+00 | partition changed | swap |
| PySCF MOs |  |  | −848.4377 | -31.3 | −23087.1660 | -0.8517 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | 2.8 | escape (same minimum) | −848.4316 | -25.1 | −23087.0000 | -0.6830 | 0.0255 | 0.0085 | 3.0e-01 | 7.4e-01 | 2.8 | escape |
| NOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| virtual-mix |  |  | −848.4377 | -31.2 | −23087.1660 | -0.8490 | 0.0496 | 0.0165 | 4.3e-01 | 1.1e+00 | 2.8 | escape (diagnostic) |  |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s} |  | purified |  |  |  |  |  |  |  |  |  |  | −848.4166 | -10.1 | −23086.5918 | -0.2748 | 12.0029 | 3.1567 | 1.0e+00 | 6.0e+00 | partition changed | swap (pur) |
| PySCF MOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| NOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s} | radial | polarized | −848.4112 | -4.7 | −23086.4449 | -0.1279 | 6.9323 | 2.3102 | 3.2e+00 | 1.2e+01 | 2.7 | gauge | −848.4069 | -0.4 | −23086.3278 | -0.0109 | 12.0000 | 4.0000* | 1.3e-03 | 6.9e-03 | 2.8 | floor |
| PySCF MOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| NOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s} |  | purified | −848.4100 | -3.5 | −23086.4122 | -0.0952 | 6.9638 | 2.3214 | 3.2e+00 | 1.2e+01 | 0.0 | gauge | −848.4055 | +1.0 | −23086.2898 | +0.0272 | 12.0000 | 3.9981* | 1.3e-04 | 5.4e-04 | 0.0 | floor |
| PySCF MOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| NOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s} | frozen window | polarized | −848.3961 | +10.4 | −23086.0340 | +0.2830 | 12.0000 | 4.8838* | 4.2e-09 | 2.3e-08 | 2.7 | baseline |  |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| NOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s} |  | purified | −848.3961 | +10.4 | −23086.0340 | +0.2830 | 12.0000 | 8.9049* | 2.8e-14 | 3.7e-13 | 0.0 | exact | −848.3948 | +11.7 | −23085.9986 | +0.3184 | 12.0000 | 1.0687* | 1.0e-13 | 4.0e-13 | 0.0 | exact |
| PySCF MOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| NOs |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |

Source: `rungs.json`.

### CAS(4,10) / {3d,4d} — planned (constrained)

No 4s: the 4s² electrons are forced into 3d/4d. C = CASSCF (single
root); S = SA-CASSCF (the M_L average). A blank factor cell inherits
from the row above; reading — README · Methods.

| start | mask | substrate | C E, Eh | C ΔE (mEh) | C E, eV | C ΔE (eV) | C ⟨L²⟩ | C ⟨L̂z²⟩ | C ‖[Ĥ,L̂z]‖ | C ‖[Ĥ,L²]‖ | C fold shift | C verdict | S E, Eh | S ΔE (mEh) | S E, eV | S ΔE (eV) | S ⟨L²⟩ | S ⟨L̂z²⟩ | S ‖[Ĥ,L̂z]‖ | S ‖[Ĥ,L²]‖ | S fold shift | S verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AOs {3d,4d} | free | polarized | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| virtual-mix |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4d} | radial | polarized | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4d} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |

### CAS(4,14) / {3d,4s,4p,4d} — planned

The physically complete valence window. C = CASSCF (single root);
S = SA-CASSCF (the M_L average). A blank factor cell inherits from the
row above; reading — README · Methods.

| start | mask | substrate | C E, Eh | C ΔE (mEh) | C E, eV | C ΔE (eV) | C ⟨L²⟩ | C ⟨L̂z²⟩ | C ‖[Ĥ,L̂z]‖ | C ‖[Ĥ,L²]‖ | C fold shift | C verdict | S E, Eh | S ΔE (mEh) | S E, eV | S ΔE (eV) | S ⟨L²⟩ | S ⟨L̂z²⟩ | S ‖[Ĥ,L̂z]‖ | S ‖[Ĥ,L²]‖ | S fold shift | S verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AOs {3d,4s,4p,4d} | free | polarized | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| virtual-mix |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s,4p,4d} | radial | polarized | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s,4p,4d} | frozen window | polarized | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| AOs {3d,4s,4p,4d} |  | purified | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| PySCF MOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |
| NOs |  |  | pending |  |  |  |  |  |  |  |  |  | pending |  |  |  |  |  |  |  |  |  |

## Open questions

⟨L²⟩ = 0 MO-window phenomenon (the U ⁵S° class; mechanism open);
semicore (deprioritized); spin–orbit.

## Others (generated between the markers; ticket 28)

<!-- study:begin --> (generated by atomterm.study)

| cell | status | E, Eh | ⟨L²⟩ | ‖[P,L̂z]‖ | ‖[P,L̂⊥]‖ |
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
