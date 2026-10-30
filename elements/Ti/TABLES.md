# Ti — titanium

cc-pVDZ. Reference **ROHF −848.406489 Eh** (3d²4s², gated good-basin).
M_S = 1. Ground term ³F. Provenance: `cells.json` (Methods compared — key-addressed, regenerated in-place),
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
(the retired `scripts/multiplet_split.py` — git history; the transform
validated on carbon and
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

## Window lattice — 3d² series

CASCI; E = the run's own total (as-run). ROHF reference −848.4065 Eh.
Reading — README · Results (window lattice). Source: `cells.json` (`*_2-0_*`).

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|
| 3d | CAS(2,5) | 10 | −848.4001 | 12.000 | −848.4065 | 0.000 | — | — |
| 3d4p | CAS(2,8) | 28 | −848.4001 | 12.000 | −848.4065 | 3.426 | — | — |
| 3d4p4d | CAS(2,13) | 78 | −848.4096 | 12.000 | −848.4065 | 5.150 | — | — |

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

Reading — README · Methods (⟨L̂z²⟩ is a gauge diagnostic, not a label). ROHF anchor: −848.4065 Eh (−23086.3167 eV).

| method | start | mask | substrate | fold | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | Span (Lowdin) | Elec (Lowdin) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **CASCI** | AO {3d4s} | — | polarized | raw | −848.3959 | +10.5 | −23086.0413 | +0.2871 | 12.0000 | 4.0000 | 4.9e-02 | 2.6e-01 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | baseline (detector) |
|  |  | — | polarized | averaged | −848.3959 | +10.5 | −23086.0413 | +0.2871 | 12.0000 | 4.8832 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | exact |
|  |  | — | purified | averaged | −848.3946 | +11.9 | −23086.0044 | +0.3239 | 12.0000 | 4.0214 | 3.4e-14 | 4.6e-13 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | exact (pur) |
|  |  | — | polarized | fractional | −848.3585 | +48.0 | −23085.0230 | +1.3053 | 12.0000 | 5.5419 | 1.0e-09 | 7.6e-09 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | different functional |
| **CASSCF** | AO {3d4s} | free | polarized | averaged | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 | s1.8 p0.2 d2.0 f0 | escape |
|  | PySCF MO |  |  |  | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 | s1.8 p0.2 d2.0 f0 | escape (same minimum) |
|  | NO |  |  |  | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 | s1.8 p0.2 d2.0 f0 | escape (same minimum) |
|  | AO {3d4s} | virtual-mix |  |  | −848.4377 | −31.2 | −23087.1779 | −0.8496 | 0.0496 | 0.0165 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 | s1.8 p0.2 d2.0 f0 | escape (diagnostic) |
|  | AO {3d4s} | free | purified | averaged | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 | s1.8 p0.2 d2.0 f0 | escape |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 0.0503 | 0.0168 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 | s1.8 p0.2 d2.0 f0 | escape (same minimum) |
|  | AO {3d4s} | radial | polarized | averaged | −848.4112 | −4.7 | −23086.4576 | −0.1292 | 6.9323 | 2.3102 | 3.2e+00 | 1.3e+01 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | gauge |
|  | AO {3d4s} | radial | purified | averaged | −848.4100 | −3.5 | −23086.4234 | −0.0951 | 6.9638 | 2.3214 | 3.2e+00 | 1.2e+01 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | gauge |
|  | AO {3d4s} | frozen window | polarized | averaged | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 12.0000 | 4.8846 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | baseline |
|  | PySCF MO |  |  |  | −848.4171 | −10.6 | −23086.6177 | −0.2893 | 0.0051 | 0.0017 | 1.6e-01 | 3.9e-01 | s0.17 p0.50 d0.33 f0 | s1.9 p0.1 d2.0 f0 | MO start was not closed |
|  | NO |  |  |  | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 12.0000 | 4.8865 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | baseline |
|  | AO {3d4s} | frozen window | purified | averaged | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 12.0000 | 7.0253 | 3.4e-14 | 3.4e-13 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | exact (pur) |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 12.0000 | 3.8854 | 5.9e-14 | 1.0e-12 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | exact (pur) |
| **SA-CASSCF** | AO {3d4s} | free | polarized | averaged | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 12.0029 | 1.8438 | 1.1e+00 | 6.0e+00 | s0 p0.17 d0.83 f0 | s0 p2.0 d2.0 f0 | swap |
|  | PySCF MO |  |  |  | −848.4316 | −25.1 | −23087.0118 | −0.6835 | 0.0255 | 0.0085 | 3.0e-01 | 7.4e-01 | s0.17 p0.50 d0.33 f0 | s1.8 p0.2 d2.0 f0 | escape |
|  | NO |  |  |  | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 12.0029 | 1.5114 | 1.2e+00 | 6.0e+00 | s0 p0.17 d0.83 f0 | s0 p2.0 d2.0 f0 | swap |
|  | AO {3d4s} | free | purified | averaged | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 12.0029 | 3.1567 | 1.0e+00 | 6.0e+00 | s0 p0.17 d0.83 f0 | s0 p2.0 d2.0 f0 | swap (pur) |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | −848.4165 | −10.0 | −23086.6002 | −0.2719 | 12.0034 | 4.5601 | 9.1e-01 | 6.1e+00 | s4e-04 p0.17 d0.83 f0 | s0.0 p2.0 d2.0 f0 | swap (pur) |
|  | AO {3d4s} | radial | polarized | averaged | −848.4069 | −0.4 | −23086.3399 | −0.0115 | 12.0000 | 4.0000 | 1.3e-03 | 6.9e-03 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | floor |
|  | AO {3d4s} | radial | purified | averaged | −848.4055 | +1.0 | −23086.3017 | +0.0266 | 12.0000 | 6.0951 | 1.1e-04 | 5.5e-04 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | floor |
|  | AO {3d4s} | frozen window | polarized | averaged | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 12.0000 | 4.8835 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | baseline |
|  | PySCF MO |  |  |  | −848.4170 | −10.5 | −23086.6151 | −0.2868 | 0.0051 | 0.0017 | 1.6e-01 | 3.9e-01 | s0.17 p0.50 d0.33 f0 | s1.9 p0.1 d2.0 f0 | MO start was not closed |
|  | NO |  |  |  | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 12.0000 | 4.8998 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | baseline |
|  | AO {3d4s} | frozen window | purified | averaged | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 12.0000 | 0.9855 | 3.4e-14 | 6.1e-13 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | exact |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 12.0000 | 2.3486 | 4.3e-14 | 4.4e-13 | s0.17 p0 d0.83 f0 | s2.0 p0 d2.0 f0 | exact |

Source: `cells.json` (this draw; verdicts — the AUDIT table of `scripts/build_methods_tables.py`).


### CAS(4,10) / {3d,4d} — planned (constrained)

No 4s (constrained partition). Pending = declared, not yet run; the CASCI anchor lives in the 3d²4s² lattice above.

| method | start | mask | substrate | fold | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | Span (Lowdin) | Elec (Lowdin) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **CASCI** | AO {3d4d} | — | polarized | raw | pending |  |  |  |  |  |  |  |  |  | pending |
|  |  | — | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  |  | — | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  |  | — | polarized | fractional | pending |  |  |  |  |  |  |  |  |  | pending |
| **CASSCF** | AO {3d4d} | free | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | virtual-mix |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | free | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | radial | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | radial | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | frozen window | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | frozen window | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
| **SA-CASSCF** | AO {3d4d} | free | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | free | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | radial | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | radial | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | frozen window | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4d} | frozen window | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |


### CAS(4,14) / {3d,4s,4p,4d} — planned

The complete valence window. Pending = declared, not yet run; the CASCI anchor lives in the 3d²4s² lattice above.

| method | start | mask | substrate | fold | E, Eh | ΔE (mEh) | E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | Span (Lowdin) | Elec (Lowdin) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **CASCI** | AO {3d4s4p4d} | — | polarized | raw | pending |  |  |  |  |  |  |  |  |  | pending |
|  |  | — | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  |  | — | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  |  | — | polarized | fractional | pending |  |  |  |  |  |  |  |  |  | pending |
| **CASSCF** | AO {3d4s4p4d} | free | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | virtual-mix |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | free | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | radial | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | radial | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | frozen window | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | frozen window | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
| **SA-CASSCF** | AO {3d4s4p4d} | free | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | free | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | radial | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | radial | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | frozen window | polarized | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | AO {3d4s4p4d} | frozen window | purified | averaged | pending |  |  |  |  |  |  |  |  |  | pending |
|  | PySCF MO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |
|  | NO |  |  |  | pending |  |  |  |  |  |  |  |  |  | pending |

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
