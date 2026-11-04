# Ti — titanium (3d²4s²)

cc-pVDZ. Reference **ROHF −848.4064 Eh** (gated good-basin).
M_S = 1. Ground term ³F. Section grammar and column keys —
README · Results.

## Term table

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| CAS(2,5) / {3d} | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 |
| CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 |

## Multiplet energies

Route free · partition averaged fold · source `report.md` (`scripts/run.py Ti`).

CASCI, CAS(2,5).

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³F | 7 | −1.4687 | — | −3…3 |
| ¹D | 5 | −1.4237 | +0.0449 | −2…2 |
| ³P | 3 | −1.4140 | +0.0547 | −1…1 |
| ¹G | 9 | −1.3982 | +0.0704 | −4…4 |
| ¹S | 1 | −1.2951 | +0.1736 | 0 |

### One run, three fields — CAS(4,6) Ms=1 per root (the split demonstration)

Partition: one column per fold · source the session record (the retired `scripts/multiplet_split.py`; validated on C and Ti CAS(2,5) — ladder exact).

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
| 0 | −3.3436 | 0 | −3.4621 | 0 | −3.3441 | 3.9036 |
| +1 | −3.343644 | 1 | −3.462172 | 1 | −3.342274 | 5.9020 |
| +2 | −3.343644 | 4 | −3.462172 | 4 | −3.342274 | 0.1927 |
| +3 | −3.343644 | 9 | −3.462172 | 9 | −3.342274 | 5.9020 |
| ⁵F −3 | −3.296365 | 9 | −3.457635 | 9 | −3.297772 | 5.9023 |
| ⁵F −2 | −3.296365 | 4 | −3.457635 | 4 | −3.297772 | 0.1927 |

## Window lattice

Route free · partition averaged fold · sector (2, 0) · source `cells.json` (as-run).

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|
| 3d | CAS(2,5) | 10 | −848.4001 | 12.000 | −848.4065 | 0.000 |
| 3d4p | CAS(2,8) | 28 | −848.4001 | 12.000 | −848.4065 | 3.426 |
| 3d4p4d | CAS(2,13) | 78 | −848.4096 | 12.000 | −848.4065 | 5.150 |

## CAS(4,·) valence lattice

Route free · partition averaged fold · source the lattice logs (session).


| window | CAS | dets | AOs, Eh | ⟨L²⟩ | PySCF MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|---|---|
| 3d4s | (4,6) | 120 | −848.3950 | 12.000 | −848.4050 | 0.000 | −848.3990 | 12.000 | — | — |
| 3d4s4p | (4,9) | 756 | −848.4290 | 12.000 | — | — | −848.4320 | 12.000 | — | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.4200 | 0.201 | **−848.438** | 12.000 | — | — |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.4390 | 12.000 | — | — | — | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — | — | — | — |
| 3d4d | (4,10) constrained | 1,200 | −848.2320 | 6.000 | — | — | — | — | — | — |
| 3d4s4d | (4,11) constrained | 1,815 | −848.4090 | 12.000 | −848.4140 | 0.010 | — | — | — | — |

## CAS ladder

Routes in rows · partition fold column · source `cells.json` (as-run) · pins `test_symm_tm.py::test_cas_ladder_*`.

<!-- lane:begin -->
| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pinned | CI | 3d4s | (4,6) | Ar | 120 | −848.3931 | 0 | −23085.9634 | +0.0000 | — | 12.0000 | 5e-15 | 2.7e-14 | 2.4e-13 |
| pinned | CI | 3d4d | (4,10) | Ar | 1,200 | −848.2156 | +177.5 | −23081.1339 | +4.8296 | — | 6.0–30.0 | 5e-02 | 3.4e-13 | 2.1e-12 |
| pinned | CI | 3d4s4p4d | (4,14) | Ar | 5,096 | −848.4374 | −44.4 | −23087.1708 | −1.2074 | — | 12.0000 | 2e-13 | 6.2e-13 | 3.4e-12 |
| pinned | CI | 3d3p | (8,8) | semi | 3,136 | −848.4520 | −58.9 | −23087.5662 | −1.6028 | — | 12.0000 | 3e-14 | 5.2e-13 | 4.1e-12 |
| pinned | CI | 3d3p4d | (8,13) | semi | 368,082 | −848.4991 | −106.0 | −23088.8489 | −2.8855 | — | 12.0000 | 1e-13 | 2.8e-14 | 2.3e-13 |
| pinned | CI | 3s3p3d | (10,9) | {1s,2s,2p,4s} | 10,584 | −848.4549 | −61.8 | −23087.6452 | −1.6818 | — | 12.0000 | 6e-14 | 7.2e-14 | 3.5e-13 |
| pinned | CI | 3d3p4s | (10,9) | {1s,2s,2p,3s} | 10,584 | −848.4529 | −59.8 | −23087.5918 | −1.6284 | — | 12.0000 | 3e-13 | 2.9e-14 | 1.8e-13 |
<!-- lane:end -->

Readings:
- the semi-core 3p⁶ correlation (−58.9 mEh) outweighs the whole
  valence 4p4d set (−44.4); with the 4d pair-correlating shell added
  (the 3p²→4d² channel) it doubles to −106.0 — the free SA's torn
  escape (−5.2 mEh) was a small corner of this;
- the 4s pair stays in the fold: back in-window it buys 0.9 mEh
  (3d3p4s vs 3d3p);
- the 3s² shell adds 2.9 mEh (3s3p3d vs 3d3p);
- every m-complete window hosts the term exactly — components
  degenerate to ≤3.3e-13 Eh;
- the planned 3d4d window as declared (fold Ar, the 4s pair in the
  dead zone) hosts D and H components with quintet admixture at
  +177.5 mEh — the campaign's constrained (4,10) cells carry their
  own semantics and are not read against this ladder.

## Method ladder, CAS(4,6)/{3d,4s}

Route free · partition the fold column · source `cells.json` (the committed draw).

ROHF anchor: −848.4065 Eh (−23086.3283 eV).

| treatment | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | Span (Lowdin) |
|---|---|---|---|---|---|---|---|---|---|
| CASCI, raw | −848.3959 | +10.5 | −23086.0413 | +0.2871 | — | 12.0000 | 4.9e-02 | 2.6e-01 | s0.17 p0 d0.83 f0 |
| CASCI, averaged | −848.3959 | +10.5 | −23086.0413 | +0.2871 | — | 12.0000 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 |
| CASCI, pur | −848.3946 | +11.9 | −23086.0044 | +0.3239 | — | 12.0000 | 8.7e-14 | 6.9e-13 | s0.17 p0 d0.83 f0 |
| CASCI, fractional | −848.3585 | +48.0 | −23085.0230 | +1.3053 | — | 12.0000 | 1.0e-09 | 7.6e-09 | s0.17 p0 d0.83 f0 |
| CASSCF, free (+mo,no) | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 1e-8/38 | 0.0503 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 |
| CASSCF, virtual-mix | −848.4377 | −31.2 | −23087.1779 | −0.8496 | 1e-9/42 | 0.0496 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 |
| CASSCF, free, pur (+mo,no) | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 1e-11/41 | 0.0503 | 4.3e-01 | 1.1e+00 | s0.17 p0.50 d0.33 f4e-04 |
| CASSCF, radial | −848.4112 | −4.7 | −23086.4576 | −0.1292 | 1e-11/17 | 6.9323 | 3.2e+00 | 1.3e+01 | s0.17 p0 d0.83 f0 |
| CASSCF, radial, pur | −848.4100 | −3.5 | −23086.4234 | −0.0951 | 0/12 | 6.9638 | 3.2e+00 | 1.2e+01 | s0.17 p0 d0.83 f0 |
| CASSCF, frozen (+no) | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 1e-11/5 | 12.0000 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 |
| CASSCF, frozen, MO | −848.4171 | −10.6 | −23086.6177 | −0.2893 | 0/4 | 0.0051 | 1.6e-01 | 3.9e-01 | s0.17 p0.50 d0.33 f0 |
| CASSCF, frozen, pur (+no) | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 1e-12/6 | 12.0000 | 8.9e-14 | 4.7e-13 | s0.17 p0 d0.83 f0 |
| CASSCF, frozen, pur, MO | −848.4171 | −10.6 | −23086.6177 | −0.2893 | 0/5 | 0.0051 | 1.6e-01 | 3.9e-01 | s0.17 p0.50 d0.33 f0 |
| SA-CASSCF, free (+no) | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 1e-9/33 | 12.0029 | 1.1e+00 | 6.0e+00 | s0 p0.17 d0.83 f0 |
| SA-CASSCF, free, MO | −848.4316 | −25.1 | −23087.0118 | −0.6835 | 1e-11/9 | 0.0255 | 3.0e-01 | 7.4e-01 | s0.17 p0.50 d0.33 f0 |
| SA-CASSCF, free, pur (+no) | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 1e-8/94 | 12.0029 | 1.0e+00 | 6.0e+00 | s0 p0.17 d0.83 f0 |
| SA-CASSCF, free, pur, MO | −848.4316 | −25.1 | −23087.0119 | −0.6835 | 1e-11/10 | 0.0255 | 3.0e-01 | 7.4e-01 | s0.17 p0.50 d0.33 f0 |
| SA-CASSCF, radial | −848.4069 | −0.4 | −23086.3399 | −0.0115 | 1e-10/8 | 12.0000 | 1.3e-03 | 6.9e-03 | s0.17 p0 d0.83 f0 |
| SA-CASSCF, radial, pur | −848.4055 | +1.0 | −23086.3017 | +0.0266 | 1e-10/8 | 12.0000 | 1.4e-04 | 5.7e-04 | s0.17 p0 d0.83 f0 |
| SA-CASSCF, frozen (+no) | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 1e-10/6 | 12.0000 | 4.2e-09 | 2.3e-08 | s0.17 p0 d0.83 f0 |
| SA-CASSCF, frozen, MO | −848.4170 | −10.5 | −23086.6151 | −0.2868 | 1e-10/4 | 0.0051 | 1.6e-01 | 3.9e-01 | s0.17 p0.50 d0.33 f0 |
| SA-CASSCF, frozen, pur (+no) | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 1e-11/6 | 12.0000 | 4.8e-14 | 4.0e-13 | s0.17 p0 d0.83 f0 |
| SA-CASSCF, frozen, pur, MO | −848.4170 | −10.5 | −23086.6151 | −0.2868 | 1e-13/6 | 0.0051 | 1.6e-01 | 3.9e-01 | s0.17 p0.50 d0.33 f0 |


### CAS(4,10)/{3d,4d} — declared, gated

ROHF anchor: −848.4065 Eh (−23086.3283 eV).

Declared, not yet run -- the optimizer cells are void until their CASCI anchors are read (the anchor-first scenario). One regen away: `regen_cells.py Ti --window 3d4d`.


### CAS(4,14)/{3d,4s,4p,4d} — declared, gated

ROHF anchor: −848.4065 Eh (−23086.3283 eV).

Declared, not yet run -- the optimizer cells are void until their CASCI anchors are read (the anchor-first scenario). One regen away: `regen_cells.py Ti --window 3d4s4p4d`.

## Studies (frozen history; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)

| cell | status | E, Eh | ⟨L²⟩ | ‖[P,L̂z]‖ | ‖[P,L̂⊥]‖ |
|---|---|---|---|---|---|
| study_ti_no_no2_casci | ok | -848.3988 | 10.6709 | 9.4e-15 | 2.6e-14 |
| study_ti_no_no2_relax | ok | -848.4003 | 10.6939 | 9.4e-15 | 2.6e-14 |
| study_ti_no_no2_sa-casscf | ok | -848.4166 | 11.0461 | 9.4e-15 | 2.6e-14 |
| study_ti_no_no3_casci | ok | -848.4324 | 8.5168 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no3_relax | ok | -848.4338 | 8.5669 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no3_sa-casscf | ok | -848.4366 | 9.2758 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no4_casci | ok | -848.4385 | 8.4219 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no4_relax | ok | -848.4398 | 8.4347 | 2.2e-13 | 4.4e-13 |
| study_ti_no_no4_sa-casscf | ok | -848.4424 | 8.0006 | 2.2e-13 | 4.4e-13 |

<!-- study:end -->
