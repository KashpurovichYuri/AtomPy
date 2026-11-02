# Sc — scandium (d¹)

cc-pVDZ. Reference **ROHF −759.736255 Eh** (3d¹4s²; M_S = ½; S
half-integer). Section grammar and column keys — README · Results.

## Term table

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |

## Multiplet energies

Route free · partition averaged fold · source `report.md` (`scripts/run.py Sc`).

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ²D | 5 | −0.342549 | — | −2…2 |

## CASCI — window survey

CI level · ΔE anchored on 3d4s pinned AO (the certified leg) · source `cells.json` (as-run) · pins `test_symm_tm.py::test_cas_ladder_*`.

<!-- survey:begin -->
| window | CAS | fold | start | route | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3d | (3,5) | Ar | AO | free | 50 | −759.4359 | +291.6 | −20665.3136 | +7.9335 | 11.9999 | 7.0e-04 | 3.1e-02 |
| 3d | (3,5) | Ar | MO | free | 50 | −759.7467 | −19.3 | −20673.7714 | −0.5242 | 0.0065 | 1.1e-01 | 2.2e-01 |
| 3d4s | (3,6) | Ar | AO | free | 90 | −759.7296 | −2.1 | −20673.3051 | −0.0579 | 5.9998 | 8.9e-04 | 4.0e-02 |
| 3d4s | (3,6) | Ar | MO | free | 90 | −759.7470 | −19.6 | −20673.7795 | −0.5324 | 0.0062 | 1.5e-01 | 7.2e-01 |
| 3d4s | (3,6) | Ar·raw | AO | free | 90 | −759.7296 | −2.1 | −20673.3051 | −0.0579 | 6.0003 | 2.8e-02 | 1.3e-01 |
| 3d4s | (3,6) | Ar·pur | AO | free | 90 | −759.7289 | −1.5 | −20673.2867 | −0.0396 | 6.0000 | 3.4e-14 | 1.7e-13 |
| 3d4s | (3,6) | Ar·frac | AO | free | 90 | −759.6995 | +27.9 | −20672.4868 | +0.7604 | 6.0000 | 3.3e-14 | 3.1e-13 |
| 3d4s | (3,6) | Ar | AO | pinned | 90 | −759.7274 | 0 | −20673.2471 | +0.0000 | 6.0000 | 3.4e-14 | 3.7e-13 |
| 3d4d | (3,10) | Ar | AO | free | 450 | −759.5553 | +172.1 | −20668.5639 | +4.6833 | 11.9999 | 2.2e-03 | 1.2e-01 |
| 3d4d | (3,10) | Ar | MO | free | 450 | −759.7482 | −20.8 | −20673.8133 | −0.5662 | 1.5780 | 4.7e+00 | 1.3e+01 |
| 3d4s4p | (3,9) | Ar | AO | free | 324 | −759.7611 | −33.7 | −20674.1629 | −0.9157 | 5.9998 | 1.4e-03 | 5.8e-02 |
| 3d4s4p | (3,9) | Ar | MO | free | 324 | −759.7481 | −20.7 | −20673.8103 | −0.5632 | 0.6272 | 3.7e-01 | 8.2e+00 |
| 3d4s4d | (3,11) | Ar | AO | free | 605 | −759.7402 | −12.8 | −20673.5943 | −0.3472 | 5.9998 | 2.5e-03 | 1.3e-01 |
| 3d4s4d | (3,11) | Ar | MO | free | 605 | −759.7484 | −21.0 | −20673.8183 | −0.5711 | 1.6033 | 5.8e+00 | 2.6e+01 |
| 3d4s4p4d | (3,14) | Ar | AO | free | 1,274 | −759.7658 | −38.4 | −20674.2909 | −1.0437 | 5.9998 | 3.3e-03 | 1.6e-01 |
| 3d4s4p4d | (3,14) | Ar | MO | free | 1,274 | −759.7510 | −23.5 | −20673.8875 | −0.6404 | 2.8678 | 6.9e+00 | 2.8e+01 |
| 3d4s4p4d5d | (3,19) | Ar | AO | free | 3,249 | −759.7660 | −38.5 | −20674.2961 | −1.0489 | 5.9998 | 4.2e-03 | 2.1e-01 |
| 3d4s4p4d5d | (3,19) | Ar | MO | free | 3,249 | −759.7696 | −42.2 | −20674.3942 | −1.1471 | 5.5905 | 1.1e+01 | 4.9e+01 |
| 3d4s4p4d4f | (3,21) | Ar | AO | free | 4,410 | −759.7697 | −42.3 | −20674.3985 | −1.1513 | 5.9998 | 1.4e-01 | 8.0e-01 |
| 3d4s4p4d4f | (3,21) | Ar | MO | free | 4,410 | −759.7720 | −44.6 | −20674.4604 | −1.2133 | 5.6353 | 1.2e+01 | 7.2e+01 |
| 3d3p | (7,8) | semi | AO | pinned | 3,920 | −759.7871 | −59.7 | −20674.8712 | −1.6241 | 6.0000 | 7.1e-13 | 4.9e-12 |
| 3d3p4d | (7,13) | semi | AO | pinned | 204,490 | −759.8409 | −113.4 | −20676.3338 | −3.0867 | 6.0000 | — | — |
| 3d3p4s | (9,9) | {1s,2s,2p,3s} | AO | pinned | 15,876 | −759.7899 | −62.4 | −20674.9462 | −1.6990 | 6.0000 | — | — |
| 3s3p3d | (9,9) | {1s,2s,2p,4s} | AO | pinned | 15,876 | −759.7900 | −62.5 | −20674.9492 | −1.7021 | 6.0000 | — | — |
<!-- survey:end -->

The semi-core window pays −59.7 mEh over the declared window at
frozen shells — the 3p⁶ correlation — with ²D hosted exactly
(components degenerate, every meter exact).  The extension legs place
the pair-correlating 4d shell and the 4s/3s pairs against it.

## Method ladder

ΔE vs ROHF anchor (−759.736255 Eh) · source `cells.json` (the committed draw).

| method | CAS | fold | route | E, Eh | ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | −759.736255 | 0 | −20673.4871 | +0.0000 | — | — | — | forced filling |
| CASCI, averaged | (3,6) | Ar | free | −759.7296 | +6.7 | −20673.3051 | +0.1820 | 5.9998 | 8.9e-04 | 4.0e-02 | reference rung |
| CASSCF, free | (3,6) | Ar | free | −759.7680 | −31.7 | −20674.3496 | −0.8625 | 0.0469 | 5.7e-01 | 1.2e+00 | escape — no term |
| CASSCF, frozen | (3,6) | Ar | free | −759.7297 | +6.6 | −20673.3075 | +0.1796 | 5.9998 | 8.9e-04 | 4.0e-02 |  |
| SA-CASSCF, free | (3,6) | Ar | free | −759.7495 | −13.3 | −20673.8478 | −0.3607 | 6.0020 | 7.9e-01 | 3.8e+00 | swap basin — totals only |
| SA-CASSCF, frozen, pur | (3,6) | Ar | free | −759.7290 | +7.3 | −20673.2891 | +0.1980 | 6.0000 | 6.5e-14 | 2.1e-13 | certified rung |
| SA-CASSCF | (7,13) | semi | pinned | pending | | | | | | | the route axis, declared |

## Studies (frozen history; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)

| cell | status | E, Eh | ⟨L²⟩ | ‖[P,Lz]‖ | ‖[P,Lp]‖ |
|---|---|---|---|---|---|
| study_sc_no_no2_casci | ok | -759.7306 | 6.0615 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no2_relax | ok | -759.7312 | 6.0602 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no2_sa-casscf | ok | -759.7495 | 8.9253 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no3_casci | ok | -759.7622 | 6.5059 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no3_relax | ok | -759.7627 | 6.5047 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no3_sa-casscf | ok | -759.7657 | 6.3497 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_casci | ok | -759.7651 | 6.6032 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_relax | ok | -759.7655 | 6.6014 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_sa-casscf | ok | -759.7710 | 6.4999 | 2.2e-13 | 3.5e-13 |

<!-- study:end -->
