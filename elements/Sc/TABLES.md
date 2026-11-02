# Sc — scandium (d¹)

cc-pVDZ. Reference **ROHF −759.736255 Eh** (3d¹4s²; M_S = ½; S
half-integer). Section grammar and column keys — README · Results.

## Contract — term table

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |

## Production — multiplet energies

Route free · partition averaged fold · source `report.md` (`scripts/run.py Sc`).

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ²D | 5 | −0.342549 | — | −2…2 |

## Production — window lattice

Route free · partition averaged fold · sector (2, 1) · source `cells.json` (as-run).

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|
| 3d | CAS(3,5) | 50 | −759.4359 | 12.000 | −759.7467 | 0.006 |
| 3d4s | CAS(3,6) | 90 | −759.7296 | 6.000 | −759.7470 | 0.006 |
| 3d4d | CAS(3,10) | 450 | −759.5553 | 12.000 | −759.7482 | 1.578 |
| 3d4s4p | CAS(3,9) | 324 | −759.7611 | 6.000 | −759.7481 | 0.627 |
| 3d4s4d | CAS(3,11) | 605 | −759.7402 | 6.000 | −759.7484 | 1.603 |
| 3d4s4p4d | CAS(3,14) | 1,274 | −759.7658 | 6.000 | −759.7510 | 2.868 |
| 3d4s4p4d5d | CAS(3,19) | 3,249 | −759.7660 | 6.000 | −759.7696 | 5.590 |
| 3d4s4p4d4f | CAS(3,21) | 4,410 | −759.7697 | 6.000 | −759.7720 | 5.635 |

## Production — the method ladder, CAS(3,6)/{3d,4s}

Route free · partition the treatment column · source `cells.json` (the committed draw).

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | cross-check (mEh) | note |
|---|---|---|---|---|---|---|
| ROHF (anchor) | −759.736255 | 0 | — | — | 0 | forced filling |
| CASCI, averaged | −759.729565 | +6.7 | 5.9998 | 8.9e-04 | 1.3 | reference fold |
| CASCI, raw | −759.729565 | +6.7 | 6.0003 | 2.8e-02 | ~0 | polarized-core detector |
| CASCI, fractional | −759.699492 | +36.8 | 6.0000 | 3.3e-14 | ~0 | different functional |
| CASSCF, free | −759.767951 | −31.7 | 0.0469 | 5.7e-01 | 1.9 | escape — no term |
| CASSCF, frozen window | −759.729655 | +6.6 | 5.9998 | 8.9e-04 | 1.3 | |
| SA-CASSCF, free | −759.749511 | −13.3 | 6.0020 | 7.9e-01 | 798.6 | swap basin — totals only |
| CASCI, pur | −759.728891 | +7.4 | 6.0000 | 3.4e-14 | 0.0 | pur substrate |
| SA-CASSCF, frozen, pur | −759.728978 | +7.3 | 6.0000 | 6.5e-14 | 0.0 | the certified rung |

## Lane — the pinned ladder (semi-core scaling)

Route pinned · partition the fold column · source `cells.json` (as-run) · pins `test_symm_tm.py::test_cas_ladder_*`.

<!-- lane:begin -->
| window | CAS | fold | dets | E, Eh | vs declared | term |
|---|---|---|---|---|---|---|
| 3d4s | (3,6) | Ar | 90 | −759.7274 | 0 | ²D exact |
| 3d3p | (7,8) | semi | 3,920 | −759.7871 | −59.7 mEh | ²D exact |
| 3d3p4d | (7,13) | semi | 204,490 | −759.8409 | −113.4 mEh | ²D exact |
| 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.7899 | −62.4 mEh | ²D exact |
| 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.7900 | −62.5 mEh | ²D exact |
<!-- lane:end -->

The semi-core window pays −59.7 mEh over the declared window at
frozen shells — the 3p⁶ correlation — with ²D hosted exactly
(components degenerate, every meter exact).  The extension legs place
the pair-correlating 4d shell and the 4s/3s pairs against it.

## Readings — open questions

Semicore measured (−59.7 mEh, the ladder pair; d-count-independent);
spin–orbit of ²D; the frac-occ septet ordering is
SCF-chaos-dependent — never pin.

## Studies (frozen history; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)

| cell | status | E, Eh | ⟨L²⟩ | ‖[P,Lz]‖ | ‖[P,Lp]‖ |
|---|---|---|---|---|---|
| study_sc_no_no2_casci | ok | -759.730616 | 6.0615 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no2_relax | ok | -759.731210 | 6.0602 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no2_sa-casscf | ok | -759.749511 | 8.9253 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no3_casci | ok | -759.762160 | 6.5059 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no3_relax | ok | -759.762734 | 6.5047 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no3_sa-casscf | ok | -759.765718 | 6.3497 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_casci | ok | -759.765064 | 6.6032 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_relax | ok | -759.765550 | 6.6014 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_sa-casscf | ok | -759.771003 | 6.4999 | 2.2e-13 | 3.5e-13 |

<!-- study:end -->
