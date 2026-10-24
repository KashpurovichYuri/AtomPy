# Sc — scandium (d¹)

cc-pVDZ. Reference **ROHF −759.736255 Eh** (3d¹4s²; M_S = ½; S
half-integer). Provenance: `results/Sc_ccpvdz_s1_ladder.json` (T3b),
the lattice logs (T3a), `results/tm_no.log`. NOTES.md: escape/swap
basins, bimodal SA, frac septet, AVAS cross-check.

## T1 — Term content (the contract)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |

## T2 — Spectrum (CAS(1,5), averaged fold)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ²D | 5 | −0.342549 | — | −2…2 |

## T3a — Window lattice at matched cost (frozen CASCI; ROHF reference
−759.736 Eh; the d-atom window already correlates)

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|
| 3d | (3,5) | 50 | −759.436 | 12.000 | — | — | — | — |
| 3d4s | (3,6) | 90 | −759.729 | 6.000 | −759.736 | 0.000 | −759.731 | 6.000 |
| 3d4d | (3,10) | 450 | −759.555 | 12.000 | — | — | — | — |
| 3d4s4p | (3,9) | 324 | −759.760 | 6.000 | — | — | −759.762 | 6.000 |
| 3d4s4d | (3,11) | 605 | −759.739 | 6.000 | −759.745 | 0.008 | — | — |
| 3d4s4p4d | (3,14) | 1,274 | **−759.765** | 6.000 | −759.745 | 0.143 | **−759.765** | 6.000 |
| 3d4s4p4d5d | (3,19) | 3,249 | −759.765 | 6.000 | — | — | — | — |
| 3d4s4p4d4f | (3,21) | 4,410 | **−759.769** | 6.000 | — | — | — | — |

## T3b — Treatment ladder at 3d4s (cross-check = the averaging effect
(1.3 mEh; ticket 25) — the old ~90 mEh window-definition offset was an
artifact of the retired raw-AO-span analysis convention; 0 on the pur
rows)

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | cross-check (mEh) | note |
|---|---|---|---|---|---|---|
| ROHF (anchor) | −759.736255 | 0 | — | — | 0 | forced filling |
| CASCI, averaged | −759.729565 | +6.7 | 6.0000 | 5.6e-15 | 1.3 | reference fold |
| CASCI, raw | −759.729565 | +6.7 | 6.0004 | 2.8e-02 | ~0 | polarized-core detector |
| CASCI, fractional | −759.699492 | +36.8 | 6.0000 | 3.4e-14 | ~0 | different functional |
| CASSCF, free | −759.767769 | −31.5 | 0.0526 | 3.2e-01 | 1.9 | escape — no term |
| CASSCF, core-relax | −759.729655 | +6.6 | 5.9998 | 2.8e-03 | 1.3 | |
| SA-CASSCF, free | −759.749511 | −13.3 | 6.0020 | 8.3e-01 | 798.6 | swap basin — totals only |
| CASCI, pur | −759.728891 | +7.4 | 6.0000 | 6.6e-15 | 0.0 | pur substrate |
| SA-CASSCF, core-relax, pur | −759.728978 | +7.3 | 6.0000 | 3.5e-14 | 0.0 | the certified rung |

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| rises everywhere | descends (−0.3 mEh @ δ = 0.1; −0.9 with four members) — JT-like |

Ticket-24 amendment: the +25 mEh radial ceiling is a window statement
(the lattice crosses below the escape via 4f, −3.06 mEh).

## Open questions

Semicore (3p; deprioritized); spin–orbit of ²D; the frac-occ septet
ordering is SCF-chaos-dependent — never pin.
