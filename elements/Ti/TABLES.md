# Ti — titanium (d²)

cc-pVDZ. Reference **ROHF −848.406489 Eh** (3d²4s², gated good-basin).
M_S = 1. Ground term ³F. Provenance: `results/Ti_ccpvdz_s2_ladder.json`
(T3b — all values verified in-session 2026-10-24), the lattice logs
(T3a). NOTES.md: the gated reference, escape/swap, JT, ⟨L²⟩ = 0 MO
windows.

## T1 — Term content (the contract)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| CAS(2,5) / {3d} | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 |
| CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 |

## T2 — Spectrum (CAS(2,5), averaged fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³F | 7 | −1.468663 | — | −3…3 |
| ¹D | 5 | −1.423747 | +0.0449 | −2…2 |
| ³P | 3 | −1.413969 | +0.0547 | −1…1 |
| ¹G | 9 | −1.398223 | +0.0704 | −4…4 |
| ¹S | 1 | −1.295050 | +0.1736 | 0 |

## T3a — Window lattice at matched cost (frozen CASCI; ROHF reference
−848.406 Eh; the d-atom window already correlates)

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|
| 3d | (4,5) | 50 | −848.054 | 6.000 | — | — | — | — |
| 3d4s | (4,6) | 120 | −848.395 | 12.000 | −848.405 | 0.000 | −848.399 | 12.000 |
| 3d4d | (4,10) | 1,200 | −848.232 | 6.000 | — | — | — | — |
| 3d4s4p | (4,9) | 756 | −848.429 | 12.000 | — | — | −848.432 | 12.000 |
| 3d4s4d | (4,11) | 1,815 | −848.409 | 12.000 | −848.414 | 0.010 | — | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.420 | 0.201 | **−848.438** | 12.000 |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.439 | 12.000 | — | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — | — |

## T3b — Treatment ladder at 3d4s (cross-check = the averaging effect
(2.8 mEh; ticket 25: the fold convention is machine-exact) — the old
~88 mEh window-definition offset was an artifact of the retired
raw-AO-span analysis convention; 0 on the pur rows)

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | cross-check (mEh) |
|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | −848.406489 | 0 | — | — | — | — | — | 0 |
| CASCI, averaged | −848.395940 | +10.5 | 12.0000 | 2.2e-14 | 4.1e-13 | 5.9e-15 | 7.6e-15 | 2.8 |
| CASCI, raw | −848.395940 | +10.5 | 12.0000 | 4.9e-02 | 2.6e-01 | 5.9e-15 | 7.6e-15 | ~0 |
| CASCI, fractional | −848.358519 | +48.0 | 12.0000 | 1.2e-13 | 7.5e-13 | 5.9e-15 | 7.6e-15 | ~0 |
| CASSCF, free | −848.437740 | −31.3 | 0.0503 | 4.3e-01 | 1.1e+00 | 2.8 | 4.0 | 2.8 |
| CASSCF, core-relax | −848.396109 | +10.4 | 12.0 | 2.3e-08 | 1.7e-07 | 3.2e-08 | 5.9e-08 | 2.7 |
| SA-CASSCF, free | −848.416563 | −10.1 | 12.003 | 1.1e+00 | 6.0e+00 | 1.4 | 1.5 | 851 |
| SA-CASSCF, core-relax, pur | −848.394753 | +11.7 | 12.0 | 5.2e-13 | 2.9e-12 | 7.2e-15 | 1.2e-14 | 0.0 |

Legend: E (as run) = the run's own variational total; ΔE vs the ROHF
anchor; cross-check = |E(our analysis fold) − E(as run)| max over
roots (ticket 25: the averaging effect + the window-definition
offset; the fold convention itself is machine-exact). Rows whose
window identity changed (escape, swap) are different problems —
totals only.

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| rises (+9.8 mEh @ θ = 0.05); the ~11 mEh gain lives on the orthogonalized s slot | descends (−163 µEh @ δ = 0.1; pair patterns vary) — JT-like |

Ticket-24 amendment: the ceilings are window statements (5d −4.36,
4f −11.46 mEh cross below the escape).

## Open questions

⟨L²⟩ = 0 MO-window phenomenon (the U ⁵S° class; mechanism open);
semicore (deprioritized); spin–orbit.
