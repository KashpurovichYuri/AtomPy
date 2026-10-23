# Sc — scandium (d¹)

cc-pVDZ (Balabanov–Peterson). Reference **ROHF −759.736255 Eh**
(3d¹4s²; M_S = ½; S half-integer — snap to the half-grid). Provenance:
`results/sc_ladder.md`, `results/*_ladder.md`, `results/tm_no.log`,
`tests/test_ladder.py` (Sc block). NOTES.md carries the basin stories.

## T1 — Term content (the contract: what exact diagonalization must return)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |

CAS(1,5) is the Sc⁺-like d¹ space on the [Ar]+4s² closed core.

## T2 — Spectrum (CAS(1,5), averaged fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ²D | 5 | −0.342549 | — | −2…2 |

## T3 — Window × treatment (energies in Eh; ΔE vs ROHF in mEh;
‖·‖ = many-body commutator; |[P,·]| = span certificate)

Optimizer treatments measured at the standard window **3d4s** (bold);
"—" = not run (ticket 20 completes the product).

| window | dets | \|[P,Lz]\| | \|[P,L⊥]\| | CASCI, Eh | ⟨L²⟩ | relax | SA |
|---|---|---|---|---|---|---|---|
| 3d | 50 | ~1e-14 | ~1e-14 | −759.436 | 12.000 | — | — |
| **3d4s** | 90 | 7e-15 | 9e-15 | −759.729 (+6.7 avg; raw fold ‖[H,Lz]‖ 2.8e-2) | 6.000 | +6.6, ‖·‖ 5.7e-3 (start-window defect) | free: swap basin (totals −13.3, ‖·‖ 0.8); frozen pur: +7.3 exact; radial pur: 0.0 (floor) |
| 3d4d | 450 | ~1e-14 | ~1e-14 | −759.555 | 12.000 | — | — |
| 3d4s4p | 324 | ~1e-14 | ~1e-14 | −759.760 | 6.000 | — | — |
| 3d4s4d | 605 | ~1e-14 | ~1e-14 | −759.739 | 6.000 | — | — |
| 3d4s4p4d | 1,274 | ~1e-14 | ~1e-14 | **−759.765** | 6.000 | — | — |
| 3d4s4p4d5d | 3,249 | ~1e-14 | ~1e-14 | −759.765 | 6.000 | — | — |
| 3d4s4p4d4f | 4,410 | ~1e-14 | ~1e-14 | **−759.769** | 6.000 | — | — |
| ROHF MO window (open+virtuals) | — | O(1), never closed | O(1) | ⟨L²⟩ 0.000–0.143 | broken | — | — |
| NO groups (term-averaged) | 90–324 | ~1e-13 | ~1e-13 | −759.731 … −759.765 | 6.000 exact | — | — |

Free single-root CASSCF at the standard window: escape basin, −31.5
mEh, L ≈ 0 window, identity-destroying (⟨L²⟩ 0.05, certs 2.3/2.7).
See NOTES.

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| rises everywhere (3d→4d, 4s→5s/6s harmful) | descends: −0.3 mEh @ δ = 0.1, −0.9 with four members — JT-like |

Ticket-24 amendment: the +25 mEh radial ceiling is a window statement —
the grown lattice crosses below the escaped value (4f −3.06 mEh).

## Open questions

- Semicore (3p) treatment (deprioritized); spin–orbit of ²D.
- The frac-occ septet ordering is SCF-chaos-dependent — never pin.
