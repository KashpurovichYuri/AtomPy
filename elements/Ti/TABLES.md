# Ti — titanium (d²)

cc-pVDZ (Balabanov–Peterson). Reference **ROHF −848.406489 Eh** (3d²4s²,
gated good-basin; M_S = 1). Ground term ³F. Provenance:
`results/ti_ladder.md`, `results/*_ladder.md`, `results/tm_no.log`,
`tests/test_ladder.py` (Ti block), `TI_EXACT_TERMS` in
`tests/conftest.py`. NOTES.md carries the basin stories.

## T1 — Term content (the contract: what exact diagonalization must return)

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

## T3 — Window × treatment (energies in Eh; ΔE vs ROHF in mEh;
‖·‖ = many-body commutator of the cell's fold; |[P,·]| = the window
span certificate)

Rows = windows; optimizer treatments measured at the standard window
**3d4s** (bold) and at the NO groups; "—" = not run (the window ×
treatment product is ticket 20's completion).

| window | dets | \|[P,Lz]\| | \|[P,L⊥]\| | CASCI, Eh | ⟨L²⟩ | relax | SA |
|---|---|---|---|---|---|---|---|
| 3d | 50 | ~1e-14 | ~1e-14 | −848.054 | 6.000 | — | — |
| **3d4s** | 120 | 6e-15 | 8e-15 | −848.395 (+10.6 avg; raw fold ‖[H,Lz]‖ 4.9e-2) | 12.000 | +10.4, ‖·‖ 2.3e-8 | free: swap basin (totals −10.1, ‖·‖ 1.1); frozen pur: +11.7 exact |
| 3d4d | 1,200 | ~1e-14 | ~1e-14 | −848.232 | 6.000 | — | — |
| 3d4s4p | 756 | ~1e-14 | ~1e-14 | −848.429 | 12.000 | — | — |
| 3d4s4d | 1,815 | ~1e-14 | ~1e-14 | −848.409 | 12.000 | — | — |
| 3d4s4p4d | 5,096 | ~1e-14 | ~1e-14 | **−848.438** | 12.000 | — | — |
| 3d4s4p4d5d | 18,411 | ~1e-14 | ~1e-14 | −848.439 | 12.000 | — | — |
| 3d4s4p4d4f | 27,930 | ~1e-14 | ~1e-14 | **−848.446** | 12.000 | — | — |
| ROHF MO window (open+virtuals) | — | O(1), never closed | O(1) | ⟨L²⟩ 0.000–0.2 at every size | broken | — | — |
| NO groups (term-averaged) | 120–756 | ~1e-13 | ~1e-13 | −848.399 … −848.432 | 12.000 exact | — | — |

Free single-root CASSCF at the standard window: escape basin, −31.3 mEh,
identity-destroying (⟨L²⟩ 0.05, span certs 2.8/4.0) — broken at every
window where tried; see NOTES.

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| rises (+9.8 mEh @ θ = 0.05); the ~11 mEh radial gain lives on the orthogonalized substrate's s slot | descends (−163 µEh @ δ = 0.1; the pair {0,1} doubles it; {0,1,2} and {0,1,3,4} rise) — JT-like |

Ticket-24 amendment: the radial ceilings are window statements — the
grown m-complete lattice crosses below the escaped CASSCF values (5d
−4.36, 4f −11.46 mEh).

## Open questions

- ⟨L²⟩ = 0 MO-window phenomenon (same class as the U ⁵S° reference;
  mechanism open).
- Semicore (3p) treatment (deprioritized); spin–orbit splitting.
