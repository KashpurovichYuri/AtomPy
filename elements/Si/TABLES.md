# Si — silicon (p²)

cc-pVDZ. Reference **ROHF −288.846437 Eh**. Sectors Ms = 1 / Ms = 0 —
as carbon throughout, with its own numbers. Provenance:
`results/Si_ccpvdz_s2_ladder.md`, `results/*_ladder.md`,
`tests/test_ladder.py` (Si block), the ticket-09 guard test.
NOTES.md carries the multi-solution record.

## T1 — Term content (the contract: what exact diagonalization must return)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {3p} | Ms=1 / Ms=0 | 3 / 9 | ³P: 3, ¹D: 5, ¹S: 1 (as carbon) |

## T2 — Spectrum

As carbon's pattern (³P < ¹D < ¹S); values in `results/*_report.md`
(the README §7.2 carbon rows are the pinned reference set).

## T3 — Window × treatment (single window {3p}, raw shell span
certificate ~1e-14; ΔE vs ROHF in mEh; full diagnostics in
`results/Si_ccpvdz_s2_ladder.md`)

| rung | ΔE | ⟨L²⟩ | note |
|---|---|---|---|
| CASCI, averaged fold | +221.5 | 2.000 exact | the confinement cost |
| CASCI, raw fold | +221.5 | 2.000 | ‖[H,Lz]‖ 5.6e-4 — the polarized-core detector |
| CASCI, frac reference | +217.8 | — | different functional |
| CASSCF, free | +0.1 | ≈ 2 | lands symmetric |
| CASSCF, window frozen | +217.5 | — | |
| CASSCF, radial restriction | +0.1 | 2.000 exact | the whole gain is the θ\* rotation |
| SA-CASSCF, free | +0.04 | ≈ 2 | lands symmetric |
| SA-CASSCF, frozen, pur core | +217.7 | — | pur family |
| SA-CASSCF, radial, pur core | +0.04 | 2.000 exact | the certified treatment |

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| gains everything: θ\* = 0.85 (3p→4p), −126.7 mEh | exactly flat |

## Window-definition note

max|dE| ≈ 9.5e-2 Eh for CASCI rungs vs the raw-window definition is the
raw-vs-orthogonalized window offset, not fold anisotropy; the CASSCF
input-window baseline ~1e-3 (the fact-17 class) is removed by the pur
substrate.
