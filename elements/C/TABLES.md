# C — carbon (p²)

cc-pVDZ. Reference **ROHF −37.682418 Eh**. Sectors Ms = 1 and Ms = 0.
The textbook atom: ³P < ¹D < ¹S. Provenance: `results/*_ladder.md`,
`tests/test_ladder.py` (C block), `tests/test_full_multiplet_ms0`.
NOTES.md carries the confinement finding and the flat-gauge census.

## T1 — Term content (the contract: what exact diagonalization must return)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {2p} | Ms=1 | 3 | ³P: 3 |
| CAS(2,3) / {2p} | Ms=0 | 9 | ³P: 3, ¹D: 5, ¹S: 1 |

## T2 — Spectrum (averaged fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³P | 3 | −1.100609 | — | −1…1 |
| ¹D | 5 | −1.012515 | +0.0881 | −2…2 |
| ¹S | 1 | −0.880374 | +0.2202 | 0 |

## T3 — Window × treatment (single window {2p}, raw shell span
certificate ~1e-14; ΔE vs ROHF in mEh; full diagnostics in
`results/*_ladder.md`)

| rung | ΔE | ⟨L²⟩ | note |
|---|---|---|---|
| CASCI, averaged fold | +250.8 | 2.000 exact | the confinement cost (NOTES) |
| CASCI, raw fold | +250.8 | 2.0000000000 | ‖[H,Lz]‖ 9.8e-6 — the polarized-core detector |
| CASCI, frac reference | +242.2 | — | different functional |
| CASSCF, free | +0.05 | 1.69–1.72 | member-tilted gauge (NOTES); ‖[P,L⊥]‖ ≈ 1.1 |
| CASSCF, window frozen | +242.0 | — | core relaxation only |
| CASSCF, radial restriction | +0.05 | 2.000 exact | the entire gain is the θ\* rotation |
| SA-CASSCF, free | +0.03 | ≈ 2 | lands symmetric |
| SA-CASSCF, frozen, pur core | +242.0 | — | pur family |
| SA-CASSCF, radial, pur core | +0.03 | 2.000 exact | the certified treatment |

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| recovers everything: θ\* = 0.61 (2p→3p), −250.8 mEh | exactly flat (five draws, one energy) |
