# C — carbon (p²)

cc-pVDZ. Reference **ROHF −37.682418 Eh** (raw core = averaged fold).
Sectors Ms = 1 and Ms = 0. Provenance: `results/C_ccpvdz_s2_ladder.json`
(T3b values), campaign logs (T1/T2/T5). NOTES.md: the confinement
finding, the flat-gauge census, the free-CASSCF tilt anatomy.

## T1 — Term content (the contract)

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

## T3b — Treatment ladder at {2p} (single window; cross-check ≈ 0:
the 2p shell is exactly S-orthogonal to the core — no window offset)

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | note |
|---|---|---|---|---|---|
| ROHF (anchor) | −37.682418 | 0 | — | — | |
| CASCI, averaged | −37.431611 | +250.8 | 2.0000 | 1.3e-15 | the confinement cost (NOTES) |
| CASCI, raw | −37.431611 | +250.8 | 2.0000 | 9.8e-06 | same run; the polarized-core detector |
| CASCI, fractional | −37.440244 | +242.2 | 2.0000 | 0.0 | a different functional |
| CASSCF, free | −37.682418 | +0.0 | 1.6864 | 2.3e-02 | member-tilted gauge (NOTES) |
| CASSCF, core-relax | −37.440394 | +242.0 | 2.0000 | 7.6e-17 | |
| SA-CASSCF, free | −37.682391 | +0.0 | 2.0000 | 1.0e-10 | lands symmetric |
| SA-CASSCF, core-relax, pur | −37.440382 | +242.0 | 2.0000 | 0.0 | the certified rung |

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| recovers everything (θ* = 0.61, −250.8 mEh) | exactly flat |
