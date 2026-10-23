# C — carbon (p²)

cc-pVDZ; reference **ROHF −37.682418 Eh** (raw core = averaged fold).
Standard sectors Ms = 1 and Ms = 0. The textbook atom: ³P < ¹D < ¹S,
the confinement finding, and the full m-symmetric radial recovery.

Provenance: `results/*_ladder.md`, `tests/test_ladder.py` (C block),
`tests/test_radial_rung_full_gain_without_l_leakage`,
`tests/test_full_multiplet_ms0`. Cross-atom verdicts: README §7.3–7.5.

## Term content

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {2p} | Ms=1 | 3 | ³P: 3 |
| | Ms=0 | 9 | ³P: 3, ¹D: 5, ¹S: 1 |

## Multiplet spectrum (our fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³P | 3 | −1.100609 | — | −1…1 |
| ¹D | 5 | −1.012515 | +0.0881 | −2…2 |
| ¹S | 1 | −0.880374 | +0.2202 | 0 |

## Treatment ladder (ΔE vs ROHF, mEh)

| rung | ΔE | note |
|---|---|---|
| CASCI, averaged fold | +250.8 | exact — the CONFINEMENT cost |
| CASCI, raw core fold | +250.8 | comm 9.8e-6 |
| CASCI, frac reference | +242.2 | |
| CASSCF, free | +0.05 | member-tilted gauge (chaotic) |
| CASSCF, window frozen | +242.0 | |
| CASSCF, radial restriction | +0.05 | exact: the ENTIRE gain is the θ\* rotation |
| SA-CASSCF, free | +0.03 | lands symmetric |
| SA-CASSCF, frozen, purified core | +242.0 | |
| SA-CASSCF, radial, purified core | +0.03 | exact |

## The confinement finding

CASCI over the ROHF's own 2p-like orbitals (canonical window)
reproduces the ROHF total exactly (ΔE = 0.00 mEh): the ³P pair carries
no net correlation the CAS can express. The +250.8 mEh "CASSCF gain"
is the recovery of the raw-shell confinement cost, not correlation won.

## Scan anatomy (§7.4)

The whole gain is ONE m-symmetric radial rotation θ\* = 0.61 over
2p → 3p (−250.8 mEh) at exact symmetry ([H,Lz] = 0 along the curve);
single-member tilts cost EXACTLY 0 (five draws, one energy).

Purified-substrate facts: the pur core closes both input windows at
machine precision; the pur family carries the window-definition offset
vs raw-fold CASCI rows (bound, don't pin).
