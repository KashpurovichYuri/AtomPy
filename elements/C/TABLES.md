# C — carbon (p²)

cc-pVDZ. Reference **ROHF −37.6824 Eh** (raw core = averaged fold).
Sectors Ms = 1 and Ms = 0. Section grammar and column keys —
README · Results.

## Term table

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {2p} | Ms=1 | 3 | ³P: 3 |
| CAS(2,3) / {2p} | Ms=0 | 9 | ³P: 3, ¹D: 5, ¹S: 1 |

## Multiplet energies

Route free · partition averaged fold · source `report.md` (`scripts/run.py C`).

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³P | 3 | −1.100609 | — | −1…1 |
| ¹D | 5 | −1.012515 | +0.0881 | −2…2 |
| ¹S | 1 | −0.8803 | +0.2202 | 0 |

## CAS ladder

Routes in rows · partition fold column · source `cells.json` (as-run) · pins `test_symm_lane.py`.

<!-- lane:begin -->
| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pinned | CI | 2p | (2,3) | {1s,2s} | 3 | −37.0678 | 0 | −1008.6665 | +0.0000 | 2.0000 | 0e+00 | 0.0e+00 | 0.0e+00 |
<!-- lane:end -->

## Method ladder, CAS(2,3)/{2p}

Route free · partition treatment column · source `cells.json` (the committed draw).

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | note |
|---|---|---|---|---|---|
| ROHF (anchor) | −37.6824 | 0 | — | — | |
| CASCI, averaged | −37.431611 | +250.8 | 2.0000 | 1.3e-15 | the confinement cost |
| CASCI, raw | −37.431611 | +250.8 | 2.0000 | 9.8e-06 | polarized-core detector |
| CASCI, fractional | −37.440244 | +242.2 | 2.0000 | 0.0 | a different functional |
| CASSCF, free | −37.682418 | +0.0 | 1.6864 | 2.3e-02 | member-tilted gauge |
| CASSCF, frozen window | −37.440394 | +242.0 | 2.0000 | 7.6e-17 | |
| SA-CASSCF, free | −37.682391 | +0.0 | 2.0000 | 1.0e-10 | lands symmetric |
| SA-CASSCF, frozen, pur | −37.440382 | +242.0 | 2.0000 | 0.0 | the certified rung |

## Basin map

Free CASSCF start dependence (CAS(2,3)/{2p}, Ms=1) · source the session record (the retired `scripts/basin_map.py`).

| start | E, Eh | ⟨L²⟩ | ‖[H,Lz]‖ |
|---|---|---|---|
| 0 (default) | −37.68241788 | 1.6864 | 2.3e-02 |
| 1–4 (perturbed) | −37.68241788 | 1.62–1.66 | 1.1e-01…1.4e-01 |

All five endpoints land on one energy (spread 0.000 mEh) with
pairwise span overlaps 0.89–0.98 and ⟨L²⟩ never 2: each start draws
its own flat-direction gauge — every endpoint is m-tilted and the
meters flag it. This is why the protocol certifies (SA + frozen
window + purified core) instead of trusting convergence.

## Studies (frozen history; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)
<!-- study:end -->
