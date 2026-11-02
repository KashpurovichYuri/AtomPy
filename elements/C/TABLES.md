# C — carbon (p²)

cc-pVDZ. Reference **ROHF −37.682418 Eh** (raw core = averaged fold).
Sectors Ms = 1 and Ms = 0. Provenance: `report.md/json` (term tables),
`scan.md`, `sa_decomposition.md`, `basin_map.md`; the Methods-compared
regeneration (`scripts/regen_cells.py C`; `cells.json/md`, key-addressed).
Records: the confinement finding and tilt anatomy — README ·
Motivation; the flat-gauge census — AGENTS fact 16a.

## Term table — the combinatorial contract

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {2p} | Ms=1 | 3 | ³P: 3 |
| CAS(2,3) / {2p} | Ms=0 | 9 | ³P: 3, ¹D: 5, ¹S: 1 |

## Multiplet energies (averaged fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³P | 3 | −1.100609 | — | −1…1 |
| ¹D | 5 | −1.012515 | +0.0881 | −2…2 |
| ¹S | 1 | −0.880374 | +0.2202 | 0 |

## The pinned lane — the declared window's term anchor

The pinned route's term energy of the declared window at the uniform
fold (the lane construction; certifying pins in
`test_symm_lane.py`/`test_symm_tm.py`; the block between the markers
is generated from `cells.json` — as-run):

<!-- lane:begin -->
| window | CAS | fold | dets | E, Eh | vs declared | term |
|---|---|---|---|---|---|---|
| 2p | (2,3) | {1s,2s} | 3 | −37.0678 | 0 | ³P exact |
<!-- lane:end -->

## Methods compared at {2p} (single window; cross-check ≈ 0:
the 2p shell is exactly S-orthogonal to the core — no window offset)

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | note |
|---|---|---|---|---|---|
| ROHF (anchor) | −37.682418 | 0 | — | — | |
| CASCI, averaged | −37.431611 | +250.8 | 2.0000 | 1.3e-15 | the confinement cost |
| CASCI, raw | −37.431611 | +250.8 | 2.0000 | 9.8e-06 | same run; the polarized-core detector |
| CASCI, fractional | −37.440244 | +242.2 | 2.0000 | 0.0 | a different functional |
| CASSCF, free | −37.682418 | +0.0 | 1.6864 | 2.3e-02 | member-tilted gauge |
| CASSCF, frozen window | −37.440394 | +242.0 | 2.0000 | 7.6e-17 | |
| SA-CASSCF, free | −37.682391 | +0.0 | 2.0000 | 1.0e-10 | lands symmetric |
| SA-CASSCF, frozen, pur | −37.440382 | +242.0 | 2.0000 | 0.0 | the certified rung |

### Basin map — start dependence (free CASSCF, CAS(2,3)/{2p}, Ms=1)

| start | E, Eh | ⟨L²⟩ | ‖[H,Lz]‖ |
|---|---|---|---|
| 0 (default) | −37.68241788 | 1.6864 | 2.3e-02 |
| 1–4 (perturbed) | −37.68241788 | 1.62–1.66 | 1.1e-01…1.4e-01 |

All five endpoints land on one energy — spread 0.000 mEh: the end
answer is start-independent within a basin. The orbitals are not: each
start draws its own flat-direction gauge (pairwise span overlaps
0.89–0.98, ⟨L²⟩ never 2 — every endpoint is m-tilted and flagged by
the meters). This is why the protocol certifies (SA + frozen window +
purified core) instead of trusting convergence. Runner: the retired
`scripts/basin_map.py` (git history).

## Frozen-scan verdicts (Motivation; curves in `scan.md`)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| recovers everything (θ* = 0.61, −250.8 mEh) | exactly flat |

## Studies (generated between the markers; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)
<!-- study:end -->
