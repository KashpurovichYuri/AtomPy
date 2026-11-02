# Si — silicon (p²)

cc-pVDZ. Reference **ROHF −288.846437 Eh**. Sectors Ms = 1 / Ms = 0.
Section grammar and column keys — README · Results.

## Contract — term table

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {3p} | Ms=1 | 3 | ³P: 3 |
| CAS(2,3) / {3p} | Ms=0 | 9 | ³P: 3, ¹D: 5, ¹S: 1 |

## Production — multiplet energies

Route free · partition averaged fold · source `report.md` (`scripts/run.py Si`).

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³P | 3 | −0.667376 | — | −1…1 |
| ¹D | 5 | −0.607963 | +0.0594 | −2…2 |
| ¹S | 1 | −0.518843 | +0.1485 | 0 |

## Lane — the pinned lane

Route pinned · partition the fold column · source `cells.json` (as-run) · pins `test_symm_lane.py`.

<!-- lane:begin -->
| window | CAS | fold | dets | E, Eh | vs declared | term |
|---|---|---|---|---|---|---|
| 3p | (2,3) | {1s,2s,2p,3s} | 3 | −288.4851 | 0 | ³P exact |
<!-- lane:end -->

## Production — the method ladder, CAS(2,3)/{3p}

Route free · partition the treatment column · source `cells.json` (the committed draw).

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | cross-check (mEh) |
|---|---|---|---|---|---|
| ROHF (anchor) | −288.846437 | 0 | — | — | 0 |
| CASCI, averaged | −288.624937 | +221.5 | 2.0000 | 2.5e-15 | 0.1 |
| CASCI, raw | −288.624937 | +221.5 | 2.0000 | 5.6e-04 | ~0 |
| CASCI, fractional | −288.628625 | +217.8 | 2.0000 | 1.3e-15 | ~0 |
| CASSCF, free | −288.846437 | +0.0 | 1.5782 | 1.1e-01 | 0.1 |
| CASSCF, frozen window | −288.628912 | +217.5 | 2.0000 | 6.1e-05 | 0.1 |
| SA-CASSCF, free | −288.846398 | +0.0 | 2.0000 | 7.4e-10 | 0.0 |
| SA-CASSCF, frozen, pur | −288.628792 | +217.6 | 2.0000 | 0.0 | 0.0 |

## Studies (frozen history; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)
<!-- study:end -->
