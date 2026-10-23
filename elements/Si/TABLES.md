# Si — silicon (p²)

cc-pVDZ. Reference **ROHF −288.846437 Eh**. Sectors Ms = 1 / Ms = 0 —
as carbon, with its own numbers. Provenance:
`results/Si_ccpvdz_s2_ladder.json` (T3b), campaign logs. NOTES.md: the
ticket-09 record, gauge draws, the fact-17 baseline.

## T1 — Term content (the contract)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {3p} | Ms=1 / Ms=0 | 3 / 9 | ³P: 3, ¹D: 5, ¹S: 1 |

## T2 — Spectrum

As carbon's pattern (³P < ¹D < ¹S); values in `results/*_report.md`.

## T3b — Treatment ladder at {3p} (cross-check ~95 mEh on the frozen
rows: the 3p|2p-core overlap gives the window-definition offset of
ticket 25; the pur rows close it to 0)

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | cross-check (mEh) |
|---|---|---|---|---|---|
| ROHF (anchor) | −288.846437 | 0 | — | — | 0 |
| CASCI, averaged | −288.624937 | +221.5 | 2.0000 | 2.5e-15 | 94.9 |
| CASCI, raw | −288.624937 | +221.5 | 2.0000 | 5.6e-04 | 95.0 |
| CASCI, fractional | −288.628625 | +217.8 | 2.0000 | 1.3e-15 | 96.7 |
| CASSCF, free | −288.846437 | +0.0 | 1.5782 | 1.1e-01 | 0.1 |
| CASSCF, core-relax | −288.628912 | +217.5 | 2.0000 | 6.1e-05 | 0.1 |
| SA-CASSCF, free | −288.846398 | +0.0 | 2.0000 | 7.4e-10 | 0.0 |
| SA-CASSCF, core-relax, pur | −288.628792 | +217.6 | 2.0000 | 0.0 | 0.0 |

## T5 — Frozen-scan verdicts (definitions in README §7.4)

| radial rotation (m-symmetric) | single-member tilt |
|---|---|
| gains everything (θ* = 0.85, −126.7 mEh) | exactly flat |
