# Ticket 44 — The PySCF-mrpt adapter (the new `perturb.py`)

## Goal

The analyzation module, PySCF-backed: one deep seam consuming the
substrate `(mol, mf, mo_full, ncore, ncas, nelecas)` + explicit
`roots`, producing per-root per-Koopmans-class rows. The engine is
`pyscf.mrpt.NEVPT` (strongly contracted; 2.13.1 native). Internally:
the multi-root CASCI rebuild on the stored MO (the readout-twin
recipe) with the BIT-EXACT stored-roots assertion; the SA guard
exercised; the 8 class functions' `(norm, energy)` extraction; the
per-class sum identity. A variant registry keeps CASPT2 (env upgrade)
and further variants row-compatible later. Spec:
`wayfinder/specs/pt2-perturbation.md`.

## Status: READY (Step 0 done 2026-11-19 — anchors recorded:
O2 −0.169785475 vs the example's −0.169785157, 3.2e-7; C2 SA pattern
green, 0.3 s)

## Inputs

- The stored cells (`elements/{C,Si}/cells.json` + chks: `mo_coeff`,
  `ci_stack`, `energies`, `ncas` per address);
- `casscf.run_irrep_casci` (the readout recipe's home);
- `pyscf.mrpt` (`NEVPT`, `nevpt2.Sr/Si/Sijrs/Sijr/Srsi/Srs/Sij/Sir`).

## Verification

- Unit tier (frozen carbon fixture, NO SCF; synthetic MC objects for
  the pure parts): per-class sum == `e_corr`; class extraction vs a
  pinned scalar run; the SA-solver refusal raised; row shaping.
- Anchors: O2 + H14 (the example values, ~1e-6 Eh band).
- Bit-exactness: rebuilt CASCI roots == stored `ci_stack` (the
  ticket-43 precedent), asserted before any PT2.
- The per-commit gate: fast tier (258) + `tests/unit/test_perturb.py`.
