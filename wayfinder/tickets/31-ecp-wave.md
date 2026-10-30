# Ticket 31 — The ECP wave: declare the ECP atom

STATUS: queued (opens on the user's go, after ticket 30 closes).

## Goal

Declare the ECP atom in `cases.CAMPAIGN` — ECP basis + pseudopotential
through the existing `AtomConfig.ecp` pairing — and retire the `run_u*`
scripts into the modern path (the U tables become the declared
campaign's output; `run_u*` stays the frozen record until then).

## Design (pre-verified)

- The seam already exists: PySCF's `get_hcore()` includes `int1e_ecp`
  when `mol.ecp` is set; h1 is built in exactly two places
  (`operators.active_space_integrals`, `operators.full_space_integrals`).
- The ncore-from-Z trap is CLOSED: `system.ncore_of` takes
  `mol.nelectron` (core electrons are removed under ECP — never
  recompute from nuclear charge).
- Optional explicitness: a named `h1_ao(mol, mf)` helper in operators
  would make the seam visible — decide during the wave.
- The ECP machinery is validated (pairing + p² and TM campaigns:
  ccECP + ccecpccpvdz; `tests/unit/test_ecp_wiring.py` +
  `tests/integration/test_ecp.py`).

## Retirement riding this wave

- `run_u*` scripts — the tables regenerate through the modern path.
- `casscf.run_casscf` + `run_sa_casscf` — pass-throughs carrying the
  stage name; the `name=` becomes explicit at the cases call site.
  `run_casci` KEEPS (real content: the CASCI class, `nroots`,
  `fci_opts`).

## References

- The campaign record: `elements/U/` + the committed `run_u*` output.
- AGENTS facts 22 (the campaign verdicts) and 15 (the open-shell
  rules); the architecture snapshot: `wayfinder/map.md`.
