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

## Open physics (runs gated on the user's enable)

- **The U reference-landscape survey**: the gated Newton reference pins
  the ⁵S° basin (a numerically pure L=0 determinant); a
  symmetry-BROKEN determinant may sit LOWER (the generic mean-field
  pattern — every other campaign atom's reference IS broken). Run
  UNGATED, from anisotropic / m-localized starting occupations, both
  Newton and DIIS paths; ask whether any converged determinant lands
  below −50.38185233, and diagnose with ⟨L̂²⟩ + the closure
  certificates. Either answer sharpens the record: a lower broken
  basin re-reads every "below ROHF" anchor in the U tables; none
  found means the symmetric eigenstate-determinant IS the mean-field
  minimum (a noteworthy exception).
- The bare-shell lemma (deterministic, pure numpy — pin it in the
  unit tier when this ticket touches the code): NO single determinant
  of f³d¹ is a rotational scalar — the minimum ⟨L̂²⟩ over all
  all-alpha 4-electron occupations of a bare 1s+1f+1d shell is 6
  (exact 0 only for the half-filled parallel controls p³/d⁵ — the
  top-exterior-power scalars). Hence the ⁵S° purity must live in the
  RELAXED l-mixed orbitals (the wedge of the occupied subspace), not
  in any m-occupation — which is why the mechanism is interesting.

## Progress

- Part 1 DONE (2026-10-31, code-only): U DECLARED in `CAMPAIGN` — the
  raw specs load from `data/` (parsed once; the explicit tag carries
  the identity), the reference spec is the gated-Newton one, the
  window tiers are declared as their exact shell subsets (this basis
  is its own exponent grid, so the record's bands ARE whole-shell
  sets: 3f = 9f10f11f13d, 3f2d = 9f10f11f12d13d = the primary, the
  other four tiers planned), sa_weights = the 17-member ground-term
  M_L set, nroots_primary = 48, run_casscf off (the multi-solution
  class). Code-only verification: the declaration resolves end to
  end (mol builds ECP-honest: nelectron 14, ncore 5; the sector and
  the label<->patterns round trip; the primary grid = 33 cells, the
  planned tiers = 132). NO SCF anywhere.
- Part 2 (next): the reference cache-warm — the classic pipeline
  expects `elements/<label>_rohf.chk`; the record's reference sits at
  `elements/U/cache/U_stuttgart78_dyall_rohf.chk` (one-time copy, a
  committed data move) or a gated fresh SCF; then the term tables
  (gated), the tier CASCI anchors (gated), the landscape survey
  (below).
- run_u* + the two stage adapters retire TOGETHER when the tables
  regenerate through the modern path (their second consumer's death).

## References

- The campaign record: `elements/U/` + the committed `run_u*` output.
- AGENTS facts 22 (the campaign verdicts) and 15 (the open-shell
  rules); the architecture snapshot: `wayfinder/map.md`.
