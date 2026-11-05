# Ticket 31 — The ECP wave: declare the ECP atom

STATUS: queued (opens after ticket 32's close; the door plan SETTLED
with the user 2026-11-05).

## Goal

Declare the ECP atom in `cases.CAMPAIGN` — ECP basis + pseudopotential
through the existing `AtomConfig.ecp` pairing — and retire the `run_u*`
scripts into the modern path (the U tables become the declared
campaign's output; `run_u*` stays the frozen record until then).

## The U14e door plan (settled with the user, 2026-11-05)

1. THE ELEMENT DIR IS `elements/U14e/`: rename the CAMPAIGN key
   `U` -> `U14e` (regen derives the element dir from the key -- that
   is the whole trick); the `tag` stays `U_stuttgart78_dyall` (the
   basis/psp label on chk names).  `elements/U/` STAYS FROZEN (the
   user's call; archive-or-keep revisited later).
2. CACHE-WARM: the record's reference chk
   (`elements/U/cache/U_stuttgart78_dyall_rohf.chk`) copies ONE TIME
   to `elements/U14e/rungs/rohf.chk` (the trusted-cache convention --
   cached chks never re-gate), or a gated fresh SCF.
3. THE UHF CENSUS: five arms at Dyall size (~30 min budget) ->
   declare `occ_e_ref_uhf` -> `scripts/uhf_anchor.py` -> the anchor
   row.  ONE FLAG: the census arms run without the Newton accelerator
   U's ROHF class needed (fact 22) -- if an arm diverges, the census
   learns `ref.conv_accelerator` (extend only on evidence).
4. THE FIRST REGEN IS ANCHOR-FIRST (~20 cells, NOT the full grid):
   declare an `anchor-first` Scenario for U14e (Ti's pattern) --
   the primary window's CI cells + the raw detector run; the four
   tier windows run their CASCI anchors only (~2/tier); SA and the
   optimizer variants stay VOID at the key (the multi-solution
   class; `run_casscf=False`'s spirit made explicit), one `--full`
   away.  U has NO scenarios today -- an unfiltered regen would
   charge into the 33-cell grid + all 132 tier cells.
5. THE RETIREMENT (G5): `run_u*` + `run_casscf`/`run_sa_casscf` die
   with the first modern U tables; `stage_facts`' `c_act` absorbs
   here (the one caller-change moment).
6. Ticket 29's T2/M batches stay user-gated.  The .chk restart
   anatomy needs NO work (every rung's chkfile IS the restart
   artifact: orbitals + CI vectors persist per address).

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
- Part 2 EXECUTING (2026-11-05, the door plan; commits 408387b,
  bca1057):
  1. THE KEY RENAME: `CAMPAIGN['U14e']` (symbol stays `U`; tag stays
     `U_stuttgart78_dyall`); `elements/U14e/` born, `elements/U/`
     frozen. CaseSpec gained `key` (the element-dir contract, one
     lookup: `campaign_spec`); CellKey admits campaign keys.
  2. CACHE-WARM: the record's chk copied to
     `elements/U14e/rungs/rohf.chk`; the trusted load reads
     E = -50.38185233 (the record's value).
  3. THE CENSUS (five arms, Newton class — the door plan's flag
     SPENT on evidence: all five plain-DIIS arms diverge,
     E ~ +2.4e3..+3.0e3, so the census learned
     `ref.conv_accelerator`): FOUR arms identical to 1e-6 —
     E = -50.38881343, conv, <S^2> = 6.2289 (exact quartet 6;
     mildly broken), sigma_min = 0.8835; huckel non-converged (not a
     second basin). SINGLE BASIN, and 6.96 mEh BELOW the ROHF — the
     open physics question ANSWERED for the UHF class: the
     mean-field minimum is the mildly-broken determinant, not the
     ⁵S° eigenstate-determinant. `occ_e_ref_uhf = -50.38881343`
     declared. (The ROHF-class broken-landscape legs stay gated.)
  4. THE FRAC REFERENCE HAS NO CONVERGED BASIN on the ECP atom
     (measured): DIIS diverges (+137); Newton collapses (-12.74 vs
     the aufbau RHF -49.04); Newton-from-ROHF does not leave the
     seed. The `frac-reference-unconverged` scenario voids the five
     frac-fold cells; MOM-anchored selection is the recorded fix
     (map Issues). The regen head builds the frac reference only
     when a non-void frac fold runs; both head guards now unlink a
     diverged reference's chk (the trusted-cache poison class).
  5. THE ANCHOR-FIRST REGEN: the gate voids the whole optimizer
     class (145 cells, all windows incl. the NO forms) + the frac
     class → 15 CASCI anchors run. IN FLIGHT (detached; priced on
     the real primary cell: 272 s, kernel-only payload class).
- Part 3 (next): uhf_anchor.py U14e → build_methods_tables.py U14e
  (the first modern U14e TABLE.md) → THE RETIREMENT (G5): `run_u*` +
  `run_casscf`/`run_sa_casscf` die, `stage_facts`' `c_act` absorbs.

## References

- The campaign record: `elements/U/` + the committed `run_u*` output.
- AGENTS facts 22 (the campaign verdicts) and 15 (the open-shell
  rules); the architecture snapshot: `wayfinder/map.md`.
