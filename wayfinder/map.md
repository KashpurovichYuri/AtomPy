# AtomTerm — the effort map

The navigator: what is live, what is open, where the worked examples
live. History is git's job — resolved tickets are deleted, survey dumps
are not kept (open findings move to Issues below before their section
dies).

## Tickets

Queue decision (2026-11-01, user): ticket 31 (U) runs AFTER the
ticket-32 close -- the C/Si/Sc/Ti tables land in their redesigned
contract first, and U inherits the settled table grammar and
declaration path as a fresh wave. (live)

- [29 — gated measurement batches](tickets/29-ti-first.md): the T2 Ti
  cells + the M store migration; runs only on the user's enable.
- [31 — the ECP wave](tickets/31-ecp-wave.md): declare the ECP atom in
  `CAMPAIGN`; opens on the user's go.
- [32 — the symmetry lane](tickets/32-symmetry-lane.md): D2h-pinned CI
  as the second route (compare, never merge) — step 8 deep in flight:
  the tables grammar + the CI window survey landed, the ROUTE AXIS
  executed on Sc (the twins run; the routes agree at SA on m-complete
  windows); next: sparse meters, saturation legs, the Ti/C/Si passes,
  then the wave close.
- [33 — the symmetry-breaking anatomy deliverable](tickets/33-symmetry-anatomy.md):
  design-first — the frozen-landscape tables that replace the retired
  scan pointers.

## Issues (open dispositions + watchlist)

- The f41bb51..HEAD review (2026-11-04, both axes; fixes in the
  round's commit, ticket 32's closing entry is the record).
  DESIGN SURVEY (the delta seams; the full src/ survey stays at the
  wave close): reference.py got DEEPER (one cache/gate skeleton +
  one acceleration dispatch behind three runners; the object
  construction `_uhf_object` is now a REAL seam -- two adapters:
  the production runner and the census arms; the declined cells-axis
  seam is documented at the module contract, where a seam refusal
  belongs); the UHF gate's declaration rides the row's
  ReferenceSpec (the declaration-owns-intent rule held).  The one
  real shallowness finding is item (1) below.  The scripts'
  vocabulary module tightened (comm_text deleted).
  RECORDED, not fixed: (1) the kernel/lane residual block is written
  twice (`kernel_cell_payload` vs `run_lane_cell` -- the apply/residual
  loop and its imports) -- the one-home candidate for the wave-close
  refactoring round; (2) multi-invariant pins (`test_pyscf_apply_is_
  the_slater_condon_action` x3, `test_method_ladder_one_table`) --
  V1 splits at the wave close, the pins are load-bearing; (3) the new
  residual bands (`< 1e-8`) carry no observed fresh-value comment
  (V7) -- capture at the next natural fixture run; (4) the "cycle
  cap" note keys on the DECLARED 200 while the payload carries only
  the count -- the proper fix stores max_cycle in the payload (a
  schema bump; rides the next store migration); (5) `cas_label`'s
  NO branch is a hypothetical seam (NO starts never reach labels) --
  delete at the wave close unless a consumer appears.

- The ladder-wave review round (2026-11-01, `701b31c..c02d0c6`, both
  axes): applied — the module-docstring lockstep, the shared
  `_substrate` builder, the V1 test split, the `_FOLD_*` renames, the
  `LadderLeg` NamedTuple.  Declined: splitting the ladder sessions
  out of ticket 32 (the ticket IS the provenance mechanism, per the
  Conventions); the map queue-decision placement nit (live prose
  until the wave closes).  Recorded partials (spec axis): the ladder
  TABLES omit a sector column (implicit in dets — add when the
  emitter generates the tables); the ours-vs-solver anchor runs on
  one leg per construction, not per leg (the uniform-construction
  argument — revisit if a leg's meters ever disagree).
- Pinned cells' window intake (the symmetry lane's cross-route
  question, 2026-10-31): the natural pinned-cell windows are
  D2h-BORN sources only (raw shells, D2h-run MOs/NOs) — free-run
  windows never cross into pinned cells; the grammar states the
  restriction at the declaration wave.  The trap it avoids:
  `label_orb_symm` labels m-mixed columns by dominance, so a
  free-run window fed to a pinned solver would appear to work and
  produce a silently wrong solve — any future intake needs an
  explicit irrep-purity check at the seam.  The laundering routes
  (NO-density α, radial-transplant β) stay diagnostic experiments,
  not grammar values.
- `run_casscf`/`run_sa_casscf` are pass-throughs carrying the stage
  name — retire at the ECP wave together with `run_u*` (their second
  consumer); the stage `name=` lives at the call site (2026-10-30).
- `stage_facts`'s `c_act` parameter: every production caller passes the
  same derived slice — deepening candidate, deferred (2026-10-30).
- The `verbose` duplication (a `run_cell` kwarg AND a budgets key) —
  the 2026-10-19 finding; confirm done-or-dead at the next cells touch.
- Ti escape comm band (test_treatments, the sa_free band): quiet
  through recent suites; capture-not-widen stands (no observed value
  recorded).
- The D2h twin's determinism certification is procedural (double-run
  sha256 in the session record), encoded nowhere in-repo; a generator
  `--verify` mode is the candidate — declined until a second regen
  cycle (2026-10-31).
- The D2h mol construction has two homes by design (the generator's
  draw, the fixture's load); a shared builder would serve one consumer
  each — revisit at a third frozen reference (2026-10-31).

## Examples (worked records in `elements/`)

- `elements/C/`, `elements/Si/` — the p² campaign: the methods grid,
  scans, multiplet tables.
- `elements/Sc/`, `elements/Ti/` — the transition-metal campaign:
  rungs, the radial axis, the d² lattice, term tables.
- `elements/U/` — the ECP-atom record (`run_u*` until declared).

## Architecture (edge-verified)

    pure math      system · detbasis · secondquant · natural
    substrate      windows (the hub)
    1-particle     angular (uses windows' shell_slots)
    assembly       operators (both PySCF boundaries) · irreps (the
                   D2h dictionary; the symm boundary)
    engines        reference · core · casscf
    analysis       analysis · stages · scans
    declaration    cases -> cells

`natural`'s mathematics is pure at module level (numpy + secondquant
only — the windows→natural edge retired with the window-kind
registry, 2026-10-30); its one driver (`no_group_windows`) lazily
imports the assembly layers. No cycles, no upward edges. Division of
labor over windows: declaration NAMES (labels), windows BUILDS (the
kinds + substrate builders), cells BINDS (the start axis picks the
builder).
