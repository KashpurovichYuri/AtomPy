# AtomTerm — the effort map

The navigator: what is live, what is open, where the worked examples
live. History is git's job — resolved tickets are deleted, survey dumps
are not kept (open findings move to Issues below before their section
dies).

## Tickets (live)

- [29 — gated measurement batches](tickets/29-ti-first.md): the T2 Ti
  cells + the M store migration; runs only on the user's enable.
- [31 — the ECP wave](tickets/31-ecp-wave.md): declare the ECP atom in
  `CAMPAIGN`; opens on the user's go.

## Issues (open dispositions + watchlist)

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
    assembly       operators (both PySCF boundaries)
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
