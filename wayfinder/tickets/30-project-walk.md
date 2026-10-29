# Ticket 30 — The project walk

Review every src and test file the way the 2026-10-28 wave reviewed
system/windows/detbasis/secondquant/angular: file by file, line by
line — purposes, goals, ways, connections, reuse, comments — per the
AGENTS `### The module walk` protocol and the Refactor / Verification
rule lists. Opened 2026-10-28 from the deep-dive session.

## Phases (one session each)

| Phase | Files | Facts homed / decisions expected |
|---|---|---|
| 1. physics core | `operators` → `analysis` | fact 3's PySCF-boundary half; fact 19's certificate; the ladder exception; the CI-ordering bridge consumer |
| 2. substrate engine | `reference` → `core` → `casscf` | facts 10, 15c, 21; 20; 6/7 |
| 3. pipeline | `stages` | the max|dE| three-contributions reading |
| 4. declaration & addressing | `cases` → `cells` → `natural` | facts 24; 23 |
| 5. retirement lens | `scans` → `treatments` → `study` → `drivers` → `pipeline` + scripts | facts 16, 18, 12, 15 — decide retire-or-keep BEFORE homing their facts (see #1 below) |
| 6. test mass | `tests/integration/*` + remaining unit files | pin migration = the retirement deliverable |
| 7. closing | AGENTS final reduction (fact index completion, layout slim); two-axis review + design survey | the navigator+rulebook end state |

DONE already (2026-10-28 wave): `AGENTS` (+ §Refactor/§Verification),
`__init__`, `system`, `windows`, `detbasis`+`secondquant`, `angular`,
tests split unit/integration, the reviewed unit test files. Facts
indexed so far: 1, 2, 4, 17.

Phase 1 DONE (2026-10-29, `456adba`/`abced03`/`87edda6` + the phase
close): `operators` + `analysis` walked (fact 3 collapsed — detbasis
had the physicist formula under a chemist LABEL, fixed; the
regex-parsing kill in shell_occupations; `_physicist_eri` one-home;
`shell_slots` walker reuse; the conftest-basename rule — campaign data
in `integration/refs.py`, unit helpers in `unit/helpers.py`; the
fast-test split with the `test_ecp_wiring` basename precedent; the
fact-24 bridge pinned inside configuration_summary at norb=4). Two
grazed bands widened with observed values (Ti escape gauge 6.40, Sc
partial-swap 0.0736). ENVIRONMENT FINDING recorded in AGENTS Tests:
free-machine suites intermittently crash the C/Si purified rung's
Davidson (NaN, cascade) while isolated reruns are green — a red suite
is never the verdict without isolation reruns.

Phase 2 DONE (2026-10-29, `ce978f2`/`20aa4c6`/`9e3015d`): `reference`
+ `core` + `casscf` walked. Headlines: the R3 one-home fix
(run_frac_occ_scf's inline S-projector metric -> span_projector, the
two occupation hooks now structurally parallel); core's STALE
purified_s_core cites killed (the function died 2026-10-27; two
docstring cites survived as documentation-that-lies) + the literal
\u-escape text fixed; casscf's `_persist` -> `_checkpoint` (watchlist
consumed), CASResult's type lies fixed, facts 6+7 owned by the module
Conventions, the adapters block states the extension rule (a new rung
= a new adapter, never a new kernel path). THE SA CONTRACT FOUND AND
DOCUMENTED (not changed): constrained_casscf returns roots-pass
energies with the SA KERNEL's ci; cells writes that SA ci as ci_stack
under the roots store -- flagged below for phase 4. VERIFICATION
PROTOCOL FINDING: cross-process bitwise A/B is INVALID for the CASCI
kernel (canonicalization rotates degenerate blocks -- the mo/ci O(1)
gauge class; energies wobble 7.1e-15 on buffer-alignment-dependent
last bits) -- the valid protocol is same-process pre-vs-post modules
(bitwise energies + core_energy). Tests: the unit split executed on
touch -- tests/unit/test_reference.py NEW, five windows pins moved to
test_windows.py, test_drivers thinned to the SCF chemistry (19->12)
with owner-direct imports. Facts collapsed: 6, 7, 10, 15c, 20, 21a.
Full suite at the boundary: 188/188 green on a free machine.

DOCSTRING REMEDIATION (2026-10-29, user-directed, same day): the
operators/analysis docstrings were re-judged against the exemplar
passes (system/windows/secondquant/detbasis — the new R1 clause) and
rewritten: opener + definitional-convention bullets, sections on
every non-trivial public function, one sentence per fact, the
Layering TOC dropped, archaeology out. Certified AST-identical
(scripts/check_docstring_only.py — which caught one dropped function
body before it could land; the certificate is part of the method
now, V6).

## Principal differences discovered while planning

1. **Retirement inverts the fact migration**: phases 5–6 may DELETE
   modules; their facts need living homes (surviving modules, tests,
   README) — decide retirement first, migrate after.
2. **The verification environment changed**: external machine load
   (~200+) makes the 13-min suite a 30+ min gamble. Working batching:
   fast tier + A/B as the per-commit gate; the full suite at phase
   boundaries (a DECLINED deviation per Verification rule 6 when
   necessary, stated in the commit message).
3. **The unreviewed mass is the tests, not the src** —
   `test_treatments` (800+ lines), `test_cells` (500+): phase 6 is a
   project in itself and IS the retirement's pin-migration work.

## Protocol

Per file: module + test file together; Refactor walked rule by rule;
Verification walked (A/B where deterministic; degenerate / coincident
/ reference pins actively sought); homed facts collapsed to index
lines; layout row updated; one lockstep commit. Sessions switch at
phase boundaries via a message.txt handoff.

## Naming watchlist (R2; spotted between phases -- each entry is
## consumed at the owning file's phase, renamed WITH its consumers)

- CONSUMED 2026-10-29 (phase 2): `_persist` (casscf.py) -> `_checkpoint`
  (the module's "checkpointed" family; one consumer, renamed with it).
- `_persist_cell` (cells.py) -- the store's language is store/load
  (`load_cell`, `store_name=`): the writer wants to be `_store_cell`,
  the natural pair of `load_cell` -- phase 4.
- (add freely as spotted; an entry dies when consumed)

## Flags for later phases (spotted during phase 2)

- SA store pairing (phase 4, cells): `constrained_casscf` returns
  roots-pass energies/core_energy with the SA KERNEL's ci; cells'
  `_persist_cell` writes that SA ci as `ci_stack` under the
  `<address>_roots` store (whose `energies` are the roots pass's).
  Exact per-root CI exists only in the unread `<address>_roots/ci`
  object-array key. Decide: switch `ci_stack` to the roots CI (changes
  committed SA payloads' config labels at the l-mass margin; needs a
  regen decision) -- documented honestly in constrained_casscf's
  docstring meanwhile.
- SA-weights default layering (phase 4/5): cells imports
  `study.ground_block` for the SA-anchoring rule -- a declaration-layer
  smell (study is a retired-lens module; the rule wants a home beside
  the SA driver or in cases).
- test_drivers.py name (phase 6): after the phase-2 rewires it tests
  the owner modules, not the facade -- rename or split when the test
  mass is reorganized.
