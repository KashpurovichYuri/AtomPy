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
