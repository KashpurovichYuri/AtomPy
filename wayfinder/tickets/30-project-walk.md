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

Phase 3 DONE (2026-10-29): `stages` walked. Headlines: the two ruler
readings HOMED (module Conventions + stage_facts Notes) -- the
reconstruction is the ANALYSIS ruler (meters/labels computed on it;
the tables' E columns carry the run's own totals; never compared),
and max|dE| decomposes into core-field anisotropy / window definition
/ bookkeeping, the last exactly 0 (the fold convention is
machine-exact vs get_h1eff, 1e-15), the window-definition term
dominating on core-overlapping shells (Ti ~85 mEh) -- fact 14's tail
+ fact 25's reading, content standing alone with no fact tags. The
SA-ci CONSUMPTION documented: cas.ci labels the stage, and on SA
runs it is the SA kernel's CI (the consumer-side half of the
constrained_casscf contract wart; both label consumers bridge the
ordering internally). StageFacts: fields regrouped (the defaulted
fields were interleaved WITH the readout properties; the "three
analysis fields" comment was a count lie over four), payload keys
declared committed data (renaming one is a payload break), the
window property's 1e-6 rounding stated (regen-diff stability), the
dead `field` import dropped. stage_facts: proper
Parameters/Returns/Notes; the ticket-25 NOTE moved out of the c_act
description into Notes where the reading lives; the converged read
commented at point of use. Tests: imports owner-direct (were through
the drivers facade), V1 docstrings, DE-SLOWED whole-file to fast
(frozen reference + norb=3 full diagonalization, 0.6 s wall,
deterministic) -- fast tier 43 -> 45. AGENTS: the stages row carries
the readings; the Tests count fixed 188 -> 189 (the phase-2 split's
edit that never landed) and fast re-worded as the machine-light
selection. Suite at the boundary: 188 green + 1 Sc complete-window
session draw, isolation-green (the recorded class, no band change).

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

## Naming review round 3 + the phase-2 review wave (2026-10-29)

The user's review lens, stated for the remaining phases: GENERALIZE,
UNIFY, and make similar approaches LOOK similar. Applied:

- MASKS: the `mask_` family, names stating the constraint --
  `mask_frozen_window`, `mask_no_core_window_mix` ("exchange" retired:
  collides with the exchange integral), `windows.mask_radial_rotations`
  ("same-shell" rejected: shell = l here and same-l rotations mix
  members). Serialized axis tokens (frozen/virtmix/radial) stay -- they
  are committed data inside cell addresses; the mapping is stated once
  at the cells layer. Round 2's three-way dictionary is SUPERSEDED.
- `run_frac_occ_scf` -> `run_frac_occ` (name-sibling of run_rohf; the
  formalism is closed-shell RHF, NOT ROHF -- run_rohf_frac_occ was
  considered and rejected). MOM caveat recorded in its Notes (the
  across-cycles anchoring weakness; the fix if a future atom
  oscillates), no code added.
- run_rohf's split acceleration dispatch unified: `solve` owns all
  three policies ("newton"/"diis+newton"/"diis"), `fresh_mf` is pure
  construction.

## Deep reuse — what new methods will need from casscf (recorded per
## the phase-2 agreement; NOT designed)

- The ONE driver: `constrained_casscf(mo_full=, mask=, weights=,
  mc_opts=, fci_opts=)` covers every start x constraint x averaging
  combination. A new method (spin-orbit, g-factors, transition
  moments, a new rung) arrives as a new `mask_` (data,
  uniq_var_indices-shaped) or a thin adapter -- never a kernel fork.
- Per-run knobs without subclasses: `mc_opts`/`fci_opts` dicts set
  attributes pre-kernel (budgets, exact-diagonalization forcing).
- The results seam: `CASResult` (results only) + the chkfile store
  grammar `atomterm/<name>/`; the SA TWO-PASS pattern (ensemble kernel,
  then the exact roots pass with canonicalization off) is the template
  for any method needing exact per-root quantities in an optimized
  window.
- The mutation discipline (fact 6): read `converged` off the live mc,
  nothing else; pass `mo` explicitly for integrals.
- DESIGN VERDICT (the deletion test): the adapter layer except
  `run_sa_casscf` (its one surviving consumer: cases) is a middle-man
  -- delete it and no complexity reappears (every caller writes one
  constrained_casscf line, as cells already does). RETIRE with
  treatments/study in phase 5; do not migrate (their store keys
  casscf_relax/casscf_virtmix are rungs-era data, replaced by cell
  addresses).

## Naming watchlist (R2; spotted between phases -- each entry is
## consumed at the owning file's phase, renamed WITH its consumers)

- CONSUMED 2026-10-29 (phase 2): `_persist` (casscf.py) -> `_checkpoint`
  (the module's "checkpointed" family; one consumer, renamed with it).
- `_persist_cell` (cells.py) -- the store's language is store/load
  (`load_cell`, `store_name=`): the writer wants to be `_store_cell`,
  the natural pair of `load_cell` -- phase 4.
- (add freely as spotted; an entry dies when consumed)

## Flags for later phases (spotted during phases 2-3)

- The ordering bridge's HOME (phase 4, cells): `pyscf_ci_permutation`
  + `ci_in_detbasis_order` live in cells, but analysis (a LOW layer)
  consumes them via lazy imports against the cells->stages->analysis
  cycle. detbasis owns the ordering statement ("ordering != PySCF FCI
  strings") and sits BELOW both consumers -- consider homing the
  bridge there when cells is walked (consumers: analysis x4, study x2,
  test_cells x7).

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

## Naming review round 2 (2026-10-29, user-directed, phase 2)

- Mask polarity verdict (do not relitigate in phases 4/5): the two
  geometry masks deliberately carry OPPOSITE grammar --
  `core_virtual_only_mask` names its single ALLOWED pair,
  `no_core_window_exchange_mask` its single DENIED pair -- because the
  allowed sets have complementary shapes; one polarity would make one
  name wrong. Symmetry lives in the docstrings (allowed + denied +
  axis token each) and the three-layer DICTIONARY lives in casscf's
  module Conventions (geometry / purpose / axis).
- `CASResult.kind` DELETED (never read anywhere; `rung.kind` in
  treatments is a different object) -- the container carries results
  only; stage identity flows via stage_facts(tag=)/store keys.
- core.py `n_core` -> `ncore` (the spelling every other module uses;
  no keyword callers existed).
- reference.py re-audited clean: `occ_e_ref`, `conv_accelerator`,
  `run_frac_occ_scf` match the established tokens (frac/pur axis
  values, prose in facts/tickets); `_forced_rohf_occ` clear.
