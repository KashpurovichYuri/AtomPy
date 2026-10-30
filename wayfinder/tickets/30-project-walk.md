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

Phase 4 DONE (2026-10-29, `549f202`/`197d8d2`/`fd48ad9`/`f2aafa0` +
the close): `cases` -> `cells` -> `natural` walked (modules + test
files together), in four commits: (1) cases -- `want` DELETED (the
R5 fake degree: 'casscf'/'sa-casscf' read nowhere; a case without
term tables now unrepresentable, the arm gates on spec.study),
public-seam asserts -> ValueError, CaseSpec -> Attributes, test_cases
moved to the UNIT tier (all SCF-free) with a parity pin replacing
the dead demands-test, test_atoms owner-direct; (2) the ONE-HOME
MOVES (the phase-3 flag + the layering flag, resolved structurally):
the CI-ordering bridge cells -> OPERATORS (detbasis died as the
candidate -- it is deliberately numpy-only and the bridge needs
cistring; operators owns the other PySCF boundary + sits below every
consumer; analysis's lazy anti-cycle imports became top-level),
ground_block study -> OPERATORS ops-first (absorbing natural's inline
duplicate found by the walk; vec returned so no consumer repeats the
eigh; kills the cells->study smell), sector_dimension homed in
detbasis (the comb formula's FOURTH copy found in natural; cases +
study + cells copies deleted); (3) cells -- `_persist_cell` ->
`_store_cell` (watchlist consumed), CellKey/CellData -> Attributes,
mol documented in run_cell for the first time, the last src 'rung'
retired; (4) natural -- the ONLY Doxygen docstring in src rewritten
to NumPy, the stale registry header deleted, the cost-stop asserts ->
ValueError, fact tags dissolved. Facts 23/24 collapsed to index
lines. FAST TIER 45 -> 48.

DECIDED-KEEP (SA ci_stack pairing, the phase-2 wart): the exact
per-root CI is unreachable from `cas` without changing
constrained_casscf's contract, no V6 certificate exists for optimizer
paths, and no consumer asks for it -- phase 2's documented stance
stands (constrained_casscf + _store_cell docstrings carry it).

PHASE-5 INPUT (from the cases walk): the treatments arm has NO
declaring case anywhere (src+scripts) -- `_treatments_arm` +
CaseSpec.study are the only cases->study bridge; decide with study's
retire-or-keep.

PHASE-5 CONSIDERATION (from the phase-4 close, user-directed):
the campaign carries TWO run vocabularies -- the classic staged
pipeline (cases.run_atom, stage tags like `casci_ms1`, the T1/T2 term
tables) and the methods grid (cells, CellKey addresses, the manifest
tables). They overlap in their CASCI stages but serve different
deliverables; whether run_atom's stages become cells (one vocabulary)
or stay a separate deliverable class is a phase-5/7 design question --
raise it when the retirement trio is on the table, do not decide now.
Extension story recorded in both module docstrings: new ATOM = a
CAMPAIGN row (cases); new AXIS VALUE = vocabulary + one resolution
branch through its owner (cells); new METHOD = a mask_ or a thin
adapter over constrained_casscf (the deep-reuse paragraph above).

Phase 5 DONE (2026-10-29, `1c2ee2f`/`40769fa`/`ad28d3a`/`4565e93` +
the close): the retirement lens. DECISIONS (retire-or-keep BEFORE
fact-homing, per the plan):

- scans KEEP + walked: the optimizer-free deterministic anatomy
  (fact 18) has no overlap with cells (which runs optimizers);
  fact tags dissolved, ScanPoint Attributes, scan_window's assert
  -> ValueError, test_scans owner-direct.
- treatments RETIRED: zero src consumers; every rung maps onto a
  CellKey address. THE PIN MIGRATION (user-directed actuality
  review): 19 of 20 pins are ACTUAL (they pin measured campaign
  physics, not the rungs machinery) and now run through run_cell --
  LADDER_ADDRESSES in test_treatments.py is the rungs->cells
  Rosetta stone; one pin DIED with its subject (_report's renderer
  flag; the payload contract stays pinned in test_stages). One
  grazed band widened with the observed value (Ti casci comm 1.1e-7,
  isolation-green).
- study RETIRED with the arm (CaseSpec.study deleted, R5): the arm
  had NO declaring case; run_cell owns the budgets; the committed
  NO-study inserts in C/Si/Sc/Ti/U TABLES.md are past one-shot
  declarations kept as frozen data. The unexercised quarantine
  dies -- if a big-window cells campaign ever needs failure
  isolation, the try/except belongs at the runner (regen_cells),
  not in the cell.
- 3 casscf adapters RETIRED (core_relax/virtual_mix/pyscf_start);
  run_casci/run_casscf/run_sa_casscf survive (cases + cells -- the
  classic pipeline's stage names; the phase-2 verdict had missed
  run_casscf's cases consumer).
- drivers + pipeline RETIRED: 28 files rewired owner-direct; two
  test files renamed honestly (test_drivers -> test_substrate_
  chemistry, test_pipeline -> test_run_carbon -- the phase-6 rename
  flag consumed early because a test file must not name a dead
  module); run_u_sa's relax stage calls constrained_casscf directly.
- sa_decomposition.py RETIRED (user: review actuality): its
  attribution (window vs core-fold) rides every cells row as the
  span-certificate + anisotropy columns; artifacts stay committed.

DECIDED (the two-run-vocabularies question, raised at phase-4
close): KEEP BOTH. The classic staged pipeline (stage tags, the
T1/T2 term tables: all-Ms sectors, full spectra, human reports)
and the cells grid (addresses: method comparisons on one window)
are DIFFERENT DELIVERABLES at different speeds -- term tables are a
stable contract (new atom = new row), the grid is under active
development (ECP wave will add axes); unifying couples them so
every cells change risks the committed tables' reproducibility.
Cost of both: ~50 shared lines and two store vocabularies --
bounded, stated in the division-of-labor docstrings. README's
pipeline diagram de-lied (it still showed the phase-4-deleted
`want` dispatch).

SUITE AT THE BOUNDARY: 184 collected, 183 green + 1 Ti escape-draw
flake, isolation-green on rerun -- the recorded rotating-flake
class, no band change (phase-3 precedent). Fast tier 115 (measured;
the handoff chain's "48" was a unit-tier figure). AGENTS: layout
rows for the four dead modules deleted, facts 9/11/12/19/20
re-pointed (drivers->owners, treatments->the migrated pins),
reference.py's run_rungs cite -> run_case.

POST-PHASE-5 ROUND (2026-10-29, user-directed, `0909b64`): (1) the
method axis takes the FULL name -- 'sa' -> 'sa-casscf' everywhere
(vocabulary, addresses, payload kind, the classic stage tag, display
maps), committed data migrated in lockstep (4 cells.json stores,
C's report, 12 chk caches via h5py subtree copy, round-trip
certified); identifiers keep underscores. The rename EXPOSED one
behavioral miss (a budgets comparison still matching 'sa') -- and
that miss found a REAL TRAP, now run_cell's contract: the
measured-block sa_weights default presupposes a DEGENERATE start
window; on a defective (fact-17-split) window ground_block anchors
ONE state -> single-root SA -> the TM collapse basin (Sc measured:
nroots=1, L2=0.047). (2) cells' comment debt brought to the
accepted form (arrival-history clauses content-ified; the
archaeology grep across src is clean). GRAZE LOG: the Ti escape
comm band (test_treatments line ~817) failed in 2 of 4 full
sessions on 2026-10-29, isolation-green 4/4, value uncaptured --
no widening without the observed value; capture it if it recurs.

ARCHITECTURE SNAPSHOT (2026-10-29, edge-verified by grep over the
top-level imports -- the phase-7 design-survey input; corrects the
conversation diagram that omitted windows/system and misplace
natural at declaration level):

    pure math      system · detbasis · secondquant · natural
    substrate      windows (the hub: angular, operators, core,
                   reference, stages, analysis, cases all import it)
    1-particle     angular (uses windows' shell_slots)
    assembly       operators (angular · detbasis · secondquant ·
                   windows; both PySCF boundaries)
    engines        reference · core · casscf
    analysis       analysis · stages · scans
    declaration    cases -> cells

windows -> natural is correctly directed (the kind registry's 'no'
kind needs the NO mathematics; natural is pure math, imports only
secondquant). Division of labor over windows: declaration NAMES
(labels), windows BUILDS (five kinds + substrate builders), cells
BINDS (the start axis picks the builder). No cycles, no upward
edges.

Phase 6 DONE (2026-10-30, bd2dbb6/6034c08/b6bc7f2/b1ca436): the
scripts surface + the test mass, one session (with the two review
skills pulled in as directed -- the phase-7 review material executed
early).

SCRIPTS (the agreed classification executed):
- run_c.py -> run.py [atoms]: the per-atom launcher over CAMPAIGN
  (no args = every declared atom in CAMPAIGN order; fail-fast on
  undeclared names). run_c's vestigial [workdir] arg dropped (its
  'results' default contradicted the README contract).
- RETIRED: run_tm_lattice / run_tm_methb / run_tm_no / multiplet_split
  (actuality review: findings institutionalized in facts 22/23/25,
  artifacts committed, ticket-29's gated batches run through
  regen_cells -- none is the gated tool). Cites now read
  'retired -- git history' (the sa_decomposition pattern).
- run_u* untouched (the U fork = ECP wave: declare U in CAMPAIGN).
- Live tools: the sys.path bootstrap idiom UNIFIED on run.py's robust
  form -- radial_scan had NO bootstrap (its usage line failed); the
  freeze script's PYTHONPATH=src note replaced by the same bootstrap.

TESTS (the never-walked files, walked):
- test_treatments: LADDER_CASES archaeology trimmed; the duplicated
  escape-comment collapsed; the Ti sa_pyscf comment's Sc copy-paste
  ('restored 5.99+ class' -- Ti restores 12) fixed; the divergence
  guard MOVED to test_stages (its true home; fast tier 115 -> 116).
- test_atoms: the stale test_pipeline cite -> test_treatments.
- test_run_carbon: function renamed with the file (phase-5 flag
  fully consumed); the SA-breaks comment de-lied (single-root CASSCF
  breaks; SA RESTORES -- the C pins say so).
- test_ecp: redundant local imports -> module-level operators.use;
  unused build_mol dropped. conftest/ecp_wiring: import hygiene.
- DECLINED (argued): carbon-fixture cfg removal (family symmetry with
  sc_rohf/ti_rohf whose cfg IS consumed); test_cells fast-split (its
  docstring states the 1:1 mixed-tier policy -- AGENTS sentence
  rewritten accordingly, FLAGGED to the user as a working-agreement
  change).

THE BOUNDARY SUITE CAUGHT A REAL ONE: test_atoms'
results['sa_casscf'] -- the 0909b64 rename's straggler (slow-only
lookup, invisible to -m fast; 2 failed 182 passed, both the one
KeyError, green in isolation after the fix). Lesson recorded in map:
a rename round touching stage tags owes a full suite run.

THE SKILLS (phase-7 material, executed early per the user):
- code-review (fixed point 30ab209): the R6 lockstep misses fixed
  (diag_fold/basin_map dangling cites in README + stages.py + C
  TABLES, map item-5 supersession, README counts 130/150 -> 184);
  the spec axis challenged the test_cells decline -- the stated
  policy stands (module-level fast WOULD leak; module-slow +
  function-fast double-selects into both -m tiers).
- codebase-design over src/: save_report retired (zero consumers);
  run_case's missing return fixed; the SELECTORS/window_* registry
  found to be a HYPOTHETICAL SEAM (zero production consumers -- the
  dispatch reappeared in run_cell's branches) -- decision recorded
  (retire vs consummate) for phase 7 / the ECP wave; run_casscf and
  stage_facts' c_ast noted borderline/deferred. map.md carries the
  full 2026-10-30 survey section.

REMAINING for phase 7: AGENTS final reduction, the SELECTORS
decision, the recorded candidates, the closing two-axis pass.

POST-PHASE-6 ROUND (2026-10-30, user-directed, bdb4b76/273ad46/
265cb27/98d4c7b -- the plan agreed in conversation):
- freeze_carbon_reference -> generate_fixture_carbon (the doing-name;
  cites lockstep). Lint DECLINED (the noqa comments stay as honest
  annotations; no toolchain).
- regen_cells GRANULARITY: key-identity filters (--method/--window/
  --mask/--start/--substrate/--fold/--sector, '-' = maskless) +
  --group grid|spine|diagnostics|planned; filtered runs MERGE into
  cells.json (untouched cells keep payloads), the full run replaces;
  empty-filter queues skip before any SCF; module-level body ->
  main(); the dead rungs.json cleanup block deleted. Certified: the
  selection matrix by direct call + the end-to-end C --method
  sa-casscf run (exactly the 14 cells refreshed, 19 preserved --
  store reverted: the refresh carries the electrons-field payload
  evolution that belongs to a deliberate full regen, ticket 29 M).
- README: 'The variety axes' section after the pipeline diagram (the
  axes table, two runners/two deliverables, filter examples, merge
  semantics, undeclared-exploration paths).
- Scripts docstring pass to the walked-src standard (what+why ->
  Usage -> Notes; function one-liners); AST-certified on 4 of 5
  files, renderers byte-idempotent on Ti; the one code change:
  generate_fixture's shadowing inner import deleted.
- Sequencing confirmed: phase 7 = small close; U-declaration opens
  the ECP wave (run_u* retire then).

THE GENERALIZATION WAVE (2026-10-30, user-directed, c2f7fd4/bdd970a/
20a4400/0264796/6f71775 + round B): the user's lens -- 'repeating
information or providing it non-generally might mean some issues' --
executed as three PROMOTIONS + two strict docstring rounds.

POLICY AMENDMENT (AGENTS R1): definitional slots carry the RULE
stated generally -- no atom names, no measured campaign values, no
fact/ticket tags; instances live in Examples/Notes sections; the
measured record lives in AGENTS facts / tests (V7) / README / TABLES.

A0a LATTICE PROMOTION: CaseSpec.lattice_windows + cells.lattice_keys;
regen's --group lattice + the 'if name == Ti' attach die;
build_lattice derives spans/sector/CAS labels (D2_SERIES deleted).
The A/B pin CAUGHT an order delta on first try (per-span triplets vs
the original ao/mo-then-NO blocks) -- byte-equality demanded the
original layout; the pin is now test_cells' permanent derivation test.

A0b SCAN PROMOTION: scans.ScanSpec + CaseSpec.scan; radial_scan
becomes the thin CAMPAIGN runner (SPECS dies). Certified: C's scan.md
BYTE-IDENTICAL from the campaign cache; Si fills a GAP (no artifact
had ever been committed; values = the README-cited theta*=0.85,
-126.7 mEh); the OUTPUT PATH BUG fixed en route (the old script wrote
untracked dash-named files, never the committed element-dir artifact).

A0c VERDICT ELIMINATION (supersedes ticket 29 T1b, noted there):
build_methods' AUDIT dict deleted, the verdict column gone, the
pending gate payload-only; the anchor line COMPUTED from the store --
which EXPOSED A STALE LITERAL (the committed caption's eV value was
wrong by 0.012 eV); README's 'Verdicts' paragraph became the READING
KEY (every class + its row-level discriminator); the hand summary
table stays as the frozen worked example.

ROUND A (scripts strict): build_methods/build_lattice conventions
bullets deleted (verified covered by README line-by-line first); the
anecdote, arrival clauses, fact tags out; render() documented.
AST-certified; both renderers byte-stable on Ti.

ROUND B (src strict): 9 edits across analysis/core/reference/casscf/
stages/windows -- fact tags dissolved (the phase-1 pass had missed
analysis's three), measured values re-homed (Ti's +123/+11.2 mEh,
U's 38 mEh and 14.5 Eh enumerations, the ~85 mEh example; kept as
classes where the caller needs the scale), cells.py's module
docstring sectioned (structure-only). All AST-IDENTICAL.

KEPT AS RULES (declined dissolution): scans' 'exactly 0 Eh' (the
invariant, not a measurement); scales stated as classes (cEh-scale
wobble, m-scalar by construction).

WAVE CLOSE (2026-10-30, the user's cells-opener review): the opener's
module enumeration was WRONG AS A MAP (stages = the analysis seam,
not an axis owner) and its example address was INVALID --
from_address rejects it (a plain window under start=no; 'free' is
not a mask token, '-' is). Fixes: opener rule-only; the grammar
section carries ONE REAL address decoded slot by slot
(ti_3d4sno2_3-1_sa-casscf_no_-_pol_avg); the docs-example
round-trip pin in test_cells keeps it honest forever. core's Si/Sc
measured scales (a round-B grep gap: bare scales, no units) -> class
statements + record pointers; a fact-17 tag dissolved with them.
Boundary suite at the wave: 185 green. Fast tier 118.

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
- CONSUMED 2026-10-29 (phase 4): `_persist_cell` (cells.py) ->
  `_store_cell` (the store/load family; one caller, renamed with it).
- (add freely as spotted; an entry dies when consumed)

## Flags for later phases (spotted during phases 2-3)

- RESOLVED (phase 4): the ordering bridge lives in OPERATORS now
  (detbasis was the candidate but is deliberately numpy-only; operators
  owns the other PySCF convention boundary and sits below every
  consumer -- analysis's lazy anti-cycle imports became top-level).

- RESOLVED (phase 4): SA store pairing DECIDED-KEEP -- see the
  DECIDED-KEEP paragraph in the phase-4 entry above.
- RESOLVED (phase 4): `ground_block` lives in OPERATORS (the
  measured-block family beside window_closure/core_anisotropy),
  ops-first with vec returned -- natural's inline duplicate absorbed,
  the cells->study lazy import gone.
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
