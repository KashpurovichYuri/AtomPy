# Ticket 32 — The symmetry lane: D2h-pinned CI as the second route

## Goal

Add the independent route to L²-symmetric states: PySCF's
symmetry-pinned CI (D2h `wfnsym` solvers, `state_average_mix_` over a
term's components) evaluated and certified by our exact meters.
Two routes, one comparison — they never merge:

- the free route (ours): diagonalize the symmetry-blind Ms sector, MEASURE
  symmetry with the exact operators, restore by SA over roots;
- the pinned route (theirs): CONSTRUCT symmetry by group theory — one solver
  per D2h component, weights from the SO(3)→D2h branching.

The shared surface is the existing seam: the substrate signature
`(mol, mf, c_act, ncore, na, nb)`, the full-space CI-vector bridge,
the analysis meters. No new kernel path.

## Why D2h (probed 2026-10-31)

- An atom with `symmetry=True` gets SO(3), but PySCF's CI symmetry
  machinery is Abelianized there (m-sum blocks) — `wfnsym` cannot pin
  L. D2h pins the |M_L| classes properly; it is the working dialect.
- PySCF's D2h irrep ids are the Z2³ XOR encoding (Ag=0 … B3u=7;
  irrep product = XOR) — a determinant's irrep is the XOR over its
  occupied spatial orbitals; doubly-occupied pairs cancel.
- `wfnsym` solvers return FULL-space CI vectors (zeros on forbidden
  determinants) — the detbasis bridge consumes them unchanged.
  Demonstrated live: C 2p² triplet sector, the three g components
  degenerate to 1e-9, our meters reading L(L+1)=2 and S(S+1)=2
  exactly on every one.
- A configuration's inversion parity gates which irreps are reachable
  at all (p² → g only; s¹p¹ → u) — reachability is a contract, not a
  guess (PySCF raises `WfnSymmetryError` on the rest).

## Design

- New module `irreps` — the SO(3)→D2h dictionary, three functions:
  - `branch_d2h(L, parity)`: which irreps a term's 2L+1 components
    occupy, with multiplicities (character algebra; anchored at real
    shells). For L≥2 an irrep can host TWO components (e.g. Ag holds
    both z² and x²−y²) — the multiplicities are the SA weights and
    per-solver `nroots`.
  - `det_irreps(orbsyms, na, nb)`: per-determinant irrep id in
    detbasis order (the block identification for our own matrices).
  - `orb_irreps(mol, mo_coeff)`: the PySCF symm boundary (labels →
    XOR ids).
- The lane evaluator (pyscf-native CASCI/CASSCF +
  `state_average_mix_`), the cells method axis, and the cases
  term-targets come in the NEXT parts — after the dictionary is
  pinned. `AtomConfig` gains its symmetry knob when the production
  path needs it, not before (the fixture builds its own mol).

## Verification ladder (unit first, integration next, campaign gated)

1. Branch multiplicities sum to 2L+1 (L=0..8, both parities);
   hand-known p/d/f patterns; agreement with `label_orb_symm` on real
   shells. [unit]
2. The flip-derived character table equals PySCF's; the irrep
   product structure equals id XOR (all pairs). [unit]
3. `det_irreps` support equals the nonzero support of the
   `wfnsym=Γ` CI vector, per reachable irrep. [unit]
4. Doubly-occupied spatial orbital = the identity irrep (the (1,1)
   cancellation). [unit]
5. THE ANCHOR, on a frozen D2h-C fixture (repr-round-trip source, no
   SCF in the unit tier): per-irrep PySCF CASCI energies ==
   eigenvalues of our H restricted to the irrep block (to solver
   tol); our meters on their states read exact L(L+1), S(S+1); the
   term's components degenerate to ~1e-9; the full sector spectrum
   == the union of blocks. [unit]
6. SA head-to-head: our two-pass SA vs `state_average_mix_` CASSCF —
   ensemble energies agree, both densities pass closure certificates
   at the same floor. NEVER bit-compare orbitals across the routes
   (degenerate gauges — facts 10/15). [integration, part 2]
7. The d-shell stress case (Sc/Ti): the pinned route's guaranteed component
   set vs our nroots-based SA. [integration, part 2]
8. Campaign tables via the lane + the TABLES.md contract redesign,
   both records in hand. [gated on the user's enable]

## Parts

- Part 1 DONE (2026-10-31, commits e22260e + 16c3d4d): the `irreps`
  dictionary (branch_d2h/det_irreps/orb_irreps, character algebra
  anchored at PySCF's table and at real shells), the frozen D2h-C
  twin (generator `d2h` mode, double-run sha256; the C1 fixture
  untouched), and the unit pins 1–5 — THE ANCHOR at machine
  precision: pinned roots == our restricted-H eigenvalues (4e-16),
  term purity exact (2e-15), components degenerate (0.0), the
  spectrum == the block union. Review round done (two-axis); full
  suite 198 GREEN.
- Part 2 DONE (2026-10-31, the lean scope): the character-algebra
  and XOR math land in the docstrings (Notes carry the four chi
  formulas + the projection; det_irreps Notes carry the Z2^3
  argument); the lane evaluators in casscf — `run_irrep_casci` (one
  pinned block solve per component, the shared folded h1) and
  `run_irrep_sa` (state_average_mix_ + the two-pass readout); pins:
  the runner == the anchor on frozen C, and THE HEAD-TO-HEAD (C in
  the unit tier, Si in integration).
  - The head-to-head FINDING: ensembles identical (C 2e-12, Si
    ~3e-11 over three fresh draws); the pinned route's optimized windows
    close orders tighter (C: lz 1e-15 vs 2.3e-6; Si: 5.1e-8 vs
    1.5e-5) — the irrep-blocked rotations cannot leak the m-set.
  - PROBED PySCF fact: state_average_ on a symmetry-adapted molecule
    forces one wfnsym per solver (its wrapper passes wfnsym through,
    and None is GUESSED to one irrep) — the symmetry-blind SA exists
    only on C1, so the head-to-head CROSSES SUBSTRATES: the free route
    native (C1, the production driver) vs the pinned route native (D2h twin).
    Recorded in the test's comment; the lane evaluator docstrings
    carry the never-bit-compare-orbitals rule.
- Part 3: the evaluator lane is DONE (lean); what remains of the
  original part-3 scope: declaration (cases term-targets, cells
  axis) + step 8 (gated campaign tables + the TABLES contract
  redesign), plus ladder step 7 (Sc/Ti stress) whenever the TM
  references are enabled.

## Step 7 spec (the next session -- agreed 2026-10-31)

DONE 2026-10-31 (commits faf1adb + ff62a3d, review round included).
`tests/integration/test_symm_tm.py`, Sc + Ti parametrized, the Si
pattern (native substrate pairs: C1 the free route + a D2h-twin the pinned route,
fresh ROHFs under pytest tmp; any converged reference works -- the
whole active set sits in the window, so the fold sees only the
closed noble-gas core).  Windows/sectors/weights from the CAMPAIGN
rows (sector_na_nb; spec.sa_weights); the ground-term table is test
data (the test_atoms kind of row).  The spec as written:

One file `tests/integration/test_symm_tm.py` (unique basename), Sc +
Ti sections, the Si pattern (native substrate pairs: C1 the free route + a
D2h-twin the pinned route, fresh ROHFs under pytest tmp; the D2h constraint
makes the degenerate-shell reference deterministic, and any converged
reference works -- the m-scalar fold + the m-complete window absorb
the broken sphere).  Windows/sectors derive from the CAMPAIGN rows
(the one-declaration rule).  The stress content, beyond C/Si:

- **The multiplicity-2 branching goes live**: Sc's ground term ²D
  (L=2, g) branches `{Ag: 2, B1g: 1, B2g: 1, B3g: 1}` -- the first
  per-solver `nroots=2` through the runner; Ti's ³F (L=3, g) goes
  deeper (7 components).  The anchor pin covers it via block sizes.
- **nroots-hoped vs guaranteed components**: the free route's Ti SA hopes
  the 7 lowest (3,1)-roots are the ³F set (¹D can interleave); route
  B guarantees each.  The head-to-head is CONDITIONAL and honest:
  measure the free route's roots with our meters FIRST; assert the ensemble
  identity iff the root sets coincide; a differing ensemble is
  itself the finding (document, never force).
- **Certificates on the hostile atom**: the pinned route's mix window at/near
  the floor vs the free route at tolerance -- the C/Si finding re-tested
  where member tilts break free runs (fact 18).
- Pins per atom: runner anchor (blocks == our H, meters exact),
  component degeneracy, the conditional head-to-head, closures.

Declined up front: pinned-route single-irrep runs (fact-18 territory
already covered for the free route by test_treatments).  Suite work, NOT
gated; full suite at the boundary; the two-axis review after.

## Step 7 record (what was found)

- THE ANCHOR holds at machine precision on both atoms, and the
  multiplicity-2 branching went live: Sc's ²D Ag block hosts TWO
  components, Ti's ³F each B_g two (per-solver nroots=2, the first
  live exercise).  Pinned lowest-per-block roots == our restricted-H
  eigenvalues; components degenerate to ~1e-14 on the m-scalar fold.
- THE HEAD-TO-HEAD SPLITS BY ATOM, and the finding refines the
  spec's gate: MEMBERSHIP ALONE IS NOT ENOUGH.
  - Sc (d¹): the routes agree as on C/Si -- ensembles to 1.1e-8,
    multiplet formed (ptp 3.1e-7).  Bonus reading: the symmetric
    free-route window's certificate runs 5.2e-3 while the ensembles
    agree to 1e-8 -- the one-particle certificate sees f-direction
    drifts the energy cannot feel (hence _CERT_A = 2e-2).
  - Ti (d²): all 7 free-route roots still READ ⟨L²⟩=12 while the
    free optimizer had performed a CORE<->ACTIVE SWAP (see the
    mechanism addendum below): certificates O(1), multiplet split
    3.1e-3, span identity vs the raw shell 1e-4, ensemble 5.2 mEh
    BELOW the pinned route's.  The production Ti table records the
    same free-route basin (the swap row: comms 1.1/6.0, Elec
    s0 p2.0 d2.0).  The pinned route's D2h-blocked rotations FORBID
    the swap (it is cross-irrep: 4s is Ag, 3p is u): declared
    identity kept (s1.0 d5.0), components degenerate to 3.8e-7,
    certificates at 7e-6..1e-5 -- the guaranteed component set at
    ~5 mEh.  The symmetry-vs-variation trade, measured.
- THE HONEST GATE is the multiplet DEGENERACY (member AND formed);
  the divergence branch asserts ens_a < ens_b -- rigorous in BOTH
  shapes (free variation in the escape; displacement-below in the
  interleave: the intruder enters below the term root it pushes
  out) -- and never forces the finding's shape.  The free route's
  certificate in the divergence branch is measured, not pinned
  (basin chaos).
- Suite gates: TM file green on repeated draws; full suite 206/206
  GREEN (three documented rotating classes across the first three
  runs, each isolation-green).  Counts 206/132 in AGENTS; the
  seam-file exception now names the test_symm_* family.

## Mechanism addendum (2026-10-31, the user's 'why Sc but not Ti')

THE ESCAPE IS A CORE<->ACTIVE SEMI-CORE SWAP, measured on the
inactive block: the Ar start reads l-content s:3 p:6; the free
SA's converged inactive reads s:4 p:5.  The 4s pair rotated INTO
the fold; one semi-core 3p orbital -- with its electron pair,
0.998 pure 3p by S-metric AO content -- rotated INTO the active.
The gain is the 3p-3d semi-core correlation (the TM M-shell
effect); the price is a torn m-set (only one 3p m-function fits in
ncas=6).  The Sc/Ti asymmetry is the semi-core appetite scaling
with the open-d count (d1: the swap does not pay; d2: 5.2 mEh).
The pin (test_symm_tm head-to-head) measures the swap detector
and asserts the mechanism in every branch arm.

Lessons from the hunt (both caught by the meters/scales):
- The first hypothesis ('the optimizer wants the p shell; the
  m-complete version should pay') DIED: raw {3d,4p}/{3d,3p}
  probes came out ~300 mEh above with scrambled terms -- AND were
  malformed (fold-9 dropped the 4s pair: 20-electron problems).
  e_cas is fold-relative; only totals compare across folds.
- 'A full 3p shell in' (my first record fix) was also wrong --
  one column, not a shell; and the pair in it is the SEMI-CORE
  pair moving in with its electrons, not the 4s pair relocating.
- THE STEP-8 INPUT (gated): the m-complete symmetric version of
  the swap is the semi-core window CAS(6,8) {3d,3p} with fold
  {1s..3s,4s} on Ti (all six 3p electrons active) -- a candidate
  planned window alongside (4,10)/(4,14).  The radial-transplant
  (β) and NO-laundry (α) experiments stay queued; the cross-route
  window-intake decision is recorded in map.md Issues (D2h-born
  sources only).

## The CAS ladder (2026-10-31, the user's 'test different CAS' -- Ti, pinned, CI level, frozen shells)

All totals = component e_cas mean + closed_shell_core_energy (the
bookkeeping constant; fold-independent by construction -- fact 9's
fixed-core-block identity -- so cross-fold comparable).  One
substrate (the D2h ROHF), one fold rule; (4,·) windows PLACED
through build_full_mo (raw columns carry the documented ~85 mEh
Ti window-definition error -- the raw (4,14) run even read
<L^2>=13.66, a non-eigenvalue: the overlapping fold breaks [H,L^2]);
the (8,8) window is overlap-free by construction (its fold excludes
3p).

| window | E, Eh | vs declared | term |
|---|---|---|---|
| (4,6) {3d,4s} declared | -848.3912 | 0 | 3F exact |
| (4,10) {3d,4d} planned | -848.2045 | +187 mEh | D/H scramble |
| (4,14) {3d,4s,4p,4d} planned | -848.4349 | -43.7 mEh | 3F exact |
| **(8,8) {3d,3p} semi-core** | **-848.4520** | **-60.8 mEh** | **3F exact** |

THE VERDICT: the m-complete symmetric semi-core correlation (the
CAS(8,8) {3d,3p} window, fold {1s,2s,2p,3s,4s}, sector (5,3),
branch(3,g) -- an EXACT degenerate 3F, ptp 5e-14) is worth ~61 mEh
-- the free SA's torn escape (5.2 mEh) was reaching for a small
corner of a large pie, and the semi-core window beats even the
valence (4,14) window by 17 mEh.  Readings for the planned-cells
design: (4,10) as-declared (fold Ar(9), the 4s pair in the dead
zone) is NOT a 3F window -- its lowest g-states scramble D and H
components; the campaign's planned (4,10) cells have their own
constrained semantics and should not be read against this ladder.

NEXT for this thread (suite work + gated): CI-level pins of the
ladder (deterministic; the 17/61 mEh margins are robust), the SA
level on the winners ((8,8) first), and the campaign-declaration
question for the semi-core window (a new planned-window row).

## The definitive ladder session (agreed 2026-10-31 -- the next unit)

The user-enabled lean set: the four recorded Ti legs PLUS legs 5-7
and the Sc pair, all under ONE uniform construction so the table is
internally comparable (the first cut mixed mf-MO folds and raw
columns; the definitive table recomputes every leg the same way):

- EVERY leg: fold = Loewdin of the fold's raw AO shell columns;
  actives placed through build_full_mo(mf, core_idx, act, c_core=)
  (the seam ALREADY accepts arbitrary core columns -- no new
  machinery); pinned CI at frozen shells; totals = component e_cas
  mean + closed_shell_core_energy(fold) (the fold-independent
  bookkeeping ruler).  The uniform placement matters for ANY window
  sharing an l with its fold (4s<->3s, 3s<->4s, and even 2p<->3p
  cross-shell overlaps) -- the raw (4,14) lesson generalized.

| leg | window | CAS | fold shells | sector | branch |
|---|---|---|---|---|---|
| 1 | {3d,4s} | (4,6) | 1s..3p (Ar) | (3,1) | 3F |
| 2 | {3d,4d} | (4,10) | Ar | (3,1) | 3F (as-declared reading) |
| 3 | {3d,4s,4p,4d} | (4,14) | Ar | (3,1) | 3F |
| 4 | {3d,3p} | (8,8) | 1s,2s,2p,3s,4s | (5,3) | 3F |
| 5 | {3d,3p,4d} | (8,11) | 1s,2s,2p,3s,4s | (5,3) | 3F |
| 6 | {3s,3p,3d} | (10,9) | 1s,2s,2p,4s | (6,4) | 3F |
| 7 | {3d,3p,4s} | (10,9) | 1s,2s,2p,3s | (6,4) | 3F |
| Sc A | {3d,4s} | (3,6) | Ar | (2,1) | 2D |
| Sc B | {3d,3p} | (7,8) | 1s,2s,2p,3s,4s | (4,3) | 2D |

Leg 6 = the full M-shell classic (3s^2 included); leg 7 = does the
4s pair want back in-window on top of the semi-core; Sc B = the
d-count scaling MEASURED (never presumed -- the m-complete window
may pay at d1 even though the torn swap did not).  Declined: the
(10,14) {3d,3p,4s,4p,4d} ceiling leg (~3M dets; adds bulk, not a
new question) and 4f-class legs.

DELIVERABLES of the session: the pins (a test_tm_cas_ladder-class
pin: the energy ORDERING with margins + term exactness per leg via
the meters + the Sc scaling as measured), the TABLES records
(elements/Ti/TABLES.md gains 'the pinned CAS ladder' section beside
the two window-lattice sections; elements/Sc/TABLES.md the scaling
rows; provenance = the session + the certifying pins), gates (fast
tier per commit; full suite at the boundary).  THEN the unchanged
follow-ups: SA on the ladder winner, the campaign-declaration
question for the semi-core window (its explicit c_core fold needs
the declaration grammar to express non-prefix folds -- a design
point for the declaration wave), step 8 proper, the compression
discussion after it.

## The definitive ladder -- RESULTS (2026-10-31, the session run)

Spec correction first: leg 5 is CAS(8,13), not (8,11) -- {3d,3p,4d}
= 5+3+5 = 13 orbitals (the spec table miscounted; 368,082 dets in
the (5,3) sector).

| leg | E, Eh | vs declared |
|---|---|---|
| 1 {3d,4s} (4,6) declared | -848.3931 | 0 |
| 2 {3d,4d} (4,10) planned | -848.2156 | +177.5 mEh (D/H + quintet scramble) |
| 3 {3d,4s,4p,4d} (4,14) | -848.4374 | -44.4 |
| 4 {3d,3p} (8,8) semi | -848.4520 | -58.9 |
| 5 {3d,3p,4d} (8,13) | **-848.4991** | **-106.0 -- THE WINNER** |
| 6 {3s,3p,3d} (10,9) M-shell | -848.4549 | -61.8 |
| 7 {3d,3p,4s} (10,9) | -848.4529 | -59.8 |
| Sc A {3d,4s} (3,6) | -759.7274 | 0 |
| Sc B {3d,3p} (7,8) semi | -759.7871 | **-59.7** |

Every m-complete window hosts its term EXACTLY (components degenerate
to <= 3.3e-13; meters machine-exact; the direct block pattern
anchored against our block-restricted H at 4.3e-14).  The verdicts:

- THE WINNER is {3d,3p,4d}: the semi-core + the 4d pair-correlating
  shell (the 3p^2 -> 4d^2 channel, closed while 3p sits in the fold)
  DOUBLES the semi-core gain to -106 mEh -- at frozen shells, already
  85 mEh below the free SA's torn escape.
- THE D-COUNT CLAIM IS REFUTED at the window level: Sc's semi-core
  gain (-59.7) == Ti's (-58.9).  The 3p^6 correlation is
  d-count-INDEPENDENT; what was d2-specific was only the torn swap
  route into a 6-orbital window (the mechanism addendum stands as
  written -- it measured the swap, not the window).
- THE 4s PAIR STAYS IN THE FOLD: back in-window it buys 0.9 mEh
  (leg 7 vs leg 4) -- the escape's placement instinct confirmed, its
  route rejected.  The 3s^2 shell adds 2.9 mEh (leg 6).

Follow-up updates: the SA-level winner question is now
{3d,3p,4d} (not {3d,3p}); the declaration question carries the same
window with its semi fold.

## The step-8 workset spec (2026-11-01 — the declaration wave's design)

GOAL: the campaign tables (all four atoms) rendered from declared
cells through ONE grammar — windows carry their folds, scenarios
gate computation and demonstration, the lane tables generate like
the production tables — so TABLES.md carries the classified record
and nothing hand-copied can drift.

Three modules, each with its seam and its deletion test:

1. THE FOLD RESOLVER (windows.py, beside the selectors).  Interface:
   a fold spec string -> the explicit shell tuple; ``Ar``/``Ne``/
   ``He``... resolve through the noble-gas prefix table (universal),
   ``semi`` through the atom-class rule (open-(n)d: prefix minus
   (n-1)p plus ns), and a bracketed/leading ``1s,2s,...`` list is
   itself.  Unknown name -> ValueError (the seam owns the invariant;
   README "Fold notation" is the user-facing statement of the SAME
   rule).  Deletion test: without it, every fold consumer re-parses
   strings — it earns its keep at the second consumer (cases
   validation + cells construction).

2. THE SCENARIO LAYER (cases.py, ONE declaration).  A scenario is
   DECLARED DATA on the campaign with per-atom override (dict merge):
   ``{name: dict(cells=..., instead_of=..., reason="...<dated>",
   render="collapse", eps=1e-4, enabled=True)}``.  Two consumers
   interpret the generic vocabulary, nothing else: the RUN GATE
   (cells: a scenario'd cell is VOID-WITH-REASON at the key — still
   addressable, ``--full`` materializes) and the RENDER POLICY (the
   emitter: rows agreeing with their ``instead_of`` head within eps
   collapse into a generated note).  Unknown scenario kind -> the
   declaration is refused (fail at the seam that owns the grammar).
   A skip reason is a dated physics claim; ``--full`` re-tests it.

3. THE LANE EMITTER (scripts/, the build_*_tables family's third
   member; the TABLES edit protocol verbatim — structural asserts,
   the .md.tmp move).  The CI-level ladder legs become DECLARED
   planned cells: window strings gain ``@fold`` ("3d3p4d@semi"), the
   cells axis enumerates them like any planned cell (chk + payload
   per key; the fold is part of the address grammar), and the
   emitter renders the lane section between markers.  The existing
   test_cas_ladder_* pins REMAIN the certification; the emitter
   reports, the tests certify.  SA-level cells wait for the cost
   probe (368k dets x iterations) — not declared until priced.

THE TABLES SECTION GRAMMAR (per atom, four section types):
contract (combinatorial, hand, route-independent) / production (the
declared cells, as-run: the treatment-grid spine + generated
markers; the multiplet/excitation tables live here) / lane (the
generated ladder + the HAND head-to-head section with the mechanism
record) / readings (scan verdicts, open questions).  Every energy
row states partition (fold notation) + route (free/pinned) +
provenance (cells.json | pin | session).

THE TI COMPRESSION under the render policy: keep the certified
rungs, the basin-different rows, the detector pairs (~10 rows);
start-triplets that agree collapse to head + generated note; ALL
pending rows leave the committed table (the full grid stays one
regen away; a pointer line says so).

ROLLOUT: Sc first (smallest), one review pass, then Ti (the
compression), then C/Si.  The stale ``<!-- study:begin -->`` blocks
(the generator is retired) get a one-line frozen-history label in
the same pass.

PINS the workset owes: the resolver (each rule class + the explicit
list + the unknown-name refusal); the scenario evaluator (on/off,
the campaign-atom merge, VOID-at-key, the eps-collapse, the
unknown-kind refusal — pure data, unit-tier); the emitter's
determinism (same payloads + policy -> byte-identical section); the
``@fold`` cells (a declared semi-core cell's payload == the
corresponding ladder pin's numbers).

SKILLS: the design survey (codebase-design vocabulary) shapes the
seams as implemented; the two-axis review at the wave close.

REMAINING for the lane: step 8 + the declaration wave (cases
term-targets, cells axis, campaign tables via the lane, the Ti
TABLES compression with both records) -- GATED on the user's
enable; the ticket-31 interleave alongside.

## Step 8, first execution round (2026-11-02 — the grammar + the runner + the Sc rollout)

The workset spec's three modules LANDED, with its pins:

1. THE FOLD RESOLVER (windows): `fold_shells(name, open_d=)` — the
   noble-gas table (`NOBLE_CORES`, universal), the `semi` class rule
   (the noble core ending at np, minus np, plus (n+1)s; open_d the
   atom's open d shell), explicit lists (concatenated or
   comma-separated, must consume entirely); unknown names refuse.
   `closed_electrons(shells)` = 2 per m-slot (the lane sector
   derivation: nelec = nelec_act + 2·ncore − closed(fold); verified
   against the ladder's CAS(8,·)/(10,·) legs exactly).  Pins: 3 unit.
2. THE SCENARIO LAYER (cases): `Scenario` (name/cells/instead_of/
   reason/render/eps/enabled; render kinds 'keep'|'collapse', unknown
   kinds refuse at construction) + `scenarios_for` (campaign defaults
   ⊕ the row by name — DEFAULT_SCENARIOS empty today, the merge is the
   tested seam).  FIRST REAL USE: Ti's `anchor-first` scenario — the
   planned windows' optimizer cells (both label forms, plain + NO)
   are void-with-reason until their CASCI anchors are read; the gate
   is exact (58/66 planned cells void, the 8 anchors + the primary
   grid run).  Pins: 4 unit.
3. THE LANE (cells): the `@fold` window grammar (the tag is part of
   the address — a tagged window NEVER collides with a free-route
   label; NO windows take none), `lane_keys` (the sector derives from
   the leg's fold), `run_lane_cell`/`load_lane_cell` (the ladder's
   uniform construction productized — the certified fixture now RUNS
   the production path; its `_ladder_leg` deleted), `void_reason`
   (the run gate; `--full` materializes).  `CaseSpec.term` (L, parity)
   declared for C/Si/Sc/Ti; `lane_windows` for all four (C/Si join
   the lane with their default folds as explicit lists).  Pins: 3
   fast + the carbon lane run (term exact, store round-trip, cache
   hit, three refusals) + the emitter determinism.
4. THE EMITTER (scripts/build_lane_tables.py, the family's third
   member): the `<!-- lane:begin/end -->` block generated from
   cells.json under the TABLES edit protocol; the prose around the
   markers stays hand (the readings, the mechanism record).  regen
   gained `--group lane` + the void gate + `--full`; render() carries
   the lane/void payload shapes in their own sections.
5. THE ROLLOUTS (the first production lane runs, all four atoms):
   Sc −759.727436 / −759.787120; C −37.067793; Si −288.485053; Ti
   all seven legs — every total matches the certified record to its
   4-decimal precision.  THE DRAW-FREE PROOF, stronger in production
   than in the fixture: the Ti run's D2h substrate SCF landed in a
   DIFFERENT ROHF basin (−848.2834 vs the fixture's −848.4065, ~0.12
   Eh apart — the multi-basin class) and the ladder totals are
   IDENTICAL anyway: only mf.h_core enters the construction, and it
   is draw-independent.  All four TABLES lane blocks are now
   GENERATED (as-run provenance; the markers bound the generated
   block, the prose stays hand).  The stale study blocks in all five
   element files carry the frozen-history label (the generator is
   retired).

Declined: the Ti compression as a scenario (it is presentation —
the emitter's row policy — not a physics declaration; the spec's
"render policy" wording read as data was over-literal).  SA-level
lane cells wait for the cost probe, undeclared.

## Step 8, second execution round (2026-11-02 — the grammar amendments + the Sc lattice/lane runs)

User review of the Sc shape (2026-11-02) settled: prose minimal with
README the one home (DONE round 1); atom files carry ATOM-SPECIFIC
output only (no cross-atom tables in element files — the d-count
verdict moves to this record); scan sections leave the TABLES (the
anatomy delivery is rethought in ticket 33); the method-ladder naming
adopts the CAS abbreviation; sa_decomposition.md retires and scan.md
deletes at the C/Si pass.

1. THE LATTICE BECOMES DECLARED CELLS: `CaseSpec.lattice_sector`
   (validated; the series' own Ms — WHAT it studies is a declaration,
   not the hardcoded free-ion (2,0)).  Sc's valence series CAS(3,·) at
   (2,1): 8 spans x (AO, MO), 16 cells run, the table GENERATED from
   cells.json (the old session-frozen hand table dies; the MO column
   now full-rank — the ⟨L²⟩ control readings complete: 0.006 at {3d}
   through 5.6 at the top).  The lattice emitter generalized (per-atom
   sector + grammar heading/tag) and pinned (determinism, pending
   cells, declared-sector CAS labels, well-formed separators).

2. THE OWN-SPAN NO ARM RETIRED from the lattice (measured): its label
   (k = shell count) presupposes a TERM-AVERAGED source density.  The
   free-ion d² source delivers one (the ground block IS the whole
   7-root ³F multiplet: occupations 0.4 x5, ONE group).  The valence
   d³ source does NOT: its ground block is a SPLIT 2-ROOT ±PAIR
   (occupations [.5,.5,1,.5,.5] — a |M_L|±pair average, m-unequal at
   the fold-anisotropy scale) and the occupation groups split by noise
   ([1,1,2,1]) — the certificate was series-specific, never
   universal.  The NO window variant lives where its donor is
   declared per key (the methods grid).  Ti's three stored ti_*no*
   cells stay (never pruned); the arm's keys leave lattice_keys and
   its pin.

3. THE Sc LANE EXTENDED to the five-leg d¹ family (user-enabled; all
   legs ²D exact, components degenerate <=1.2e-13):

   | leg | E, Eh | vs declared |
   |---|---|---|
   | 3d4s@ar anchor | −759.7274 | 0 |
   | 3d3p@semi | −759.7871 | −59.7 mEh |
   | 3d3p4d@semi | **−759.8409** | **−113.4 mEh** |
   | 3d3p4s@{1s,2s,2p,3s} | −759.7899 | −62.4 mEh |
   | 3s3p3d@{1s,2s,2p,4s} | −759.7900 | −62.5 mEh |

   READINGS (atom-specific): the pair-correlating 4d shell DOUBLES
   the semi-core gain at d¹ (−113.4) — the 3p²→4d² channel; the 4s
   pair back in-window costs +2.8 mEh against the semi fold (it wants
   to stay folded); the 3s² shell adds nothing beyond the 4s
   placement (+2.9 over the anchor leg alone).  The cross-atom
   statement (the d-count independence of the 3p⁶ correlation and the
   pair-channel doubling at both d¹ and d²) lives HERE, in the record
   — not in the element tables.

4. Grammar amendments landed: 'the method ladder, CAS(3,6)/{3d,4s}'
   naming; the frozen-scan sections out; the studies blocks' broken
   ASCII-pipe headers repaired (‖…‖ — they never rendered as tables);
   the lane reading atom-specific.

REMAINING: the Ti pass (the methods compression + the renderer row
policy, headings lockstep with the emitters) and the C/Si pass (+ the
sa_decomposition.md retirement and the scan.md deletion), then the SA
cost probe on {3d,3p,4d} (gated — SA cells undeclared until priced),
the wave-close two-axis review + design survey + the refactoring
round (user-strengthened), then ticket 31 (U).

## Step 8, third round (2026-11-02 — the Ti pass: compression + grammar)

[Record repair 2026-11-02: this and the fourth round's ticket edits
silently no-opped at their commits (an unstated-em-dash mismatch in a
script replace; the commit messages carried the content until now).]

The methods renderer's ROW POLICY (the compression, renderer-local --
declined as scenario data in round 1): run cells only (pending and
void cells leave the table); within a (method, mask, substrate, fold)
group the AO head absorbs MO/NO runs agreeing within 1e-4 Eh -- named
in its start cell, "(+mo,no)"; disagreement keeps its own row (the
basin and construction classes stay visible); a window with no run
cells renders the POINTER line (declared, gated -- the anchor-first
scenario named; one regen away).  The policy lives in ONE home
(_groups/_collapsed, shared by rendering and the structural asserts);
window_section extracted for pinnability.  THE PIN:
test_method_ladder_row_policy (byte-determinism, full collapse,
disagreement, pending/void drop, pointer branch).

Ti TABLES re-sectioned to the grammar; all three emitters re-rendered
(method ladder 29 run cells -> 20 rows + 2 pointer lines; lattice
6/6; lane 7 legs byte-stable); net -91 lines.  THE STORE FACT the
compression exposed: the (4,10)/(4,14) planned windows have ZERO run
cells (the 8 'anchors' never ran); both windows render pointers.

## Step 8, fourth round (2026-11-02 — the C/Si pass + retirement)

C and Si re-sectioned into the grammar (the cross-check legends die --
README carries the key once; the basin map labeled a READING with its
session source).  Si's multiplet table carried a DEAD POINTER (values
in results/*_report.md -- the dir does not exist; the values were
nowhere in the tree): the staged pipeline ran for Si on the cached
reference and the Production table carries the real values (3P
-0.667376 / 1D -0.607963 / 1S -0.518843; report.md+json committed).
The scan sections are gone from every atom; RETIRED the artifacts:
scan.md x3 and sa_decomposition.md x3 deleted (regenerable by
scripts/radial_scan.py where ticket 33 needs them; sa_decomposition's
undelivered content -- the per-rung attribution column, the
tilt-follow curves -- is ticket 33's requirements list).  The lane
emitters re-rendered byte-stable on both rewritten files.

## Step 8, fifth round (2026-11-02 — the user's table review: A-prime)

The review's directives, executed: the OPEN-QUESTIONS sections leave
the atom files (their content is cross-atom/general -- new tables
deliver readings when their time comes); the cross-check COLUMN
leaves the method ladders (README keeps the what-it-is/how-to-deal
note; values live in the store); article-free headings everywhere;
the LANE vocabulary renamed LADDER (the section carries window
families x routes -- a route-named type would lie once the free rows
arrive).

THE LADDER TABLE redesigned (the meters speak, no verdict column):

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh |
  | E, eV | ΔE, eV | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |

- ΔE (was "vs declared") anchors on the primary leg; eV twins as in
  the method ladder;
- route/method columns carry the treatment identity (pinned/CI today;
  the free-route rows arrive with the route axis -- the user's "all
  we run with pinned we can run with free");
- ⟨L²⟩ renders ONE value when the branch components agree and the
  min-max RANGE when they scramble (Ti 3d4d reads 6.0-30.0 -- the
  scramble in meters, the verdict word retired);
- ptp = the component spread, the always-available degeneracy
  certificate;
- ‖[H,Lz]‖/‖[H,L²]‖: EXACT dense commutators, computed by
  run_lane_cell under a 6k-determinant memory cap (_COMM_CAP); past
  it they read em-dash (structural zeros for the certified fold
  class -- the note lives in README's ladder key).  All legs re-run
  fresh with the new payload fields; every total bit-stable vs the
  certified record.

## The SA cost probe -- PRICED (2026-11-02, unit B)

One timed pinned-SA solve on the ladder winner {3d,3p,4d}@semi
(368,082 dets, sector (5,3), branch(3,g) = 7 components,
state_average_mix_ + the two-pass readout), through the c_core fold
passthrough added to run_irrep_casci/run_irrep_sa for exactly this
(the non-prefix semi fold; default None keeps the prefix contract):

- WALL TIME: 6 min 52 s end-to-end (D2h ROHF ~1 min + construction +
  the SA solve) -- the solve itself ~6 min on the biggest leg;
- RESULT: ensemble -848.510765 Eh, converged -- 11.7 mEh below the
  CI-level term total (-848.4991): the orbital relaxation the frozen
  shells could not take, at exact component structure (irrep-blocked
  rotations).

THE PRICED STATEMENT: the SA level is MINUTES-CLASS everywhere in the
declared ladder families -- the 368k leg ~6 min, Sc's 204k leg a few
minutes, every other leg seconds-to-a-minute; a full route-axis batch
(pinned SA + free SA over all legs of all four atoms) is well under
an hour of compute per regen.  No hours-class cell exists.  DECISION:
the route axis is affordable -- declare it (free twins: CI / CASSCF /
SA through the free driver with the fold-tag resolution; pinned SA
rows through the lane runner's method axis).

## The route axis, EXECUTED on Sc (2026-11-02 — the twins run)

The user's directives: run the free variants ("all we run with
pinned we can run with free"), plan the full variant set in the method
ladder, ROHF at 4 digits.  Machinery: run_cell's FOLD BRANCH (the
tag names the partition; the fold's raw-shell Loewdin span is the
substrate; the SUBSTRATE axis names the route -- pol = free driver,
pur = the pinned lane class; validated BEFORE the store cache);
run_lane_cell's PINNED-SA branch (run_irrep_sa with the c_core
passthrough; meters on the two-pass readout roots); route_keys (four
twins per leg); regen --group route.

THE DENSE-BUDGET WALL, hit and fixed twice: the free twins on the
204k-det window entered run_cell's stage analysis (dense sector
matrices, O(ndet^2) -- 1 TB resident class).  Fixes: (1) run_cell
skips _facts past _COMM_CAP and the caller stores KERNEL-ONLY
payloads (the run's own energies; every dense meter NaN -- the
em-dash class); (2) analyze_cell (the CACHE path) got the same guard
plus fold-aware core sizing -- the first guard missed it and the
cached chk re-entered the dense analysis.

TWO CONVENTION BUGS the runs exposed: (1) run_irrep_sa's two-pass
READOUT folded the prefix (c_core not passed through to its internal
run_irrep_casci) -- the pinned-SA totals were ~4 Eh wrong; the
ENSEMBLES were always right (PySCF's own bookkeeping).  (2) The SA
total is the KERNEL ENSEMBLE, not readout+bookkeeping: the optimizer
rotates same-irrep core/window columns (4s is Ag; the semi windows
host Ag members) -- the readout's absolute energies are of the run's
OWN partition (the partition-bookkeeping rule); its component spread
and meters are the diagnostics.  Also load_lane_cell's fixed field
list predated the commutator fields (cache-hits silently dropped
them).

RESULTS (Sc, all 25 tagged cells run):
- SPAN IDENTITY at CI: every free twin == its pinned leg to 4
  decimals, INDEPENDENT meters agreeing (7.1e-13 vs 7.6e-13).
- ROUTE AGREEMENT at SA: pinned == free on four of five windows
  (-759.8523/-.8523 on the winner (7,13); -.8218; -.8266 x2).  The
  ONE divergence is 3d4s@ar (pinned -.7397, free -.7495): the
  escape-class window, the free SA finding a lower basin -- the
  multi-solution class measured ON the lane.
- The winner's ladder: CI -.8409 -> SA -.8523 on BOTH routes (the SA
  relaxation +11.4 mEh); single-root CASSCF reaches -.8668
  (kernel-only at 204k; its (7,8) twin prices the class: L2=5.07,
  commutators O(1) -- the JT-tilted root, no term).
- TABLES: the survey renders 29 rows (twins interleaved); the method
  ladder 20 rows with the full route plan; ROHF literals 4 digits
  everywhere (the first pass's regex missed 6-decimal literals --
  redone directly and verified).

The pinned-SA payload: total = the kernel ensemble; ensemble,
per-component readout energies, ptp, meters stored.  The route pin
(test_route_axis_twins): the span-identity certificate, the twin
derivation, every refusal class, the pinned-SA shape.

REMAINING (the next wave, agreed 2026-11-02 -- design settled in
conversation, see also the handoff):

- UNIT P (persistence, small): the Sc method-ladder Readings paragraph
  OUT (the table speaks; row-level notes stay); lane payloads gain
  mo_coeff (the pinned route's converged orbitals -- currently LOST:
  restart + orbital-level route comparisons) and cycles (CASResult
  gains the optimizer macro-cycles); store round-trip pins.
- UNIT S (sparse meters -- the dense wall dies): the matrix-free
  apply IS our _build with the N^2 container deleted -- same pair
  classification (occupation table + prefix counts), same
  Slater-Condon values, same signs; the scatter writes into the
  output VECTOR.  The zeros were never evaluated (the builder already
  skips >=3-excitation pairs by Pauli-vanishing + orthogonality) --
  their cost was storage (333 GB at (7,13)) and dense algebra
  (matvec N^2, eigvalsh N^3).  FOUR verification legs: apply ==
  _build matvec; == _slater_condon_slow (the elementwise direct
  reference, already pinned); == secondquant's explicit a/a-dagger
  application (the algebraic-direct); PySCF contract_1e/contract_2e
  the C leg.  Meters beyond the cap: L2/S2 per root + the STATE
  residual ||[H,Lz]psi|| (replaces the operator Frobenius norm at
  size; exact zero <==> symmetry eigenstate); dense stays below the
  cap, cross-certified.  README · Methods gains the matrix-free note
  + the excitation-class table.  Then re-run the big Sc cells: every
  em-dash fills.  nnz grows LINEARLY in N (dense quadratic): (7,13)
  3e8 couplings/apply sub-second; (7,20) 240 TB dense DEAD,
  apply seconds; (7,27) a minute.
- UNIT L (saturation -- the user's open question, NOT stateable from
  current data: at (7,13) CASSCF is still 14.5 mEh below SA): test
  whether a complete-enough window kills the tilt (single-root lands
  ON the term; SA == CASSCF; CI gain/shell -> 0).  Legs: S1
  {3d,3p,4d,4f}@semi (7,20) ~5.5M dets; S2 {3d,3p,4d,5d}@semi (7,20);
  CI first (minutes-tens), then ONE optimizer pair on the
  still-gaining leg (hours-class); S3 the union (7,27) ~51M only if
  needed.  AFTER Unit S (the meters are the point).
- UNIT T (Ti): the survey merge (emitter on Ti: lattice + 7 legs +
  variants + twins ~30 rows); the method-ladder RENDERER v2 (Ti is
  generated): method|CAS|fold|route columns, pending rows for
  declared cells; the Ti route batch (28 cells; the (8,13) SA pair
  ~6 min each); planned windows keep pointer lines (pending rows if
  the user prefers).
- C/Si passes, then the WAVE CLOSE (Unit C: two-axis review over the
  unreviewed span, design survey, the DEDICATED refactoring round,
  map dispositions, ticket 32 resolved; ticket 33 next).

## The review wave, applied early (2026-11-03 -- the skills run at the plan door)

The user directed the two-axis review + design survey BEFORE the next
units (plan-accurately).  Fixed point 75d3a93..HEAD.  Three commits:

- 0f457c1 (cells): _fold_span ONE home (was written three times --
  run_cell, analyze_cell without the refusals, run_lane_cell);
  run_lane_cell's CI branch calls run_irrep_casci (the inline loop
  was verbatim _irrep_solvers) and SA/CI share ONE readout/meter
  loop (the block-support assert now covers both; anchor stays
  CI-only -- the SA parts are the start window's); _LANE_FIELDS one
  field list (writer+reader); __all__ lockstep (kernel_cell_payload,
  route_keys); run_cell's Returns states the None-past-_COMM_CAP
  contract; route_keys calls lane_keys once.  Certified: the
  lane/route pins green (span identity 1e-8).
- 5604460 (cases): Scenario leans to its one consumer -- the render
  vocabulary (instead_of/render/eps/collapses) had ZERO consumers
  (the render-policy consumer was declined at round 3; the dead
  vocabulary died).
- ca98fc8 (emitters): scripts/_tables_common.py (constants + label
  helpers + the marker-form edit protocol; the run-payload predicate
  joins its family as payload_of); the lane and lattice emitters
  SKIP survey-carrying atoms (Sc was a hard sys.exit -- a declared
  emitter broken against a declared atom); methods._values guards
  the kernel-only payload class.  Byte-identical output certified
  against the committed TABLES (the pre-refactor reference).

DECLINED/DEFERRED with reasons: the lane:->ladder: marker rename
(moot -- the sections die at the C/Si migration); cells.json NaN
literals (Python-tolerated; writer and reader are both ours); a
cells.py split into lanes.py (hypothetical seam -- one consumer);
stage-adapter retirement (live consumers; correctly scheduled with
run_u*/ECP); moving the eps-collapse row policy into scenarios
(presentation lives with the emitter -- the earlier decline
re-affirmed by the deletion direction); exposing the apply's
excitation indices (per Refactor 4 the closure IS the product; a
public pair-handle widens detbasis for one consumer).

The SA-total==ensemble rule gained its pin (test_route_axis_twins).

## Units P + S EXECUTED (2026-11-03 -- persistence + the sparse meters)

UNIT P (84a1afb): CASResult.cycles (the per-cycle callback IS the
instrument -- PySCF 2.13 stores no count; mcscf.CASSCF is a factory,
mc1step.CASSCF the class, CASCI not its subclass); the lane SA
payload persists mo_coeff (the converged FULL MO -- was LOST) +
cycles; the free-route chk stores cycles; Sc's Readings paragraph
deleted (tables speak).

UNIT S (e9cdd8f, f580b65, fbb541b): the dense wall dies.
- detbasis.operator_apply: the SAME signed couplings the matrix
  builder scatters, accumulated into the output vector (bincount) --
  bitwise-certified refactor (_one_body_pairs/_two_body_pairs = the
  one home for the rule sets; _build rewired, old vs new module
  bitwise-equal on matched inputs).  One-body-only closures ride a
  SINGLES-ONLY enumeration (bitmask searchsorted, O(n*nelec*nvirt)):
  1.9s build + 0.09s/apply at 204k.  The bra/ket sides: the SOURCE
  determinant supplies the sign; the scatter lands in the partner
  rows.
- operators.pyscf_apply: H through PySCF's C kernels at any size
  (absorb_h1e at the kernel's OWN .5 factor + contract_2e; the
  permutation owned inside; complex vectors via re/im halves).
  Conventions MACHINE-CERTIFIED, recorded: contract_1e packs the
  lower triangle (silently symmetrizes -- Lz would be corrupted; the
  one-body C leg does not exist); contract_2e alone is NOT our
  two-body (the E-product carries a one-body delta); the C kernel
  serves only the Coulomb-symmetric class (L2's composite two-body
  cannot ride it -- ao2mo.restore assumes the 8-fold symmetries).
- Meters past the cap: lane payloads + kernel_cell_payload(key, mol,
  mf) carry res_h_lz AND res_h_l2 (max over component roots of
  ||[H,O]psi||, O applied through operator_apply -- Lz as its real
  form M = -i Lz; H through pyscf_apply).  The lane meters moved to
  the READOUT window (SA: the converged one).  Emitters' symmetry
  columns: dense norm below the cap, the stored residual past it,
  em-dash if neither.  PRICED: the two-body closure build at the
  204,490-det sector = 1208s once per cell (4.3s/cached apply) --
  executable at 2e5; at 5.5e6 the n^2 scaling is dead (days) --
  Unit L's L2 meters need either a smarter doubles enumeration or
  the one-particle closure certificate there.
- regen --fresh (the store-migration flag).  README: the
  matrix-free note + the excitation-class table; the de-slang round
  (leg/arm out of the one-home; the E key states comparability
  first; the meters key leads with what zero certifies).

REVIEW WAVE applied early at the user's direction (0f457c1, 5604460,
ca98fc8, f41bb51): _fold_span one-home; the run_irrep_casci adapter
+ ONE unified readout loop; Scenario lean; _tables_common.py; the
Sc lane/lattice emitter exit FIXED (the survey skip); methods NaN
guard; __all__ lockstep; the SA-total==ensemble pin.

USER DIRECTIVES (this conversation, standing): tables carry NO
slang; E's comparability rule first; SA totals are the uniform-weight
ensemble (the degenerate limit of the mean -- the readout absolutes
move under same-irrep partition rotations); the residual meters test
THE STATE (apply on the solution); Ti gets EXACTLY Sc's policy (no
atom-specific cases; renderer v2 declaration-driven, any atom --
windows are per-atom DATA, rows not code; planned windows keep
pointer lines, never pending rows); execute and price, don't
pre-decline.

Settled in conversation (the standing answers): the SA route
agreement is genuine (1e-9..2e-7 at solver tolerance on m-complete
windows; 9.77 mEh basin split on the escape-class 3d4s window);
iterations were not stored; restart = converged-state chks (PySCF
per-cycle dumps exist via mc.chkfile); lane mo_coeff was a gap.

## The store migration EXECUTED (2026-11-04 -- residuals + conv everywhere)

The batch ran chunked after TWO operational bugs (both fixed in
d380aa7/36be795): (1) the shell host kills background commands at
exactly 2h wall -- batches run in <=90-min chunks since, with a
user-directed WATCHDOG CYCLE (sleep+check: process alive, log
advancing, store fresh; kills cost only the in-flight cell now);
(2) regen wrote cells.json at group end -- any death lost the whole
invocation (twice).  The store now writes PER-CELL.

RESULTS (all committed, 238de78):
- Every certified total reproduces: 7 Ti legs + both Sc routes
  (route gap -0.00 mEh at the (7,13) SA winner) + C/Si bit-match.
- The em-dashes fill: state residuals at the apply floor on every
  pinned leg (1e-15..4e-13) and independently on the free twins;
  free <L^2> columns filled (6.0000 at the m-complete windows).
- THE JT READING: the (7,13) free CASSCF root -- res_h_lz 0.25,
  res_h_l2 1.31, <L^2> 4.86, conv 1.1e-06/200 (cap hit).  Ti's
  scramble window (3d4d@ar) carries ITS verdict in ptp 5e-02 with
  floor residuals -- eigenstate-ness (residuals) and term identity
  (ptp) are separate certificates and the tables now show both.
- THE CONV COLUMN (user directive '1e-x/77'): _cycle_trace captures
  the final |dE| beside the count; CASResult/StageFacts/chk/kernel
  and lane payloads persist it; conv_cell renders it (methods +
  ladder tables; README keys).  First texture: the escape-class
  free SA 1e-9/75 vs pinned 1e-13/32; the JT root cap-hits at 200.
- Pricing measured: the 204,490-det two-body closure 23 min in-cell
  (368k ~66 min) -- the lane SA legs carry it; fine at this scale.

REMAINING (unchanged): Unit L (saturation legs; CI legs priced then
ASK before the hours-class optimizer pair; L2 meters dead at 5.5M by
the n^2 measurement -- one-particle certificate or smarter
enumeration); Unit T (declaration-driven renderer v2 -- the methods
WINDOWS assert is the live boundary; survey merge; Ti route batch --
the anchors it needs are NOW RUN (33 cells)); C/Si; the wave close
(the full suite STILL OWED at the boundary -- physics-bearing
modules changed this wave).

## The ONE-TABLE grammar round (2026-11-04 -- the user's design review, Sc settled)

User-directed, landed (46c8229, 056e718, cc1b5f7): the element file
BECOMES the one table.  elements/<Atom>/TABLE.md = the atom heading
+ the ladder between ladder:begin/end markers; TABLE.csv the mirror
(same rows, ASCII headers, csv-quoted).  The user's directives, in
order: ONE big table (the per-CAS sections died); method first,
route second; the method column bare + restriction only (CASSCF,
radial); the START rides the window cell ('3d4s / AO', slash); NO
starts are not ladder rows (cells stay runnable, out of the render);
CASCI BESIDE its optimizer variants per fold group (CI interleaved,
the survey/lane sections die with the atom's migration); ptp
INTRINSIC (stored pinned, computed at render on free SA root sets;
dash = no set); the NOTE column (terse data-derived: 'cycle cap' at
cycles=200); THE CLOSURE COLUMNS ||[P,Lz]||, ||[P,L]||
(window_closure on the readout window -- run_lane_cell +
kernel_cell_payload carry them; the MO-slice rows read 3e+00, the
certificate working; 10 pre-field free-twin payloads dash until
their next re-run); the fold markers -raw/-frac/-pur/-lane defined
at Fold notation (the one home); the FURTHER-VARIANTS bullet list
(the diagnostic levers, e.g. on bad convergence); eV twins KEPT (the
user's call); the collapse row policy DELETED (it hid conv/meters).
Sc: 66 rows (29 CI + 16 CASSCF + 20 SA + anchor), 17 family rules
(8 lattice spans + 4 fold groups + 5 lane families), deduped by
address.  The generators skip one-table atoms by file; the ladder
builder refuses atoms without TABLE.md.  The splice eats everything
ladder-heading-to-boundary (historical forms cannot survive).
Bonus reading: the pinned (9,9) SA's converged span is NOT perfectly
m-closed (PLz 3e-05) -- irrep-blocked rotations keep the state
symmetric while the span drifts; exactly what the new columns read.

NEXT (user-ordered): UNIT U (UHF as a reference class) BEFORE the
Ti/C/Si migration -- run_uhf (character-forced alpha/beta
occupations, the basin gate -- the worse multi-basin class --,
guards, chk trust), spin-paired cores first cut (broken-core fold
generalization a separate declared step), the <S^2> contamination
meter on the reference itself (the exact S2), the conv trace; open
seam question FIRST: how runs declare the UHF reference (a
frac-class fold variant vs its own substrate class).  THEN: Ti's
route batch + the TABLE.md migrations (Ti, C, Si), the survey/lane/
lattice emitters dying per atom, the boundary suite at the close.

OPERATIONAL (this session's lessons): the shell host kills
background commands at exactly 2h -- chunks <=90 min, per-cell
store writes, the user-directed watchdog cycle (sleep+check on the
IN-REPO log); pkill -f self-matches its wrapper (use [b]racket
patterns); detached suites via setsid survive the cap.

## Unit U EXECUTED (2026-11-04 -- UHF as a reference class, user-ordered before the migration)

THE SEAM DECISION (settled first, as planned): UHF is a
REFERENCE-LAYER class, NOT a cells axis.  Both candidate axes die on
shape: every fold path consumes ONE closed-shell core density
(active_space_integrals: d = 2 C C^dag) and one MO set
(occupied_core_idx/build_full_mo/constrained_casscf all index
mf.mo_coeff[:, ...]); the frac precedent does not carry (the frac mf
is RHF-shaped and drops in unchanged, the UHF mf is (2, nao, nao)).
Consuming a UHF reference downstream = the declared BROKEN-CORE FOLD
step (the alpha/beta c_core generalization).  Recorded in the
module docstring + run_uhf Notes + AGENTS row.

WHAT LANDED (src/atomterm/reference.py):
- _forced_uhf_occ: the per-spin character-forced hook (the same
  ReferenceSpec vocabulary; at maximal-Ms sectors the beta open
  block is empty -- n_pair pairs + mol.spin alpha singles exhaust
  the electrons; the 2:1 premium unchanged);
- _accelerated + _converged_cached: the cache/gate skeleton
  EXTRACTED from run_rohf (the deletion test: ~25 lines x2; the
  SCF call sequence identical -- behavior-preserving);
- run_uhf: the UHF runner on the shared skeleton; the gate note:
  an occ_e_ref from an ROHF census gates the UNBROKEN basin only --
  a broken UHF sits BELOW it and trips BY DESIGN (anchor on a UHF
  census or leave unset while census-ing);
- uhf_contamination(mf) -> (s2, sigma): the overlap formula
  <S^2> = Sz(Sz+1) + N_beta - sum|<phia_i|phib_j>|^2 + the
  corresponding-orbital pairing spectrum (singular values of the
  occupied overlap; 1 = paired, 0 = broken).

CERTIFICATION (Verification rule 4, the three-layer chain):
- (a)+(c) unit pins, no SCF: the formula == the DIAGONAL element of
  the exact determinant-space S2 (s2_operator_parts through
  secondquant_matrix) on 0/1-overlap patterns AND the continuous
  broken pair (angle theta; the exact side = the one-body-rotated
  paired determinant, expm; sigma reads cos(theta); analytically
  1.75 - cos^2(theta) at the (3,2) norb=4 case);
- the paired fixed point: exactly S(S+1), sigma == 1;
- (b) integration (tests/integration/test_reference_uhf.py): C
  converges with the per-spin forcing, E <= E_ROHF + 1e-9 (UHF
  contains ROHF), s2 >= the pure-triplet floor, sigma in (0, 1].

MEASURED (the demo runs, in-repo temp chks, removed after):
  C  <S^2> = 2.006315 (pure 2); sigma = [1, 0.9968] (the 2s pair
     polarizes, 1s exact); dE = -4.127 mEh below ROHF.
  Sc <S^2> = 0.753166 (pure 0.75); sigma = [1x5 (the [Ar] core,
     EXACTLY paired), 0.999988/974/975 (3p), 0.999572 (3s),
     0.998908 (4s -- the MOST broken)]; dE = -1.498 mEh.
The pairing spectrum reads the spin-polarization profile SHELL BY
SHELL -- the [Ar] core intact, the valence s pair polarizing around
the 3d1 hole (the classic exchange polarization); both campaign
atoms break at the maximal-Ms forced occupations.

COUNTS: 245 tests (160 fast) -- the AGENTS 229/146 were STALE (the
one-table round added pins without the count update; 237/154 at
HEAD before Unit U).

REMAINING for U (the consumption wave, user-gated): the UHF anchor
rows in the tables (the census + gate numbers per atom), the
broken-core fold design (alpha/beta c_core), Ti/Si readings.

## The two-axis review + the consumption round (2026-11-04 -- the skills over f41bb51..HEAD)

The user-directed consumption round LANDED (6197f68): the census
(uhf_basin_census -- five arms, any atom one call, no chkfile),
occ_e_ref_uhf (the UHF gate anchor; the skeleton's gate energy
parameterized; _uhf_object the ONE construction), scripts/
uhf_anchor.py (the gated run + the store's uhf record), the builder
renders UHF (anchor) beside ROHF (DeltaE against the same zero,
<S^2>/sigma_min in the note).  Sc DECLARED: occ_e_ref_uhf =
-759.73789123 -- the census reads a SINGLE basin (all five arms
identical to 9 decimals, 13 s), <S^2> = 0.7620, sigma_min = 0.9945.
The per-commit TEST GATE re-settled (user): fast tier + targeted
files, 5-10 min MAX; the heavy fixtures once at the wave close,
detached, OMP_NUM_THREADS capped (the uncapped 255-thread run hit
an hour; the in-flight suite was killed under the new policy).

REVIEW FINDINGS APPLIED (this commit): comm_text deleted (consumer-
less); the redundant lattice_keys local import dropped; the builder
module docstring rewritten (it still documented the DELETED
row-collapse policy -- R6); the 63-bit mask limit GUARDED in
detbasis._singles_couplings (int64 shifts silently truncate past 63
spin-orbitals -- exactly operator_apply's advertised regime; the
limit is in operator_apply's contract); the displaced S(S+1) test
comment fixed; the lane_families helper's identical pinned/ci
branches merged; the anchor-row pin asserts the RENDERED row (V2);
the cycle-cap comment names the contract (200 = a DECLARED cap;
production runs PySCF's 50 default -- the payload carries only the
count, the cap-in-payload proper fix recorded in map.md); the
anchor-row energy formatting... SHARED via the minus loop (the
_quad extraction declined: two call sites, the deletion test says
pass-through).

DISPOSITIONS RECORDED (spec axis):
- Unit U's "the conv trace" clause: the CAS optimizer's _cycle_trace
  instrument does not apply to an SCF (no macro-cycle trace); the
  UHF record carries converged + conv_tol=1e-12; the anchor row's
  conv column dashes.  The clause is DISPOSED, not missing.
- kernel_cell_payload's residuals are GROUND-ROOT only (the ticket's
  recorded "max over component roots" is the lane's behavior); the
  max is a x nroots cost on the big sectors -- disposition at the
  wave close (narrow the wording or pay the cost).
- The SA budget guard (3adfe51) was never ticket-recorded -- THIS
  entry is its record (a declared block longer than the sector
  anchors on the measured one; PySCF crashes the never-run class
  otherwise).
- The ticket's "[P,L^2]" closure-column wording: the settled
  vocabulary is [P,L_perp] (fact 19) -- the ticket line is amended
  here, not in place (history stays).

## The migration sessions closed (2026-11-05 -- C/Si/Ti on the one table)

All four atoms now carry TABLE.md + TABLE.csv (Sc the exemplar):
- Ti's ROUTE BATCH ran detached (~5 h wall; the two 368k-det
  optimizer pairs the long poles; per-cell store writes carried
  every kill) -- 28 twins fresh, the lane families fill (5 rows per
  leg).  THE RECORDED VERDICTS REPRODUCE IN THE RENDER: CI span
  identity pinned==free on 6/7 legs; the EXCEPTION IS THE RECORD --
  3d4d's scramble window reads DIFFERENT STATES per route (the
  pinned blocks report the targeted term components, <L^2>=6.0-30.0,
  21.5 mEh from the free CI's sector ground: the recorded
  not-a-3F-window verdict); SA route agreement 5/7 legs at solver
  tolerance; the other exception: 3d4s@ar's pinned SA +7.2 mEh from
  the free route's escape (the core<->active swap is cross-irrep);
  the winner (8,13): CI -92.6 -> SA -104.3 both routes, CASSCF's
  tilt at <L^2>=10.20.
- The UHF anchors: Ti TWO census basins (default/huckel land +7.8
  mEh ABOVE the ROHF anchor; atom/broken land the spin-broken
  ground -3.4 mEh) -> occ_e_ref_uhf=-848.40991227, the gate cycles
  onto the winner exactly; C and Si single-basin (all five census
  arms identical), anchored ungated: C -4.1 mEh (<S^2>=2.0063,
  sigma_min=0.9968), Si -3.6 mEh (<S^2>=2.0110, sigma_min=0.9945).
  uhf_anchor's contract: the gate is the multi-basin GUARD --
  a single-basin atom runs ungated, the census verdict being its
  cover.
- FIX constrained_casscf (physics-bearing, pinned): a ONE-root
  weights list IS the single root (PySCF's SA wrapper einsums the
  energies; nroots=1 returns a scalar) -- Si's sa-casscf_no cells
  crashed exactly there; tests/integration/test_casscf.py pins it.
- The stores completed before the renders: C/Si gained their 16
  missing grid payloads + 4 route twins each (mostly cache hits).
- THE EMITTERS DEAD: build_{lane,survey,lattice}_tables.py deleted
  with their pins (the one-table builder is the only render path);
  the stray old-stage chk files deleted (elements/Si_*.chk at the
  root, Sc_rohf.chk/Ti_rohf.chk at the element roots -- the
  pipeline path is rungs/<address>.chk); elements/ecp/ STAYS (the
  ECP tests' live workdir).
- Docs lockstep: AGENTS conventions (the TABLE.md grammar), README
  (the worked example REFRESHED from Ti's render -- the old excerpt
  still showed the dead verdict-column format), test_symm_tm's
  record pointer, map.md's sector-column Issue resolved.

Counts: 253 total, 164 fast (the three emitter pins died fast-
marked; the casscf pin is slow-tier).

## The wave close (2026-11-05 -- review, survey, refactoring round)

- THE TWO-AXIS REVIEW DELTA (e165240..70d2568, applied 12e4773 +
  e9952b9): the slow mark on test_casscf; the empty-weights R5
  refusal + pin (observed |dE| = 1.4e-14 cited); the ladder-
  protocol wordings; the separator derivation; THE RECORD CORRECTED
  (CI span identity 6/7 legs -- 3d4d's scramble reads different
  states per route; SA agreement 5/7).
- THE FULL src/ DESIGN SURVEY (three lenses: cells/cases deep-
  dive, the 14-module health skim, the test/scripts hygiene pass):
  the architecture HELD -- no upward edges, windows the substrate
  root, natural pure at module level, reference's _uhf_object a
  real two-adapter seam, every cells/cases export consumed.  The
  findings fed the round.
- THE REFACTORING ROUND (6395de7; the full list in the commit +
  map.md's resolved entry): the two cells one-homes (_key_core,
  _apply_family), run_carbon reads the row, void_reason's
  double-check gone, the dead main() deleted, the max_cycle schema
  bump (identity flagging, legacy fallback), the verbose kwarg
  dead, the scripts one-homes + docstring truths + the cas_label
  NO branch deleted, the V1 splits (3-way + 6-way), the V7
  citations captured -- AND A GRAZE CAUGHT: test_symm_tm's l2 band
  widened to 1e-6 (observed worst 8.2e-9, Ti 3d3p4s SA leg; 1.2x
  headroom was the flake class).
- THE DELTA REVIEW OF THE ROUND (both axes over e9952b9..6395de7):
  R6/V6/V7 verified clean by the sub-agents (the extraction order
  reproduced verbatim; the max_cycle chain consistent store to
  builder; the SA roots path carries the optimizer mc).  One delta
  found + dispositioned: run_carbon's artifact stem aligns to the
  declaration ("C"; no committed artifact affected -- recorded in
  map.md).  map.md dispositions landed with this entry.

Counts after the round: 253 total, 164 fast.

## The suite verdict + the policy amendment (2026-11-05, user-directed)

The wave-close full suite ran detached and was KILLED at the user's
decision ~30 min in (the heavy tail re-priced live: the test_symm_tm
fixture carries the (8,13) SA-leg optimizer pairs, 1-2 h class --
the stale "~15 min" record predated their admission).  PARTIAL
VERDICT: 91 tests completed before the kill -- one F (test_symm_si's
head-to-head, the RECORDED rotating graze class) GREEN on its
isolation rerun minutes later; the other 90 passed.  The REMAINDER
WAIVED: the wave closes on the gates (the fast tier + the targeted
files, all green today, including every file the round touched).
THE POLICY AMENDMENT (AGENTS Tests, user-directed 2026-11-05): the
full-suite run RETIRED -- no wave-close full run, the gates close a
wave; heavy pins run isolation on evidence; THE BOUNDARY: tests
verify constructions at pinned small sizes, heavy compute belongs
to the atomic calculations (the campaign runs -- the project's
goal); test_symm_tm's SA legs re-scope to store-anchored
verification (the map carries the item).

## References


- The live probes and the mini-anchor of 2026-10-31 (session record;
  the demo script is the anchor test's skeleton).
- Block's with-pyscf docs (the `state_average_mix_` pattern); PySCF
  `fci.direct_spin1_symm`, `symm.param` (the XOR id convention).
- Ticket 31's gated U batches re-scope in this lane's light once the
  evaluator lands.
