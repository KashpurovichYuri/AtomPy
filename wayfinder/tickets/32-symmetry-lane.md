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

REMAINING: the SA cost probe on {3d,3p,4d} (user-enabled in
principle -- the priced statement first), then the route axis (the
free-route rows: CI / CASSCF / SA-CASSCF per leg -- the user's
"all we run with pinned we can run with free"; machinery = the
fold-tag resolution in the free driver + the c_core passthrough for
pinned SA), the wave-close two-axis review + design survey + the
refactoring round (user-strengthened), then ticket 31 (U).

## References


- The live probes and the mini-anchor of 2026-10-31 (session record;
  the demo script is the anchor test's skeleton).
- Block's with-pyscf docs (the `state_average_mix_` pattern); PySCF
  `fci.direct_spin1_symm`, `symm.param` (the XOR id convention).
- Ticket 31's gated U batches re-scope in this lane's light once the
  evaluator lands.
