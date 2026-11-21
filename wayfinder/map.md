# AtomTerm — the effort map

The navigator: what is live, what is open, where the worked examples
live. History is git's job — resolved tickets are deleted, survey
dumps are not kept.

## Tickets

The wave spec: [`specs/pt2-perturbation.md`](specs/pt2-perturbation.md)
— the PT2 analyzation's contract, REWRITTEN PySCF-mrpt-primary
(2026-11-19): the window-ladder diagnostic; the symmetry probes; the
delivery chain.

- [45 — the PT2 runs + PERTURBATION_TABLE.md]: COMPLETE 2026-11-20
  and RETIRED to git history (see the resolved record below); the
  declared-out classes' gated recoveries moved to the Disposition
  candidates; the 2026-11-20 review + the V&V instrument recorded in
  Issues.
- [46 — the TABLE.md physics pass](tickets/46-table-physics-pass.md):
  DESIGN-FIRST, user-gated — the render-level significance policy
  (the stores keep everything; the table keeps the physics).

Resolved (git history the archive): 45 the PT2 wave, COMPLETE
2026-11-20 (2026-11-19..20): the lane seam (`cells.readout_lane` --
the store's readout twin, kernel-ready per-block `mc` out;
`perturb.nevpt2_rows_on` -- the pricing split; `casscf.irrep_solvers`
-- the pin's one construction; `cases.lane_substrate` moved in) +
step 5's ladder on ALL SEVEN atoms (C 29, Si 37, Sc 52, Ti 46, the
uncs 26/26, U14e 32 cells; the combined gate logs/validate_pt2_all.log
-- 1259 roots certified) + the flagship CAS-size read (the corrected
ladder PREFERS (4,12): the (4,15) promotion's CAS gain is dynamic
double-counting, E(2) −341.7 → −336.0 mEh, the corrected total
worsens −50.7655 → −50.7624) + the seam's measured bands
(`readout_mc(ci0=)`, e_noise 1e-6, p_tol 1e-4 — the per-atom tail
sweep Si 2.2e-7 / Ti 3.9e-6 / Ti-unc 5.0e-5; the meters matrix-free
— the dense build hit ~1.35 TB on the M-shell pol twins) + the
wave's reading discipline in README · Analysis.  THE WAVE'S ANSWER:
the optimized ladder SATURATES (the optimizer repairs the poverty
class; the thin spans fuel and PT2 does not repair them — the
remainder is static); E(2) is orbital-sensitive,
rotation-insensitive.  CLOSED OUT with the U14e V&V instrument
(`scripts/u14e_vv.py`, the run instrument; logs/u14e_vv.log green
6/6 — the wave's answers re-derived standalone at machine precision:
E(2) delta 0.0000 mEh, ⟨L²⟩ 72 at 8.9e-15 free / 1.9e-07 pinned,
route delta 5.3e-09 Eh) and the 2026-11-20 two-axis review + design
survey (the fixes + declinations in Issues).  45 step 3 — THE SYMMETRY
PROBES (2026-11-19, `scripts/probes/probe_pt2_block.py`, the certified pur
cell `c_2p_2-0_sa-casscf_ao_-_pur_avg`): E(2) = −58.1340 mEh
(−1.5819 eV) IDENTICAL across the ³P Sz=+1 block; probe A (gauge
redraw, 8 unitaries × 3 roots): controls 1.2e-14 (⟨L²⟩ = 2 exact,
vHv = e_cas exact), spread 2.7e-11 mEh = **7.3e-13 eV** — the
engine's wobble floor; probe B: pure M_L states (−1/0/+1) exact
through the complex-Hermitian Lz block (Lz is anti-Hermitian in the
real p gauge — the .real symmetrization zeroes it; the machine
lesson), the REAL ML=0 component reads the same constant, ±1 are
complex-gauge (covered by A's transitivity).  VERDICT: clean —
per-root rows stand as stored; the degenerate-block gauge class is
CLOSED for SC-NEVPT2 at L=1.  THE NUMBERS: E(2) is ~7.5× the
EN-era preview (−7.7 mEh) — the H₀-choice axis is a big effect
(EN's exact-configuration denominators vs Dyall's orbital gaps);
the class breakdown: Sir (0)' −27.4 mEh leads, then Srsi −10.7,
Sijrs −7.2, Sij −4.9, Srs −4.5, Sijr −3.4; Sr/Si read EXACTLY zero
(perturber norms 0 — no N±1 coupling at this 2e/3o CAS on pur); the
mechanism (structural vs selection-rule) rides the ladder read
across cells/atoms.  The stored block sits 26.6 µEh
below ROHF at the same partition (fact 13's near-equality + the SA
orbital gain); the guard's first real bite: store-flat vs
reader-column CI shapes (fixed: amplitudes are the assertion).
44 the PySCF-mrpt adapter
(2026-11-19): `perturb.py` rebuilt at the PySCF boundary —
`readout_mc` (the run_casci twin without the checkpoint) +
`assert_stored_roots` (fail-fast; energies strict, CI vectors
gauge-aware through the degenerate block's joint projector) +
`nevpt2_rows` (the kernel's canonicalization prelude + the 8
Koopmans classes via mrpt's public functions; Sijrs/eris and the
per-class call-shape dispatch mirrored).  Certified: the fixture's
exact ³P block on the pur core; the RAW-core split pinned at unit
scale (4.4e-5 Eh, the anisotropy detector class); the gauge-tolerant
refusal (a random unitary in the block passes, a zapped root
refuses); the plumbing-vs-kernel twin — 1e-10 on the small CAS, and
AT CAMPAIGN SIZE THE ENGINE'S OWN WOBBLE FLOOR ~3e-7 Eh RUN-TO-RUN
(PySCF's kernel three-run record on O2: −0.16978547 / −0.16978542 /
−0.16978512 — thread-parallel reductions in the big contractions;
the pins sit at that floor); the O2/H14 anchors at 1e-6 of the
examples' retained values; two-root H2 root selection; the SA guard
(kernelled).  13 tests; counts 365/271.  Two unit-tier facts the
build caught: PySCF returns CI vectors as (n,1) column arrays (the
store convention handles it), and the kernel loop WRITES the
canonicalized CI into the shared mc.ci — our plumbing detaches the
list, the twin pins run on fresh readouts.  the PT wave REPLAN — DELETION
(2026-11-19): the spec's "external-anchor-absent" premise was FALSE —
pyscf 2.13.1 ships `pyscf/mrpt` natively (`NEVPT` = sc_nevpt,
strongly contracted; `dfnevpt2`; NO caspt2 in 2.13.1, master-only —
the env fact the original spec missed). The home-grown ladder DELETED
in one commit: `perturb.py` (the FOIS/EN rungs 0–1, the Dyall rung,
the density layer), `test_perturb.py`, `validate_perturb.py`, the
committed `logs/validate_perturb.log` — git history the archive. THE
BUG that forced honesty: `dyall_h0`'s virtual-range `eigh` mixed α/β
virtuals through degenerate eigenspaces (coupling norm 0.016 →
5.4e-11 measured on the carbon fixture) — the committed Dyall script
numbers were INVALID; the 44b density certificates stood but their
consumer was gone (the deletion test passes every rung). EN's one
future perspective — the EXPLICIT first-order wavefunction
(c_I = V_I0/(E0−H_II)) — archived with it; revival = the
wavefunction-level analyzation if the probes demand it. Architecture
now: PySCF owns the PT2 algebra (`mrpt.NEVPT` per root on the stored
roots — the bit-exact readout precedent); we own the adapter + the
meters + the probes + the render (README's Analysis chapter, 1adb82a,
is the reading discipline: the gauge ledger, the two probes, the
per-L expectations). Step 0 anchors: O2/6-31G CASCI(6,8) NEVPT
−0.169785475 (vs the example's retained −0.169785157, 3.2e-7), C2
SA-CASSCF(8,8) → 4-root CASCI → NEVPT(root=1) green in 0.3 s.
Counts 372/278 → 352/258.  43 the validation script (2026-11-18;
the artifacts deleted at the 2026-11-19 replan — this record
describes the DELETED chain, git history the archive; its three
findings carry: the ruler + the gauge class (now formalized in
README · Analysis) + the reach class):
`scripts/validate_perturb.py` — the run instrument (the run.py
precedent), the three-layer chain printed as certificate lines on the
declared C/Si SA scope (grid + planned, ao start, free-route: 14+18
runnable cells, 3 roots each): (a) twins at campaign size (production
coupling vs `pyscf_apply`, EN vs the dense Slater–Condon twin), (b)
the e0 external anchor (the store's PySCF CI energy; the E(2) anchor
absent — the pyscf-forge note prints), (c) whole-space exactness +
completeness (production == dense on a synthetic full-support CI) +
reach containment.  GREEN: C 197 / Si 253 PASS, 0 FAIL
(logs/validate_perturb.log committed).  THREE FINDINGS the
certification caught: (1) THE RULER — the analyzation must fold the
RUN'S OWN core field RAW (`spherical_core=False`, `c_core =
mo[:, :ncore]`): the stored energies are the kernel's values on that
field, and the optimized core's gauge p-admixture moves e0 at the
5e-5 Eh scale between the raw and averaged folds (fact 14 applied to
the post-hoc analyzation; with it the anchor sits at 1e-15).  (2) THE
GAUGE CLASS — the Sz=+1 p² CAS block is exactly degenerate: the
stored SA roots are arbitrary block gauges (support 1/3 dets on the
pur substrates — Davidson's draw), disclosed per root by the script's
gauge line; the meters read the stored draw.  (3) THE REACH CLASS —
on symmetry-adapted tensors (pur) selection rules zero couplings
INSIDE the brute reach (43 of 63): completeness is certified as
production == dense (independent constructions on the same tensors),
the brute count is the disclosed upper bound.  THE LADDER READ
(preview of 45): ground-root E(2) is FLAT across the ladder — C
−7.77..−7.79 mEh (compact/@r1/@s1), Si −6.9 → −5.5 (the @s2 decline)
— against the +250 mEh CAS spread: the external space supplies the
same small dynamic share on every window (cc-pVDZ's valence is thin —
the reading key must carry the basis-thinness caveat); 42 the module born
(2026-11-18; DELETED at the 2026-11-19 replan, git history the
archive):
`perturb.py` rungs 0–1 — the FOIS (the exact apply-pattern
enumeration beside `fois_indices`, the brute spin-preserving
completeness certificate) and the EN E(2) in two independent
constructions (`epstein_nesbet` matrix-free production vs
`epstein_nesbet_dense`, the dense Slater-Condon twin) + the
occupation-sum `diagonal_energies` (the antisymmetrized pair sum —
the symmetric direct form is WRONG on general arrays; the random
non-symmetric integrals caught it, the verification rule's
purpose). Meters: floor + w_ref in `PT2Result`. 7 tests: whole-space
→ 0, one-electron closed answer, the zero-denominator fail-fast, the
twins elementwise (complex non-symmetric, 3 shapes), the diagonal
twin, the cardinality meter, the carbon pin (frozen fixture, NO SCF:
e0 −1.35139 valence-space ruler, e2 −7.689 mEh, floor 0.054,
w_ref 0.996 at n_fois 140); 33 the symmetry-anatomy
deliverable — RETIRED UNEXECUTED 2026-11-18 (the user's call); the
idea stays specified: per-atom Readings sections with (1) the
m-symmetric θ-grid (θ | ΔE | ‖[H,Lz]‖), (2) the member-tilt ladder
exposing the ±pair blind spots (descend-with-zero-commutator = the
directions the sector commutator cannot see), (3) one README reading
key; everything regenerable from `scans.py`, verification = renderer
pins. Would have delivered: the explanatory record behind the
disclosure columns. Unneeded now: the verdicts it would visualize
already live in the map's settled facts, the element files, and the
README rules. Revival = re-open a ticket citing this entry (shape
decisions recorded: descent ladder ~10 rows over the full 31-subset
map; per-atom Readings; the Readings grammar home); 29 the Ti-first protocol; 30 the
module walk; 31 the ECP wave; 32 the symmetry lane (C/Si/Sc/Ti one
table, the route axis, the UHF anchors; closed 2026-11-05); 34 the
folded reference (steps 1–5 + the gap closures; the After's T-C done
as 36A/38, T-D as the 2026-11-16 docs round); 35 the
basis-preparation layer (the design landed, the honest prior table,
the census + the Newton-disease certificate); 36 the runs wave
(batches B/C/D/E ran — the stores carry 66–208 prior-treatment legs
per atom; batch A closed as ticket 38); 38 the U14e close-out (T1 +
T3a/b/c + T4, 2026-11-14..16); 37 the basis-repair bench, the @s
wave, the subspace-ROHF (Stages 0–3, 2026-11-12..16: the seeded-MOM
repair certified + integrated, the twins' series + @s legs 36/36 per
twin, `_unc_s_block` retired, the UHF anchor re-landed, the frac
class measured-refused, `subspace_rohf` integrated with the battery
green and the per-contraction ΔE certificate measured — the record in
the ticket's git history + elements/<twin>/logs/subspace.log); 39 the
U14e prior axis (2026-11-17: the series census measured-REFUSED the
@s half — the ECP substrate holds no member basin, every converged
draw is the head f³d¹; the @r1 legs ran and CURED the band-identity
poverty — +146.9/+149.3 vs the incumbent's +923.1, the exact ⁵L
(72.0000) at machine-zero closure on pur; the optimizer follow-on is
the disposition candidate below); 40 the U14e series REVIVAL
(2026-11-17: the refusal re-attributed — the weight hook's
character-only premium, not the landscape; the members landed under
the self-seeded MOM (`run_rohf(mom_seed=)` + `SeriesSpec
seed="self"`/`open_counts`, certified at ~4e-9, f²d² under
diis+newton); @s1 ran — +141.3/+143.7, the series average BEATS the
head's own block by ~5.6 mEh, exact ⁵L on pur; @s2 MEASURED-REFUSED
at construction — no second f radial, the 5f radial
configuration-inert).

## Issues

### Open

- **The u14e V&V family + the scripts regroup** (2026-11-21, the
  user's directive: several INDEPENDENT scripts, one per CAS, each
  readable top-to-bottom, PySCF-direct; metadata from one shared
  place only): `scripts/u14e/` — `_vv_common.py` (metadata, the
  trusted-rung loader, anchors, meters, ledger, the row store + the
  one `VV.md` table), the 4_12/4_15/6_13/uhf legs, and the old
  orchestrator re-formed as the thin driver.  The free routes are
  fully independent of the production construction chain (the
  window written out, certified elementwise against the production
  builder at 2.2e-16; SA/readout/PT2 called directly); the meters
  are the kept analysis; the pinned route and the 6_13 7s radial
  are marked production calls.  Green: (a) −0.0029 mEh, (b) 3.0e-06,
  (d) 3.3e-09 Eh, (a15) −0.0006, (c) fresh-vs-fresh +3.087 (stored
  +3.086) — logs/u14e_vv_*.log + elements/U14e/VV.md.  THE
  REGROUP: `scripts/probes/` (the declaration-provenance
  instruments; `probe_pt2_block` stays the gate's live dependency),
  `scripts/generators/`, `complete_u14e_table.sh` retired.  TWO
  DISCLOSURES recorded: the UHF-seeded SA from a RAW start descends
  broken (⟨L²⟩ dev 7.2e+01; the committed converged-start reading
  holds — the repair holds the converged gauge, a raw spin-broken
  start escapes); the production-distilled 7s radial is 0.52
  core-lean, its pur remainder 0.0092 of the window gram — the
  (6,13) window's 7s direction is compromised on the pur substrate,
  the class behind the stored (6,13) SA refusals.
- **The 2026-11-20 wave review + design survey** (ticket 45's close;
  the two axes + the codebase-design vocabulary over `src/`, fixed
  point a0fe3fe): fixed = the renderer's ground-block spread got one
  home (`_ground_spread_meh`); the dead `size` branch retired (no
  manifest declared it); the lane block-support `assert` → a
  fail-fast raise (rule 5); `key_core`/`lane_core`/`window_columns`
  promoted public (the store-free V&V instrument is the third
  consumer; the docstring refs already pointed there; `__all__` +
  the layout row + `kernel_cell_payload`'s row gap, one lockstep
  commit); test_perturb's band labels run (a)-(d); the committed
  PERTURBATION_TABLEs were one blank line short of the renderer's
  own bytes (pre-existing drift, surfaced by the re-render, lockstep
  restored); AGENTS/README test counts lockstep (376/282 — the
  handoff's 374 was stale at 67ab1da).  DECLINED (judgement calls,
  recorded): `lane_substrate`'s docstring carries measured instances
  (the repo's exemplars land measured values in Notes-style prose —
  `run_rohf`'s own Notes cite Ti numbers; the slot rule reads those
  as Notes, not contract); `lane_substrate`'s basin-steer retry
  branch untested (campaign chemistry, run-verified — a pin would be
  an hours-class test costume); `_meters`/`EH_EV` duplicated across
  the standalone instruments (the run.py precedent: scripts are
  self-contained run instruments); `_lane_block`'s 13-parameter
  clump (a block-context bundle type — a later refactor, not this
  close); `validate_pt2._probe`'s string-parse coupling to
  probe_pt2_block.py's stdout (gate + probe are committed together,
  one home); `validate_pt2`'s `FAILS` module global (script style);
  the (6,13) leg print-only in the V&V (the handoff's assertion list
  is (a)-(e); the (6,13) meters carry the certified L=9 disclosure);
  the UHF-seeded start MO is the symmetric window (the handoff's
  literal wiring — the reference feeds PySCF's kernel; a
  broken-orbital-start variant is a follow-up if the spin-broken
  START class is ever wanted).  THE V&V RUN: logs/u14e_vv.log green
  6/6 — E(2) delta 0.0000 mEh, ⟨L²⟩ 72 at 8.9e-15 free / 1.9e-07
  pinned, the (4,15) worsening +3.081 mEh (stored +3.086), route
  delta 5.3e-09 Eh, the UHF sibling printed (E −6.961 mEh, ⟨S²⟩
  6.229, all five occupied pairs broken); the (6,13) meters
  reproduce the certified L=9 disclosure.
- **The 2026-11-19 review's declinations** (the two-axis review over
  the PT2 wave, fixed point 7ab52df): fixed = probe_tmp.py scratch
  deleted; the AGENTS cases.py layout row (lockstep); the measured
  values moved out of assert_stored_roots' contract slot; the
  shape-mismatch refusal message (a numpy broadcast error used to
  replace the declared text); e_tol/c_tol force-pinned; the
  compression meter (ClassRow.e_per_norm) landed; the absolute
  fixture E(2) pin (observed −0.26575851 Eh — which reproduced the
  campaign poverty-class cell's E(2) from an independent substrate
  rebuild, ~1e-7 Eh); H2 build one home; probe-A onto _expect; the
  spec letter amended to the declared stored-roots form + the mrpt
  internal-prelude honesty.  DECLINED (judgement calls): the probe's
  _root_row private reach (one instrument; a promoted seam would
  widen nevpt2_rows' interface for a single consumer — the deletion
  test fails it); cell_substrate returning a bundle type (fights the
  declared substrate-signature convention); the _VARIANTS guard
  (the registry deferral is documented); readout_mc's unused `mol`
  (the substrate-signature convention).
- **The future-wave queue** (ticket 34's After + the standing
  candidates): the contraction-recipe provenance review (our radials
  vs the basis authors' coefficients); then spin–orbit, transition
  moments, g-factors.  (The term-ID round — splittings vs NIST —
  DECLINED 2026-11-18, the user's call.)
- **The C/Si pinned-lane legs — APPROVED 2026-11-18** (the 2+2
  pending rows, the last non-twin pendings): `--group lane` on C and
  Si; minutes-class; closes the lane families completely.
- **The queued CASCI review — EXECUTED 2026-11-17** (the three
  classes classified; no analysis bugs found): (1) the window-poverty
  class (+923 mEh) is the correct declared reading — U14e's
  band-identity 5f6d holds the compact radial, poverty is radial not
  symmetry (pur payload: exact closure + exact ⟨L²⟩ at +923.1), and
  the class is the normal CASCI shape campaign-wide (C/Si/Sc rows at
  +180…+615 with exact certificates); fixed by the README ΔE rule
  (the partition's mean-field price) + the D-window band clause.
  (2) the reference-echo roots (ΔE ≡ 0: the reference determinant's
  singles Brillouin-dead, the CI reads the mean field back) are
  presented honestly — Ti-unc's (2,5) rows carry them, Sc-unc has
  none (d¹ lacks the same-l double channel); noted in the twin's
  element file.  (3) the raw-fold ⟨L²⟩ mixing is the declared
  (avg, pol)-ruler anisotropy — U14e's 41.9649 is the honest ground
  root (⁵I on all three fields, no root swap; pur exact 42.0000);
  noted in U14e's element file.  Candidate left open below (the
  renderer's note-class).
- **The folded-field campaign identity** (ticket 34's open
  decision): the substrate tag in the key grammar, the free-route
  dominance declaration, the table re-anchoring — user-visible
  before any table carries folded rows.
- **The pinned lane payload seam — RESOLVED 2026-11-18** (7e4ab14):
  the renderer's row gate reads all three payload shapes (the lane
  payloads' ``total`` beside the free cells' ``e_ground``/
  ``energies_pyscf``); the stored pinned legs render with the
  bookkeeping-form disclosure.  Ti/Sc 4→0 pending each, the twins
  54→51, U14e 26→24.
- **test_symm_tm's SA legs** re-scope to store-anchored verification
  (the AGENTS Tests boundary, user-directed 2026-11-05): fresh CI
  legs stay, the SA verdicts load the committed route-batch stores —
  lands at the file's next touch.
- **The char-selected fold seam** (the m-complete 6p promotion,
  CAS(10,15)): the fold-tag grammar cannot express multi-band folds
  and `closed_electrons` assumes one band per shell; the fix is
  character-selected fold spans with reference-occupation counts —
  a real seam extension when wanted.
- **Review declinations, ticket-40 wave** (2026-11-17): (a)
  `SeriesSpec.seed` stays `str | None` + runtime refusal (NOT
  `Literal["self"]`) — consistency with the declaration-style
  validation (`ReferenceSpec.accel`'s own pattern), the battery pins
  the vocabulary; (b) the `_cap` log-line closure stays nested in
  `_window_columns` — display shaping at its one call site, a module
  level name would widen the interface for a print; (c) the probes'
  pick twins (`probe_u14e_mom.seed_columns`,
  `probe_u14e_series.seed_density`) stay un-retired to the
  production `mom_self_seed` — committed historical instruments whose
  logs are the provenance of what ran (their duplicate BLOCK, a
  real bug, was fixed); (d) the seeded-refusal's surviving
  `<chk>.<seed_basis>.source.chk` cache is CORRECT, not an
  asymmetry to fix — the source is the small-basis run certified
  independently at its own creation (its own post-check); the
  target's refusal does not impugn it, the cache convention applies.

### Disposition candidates (the user's future call)

- **The PT2 wave's declared-out classes** (gated recoveries; each
  carries its measured reason in ticket 45's git history): the uncs'
  lane legs (the lane substrate's basin steer diverges on today's
  box); Ti's 3d4d@ar / 3d4s4p4d@ar lanes (the same-l distilled wall
  — `window_distilled`'s future seam — plus schema-stale stores) and
  the semi lanes (the Löwdin needle on fresh construction); Sc's
  semi lanes (the same Löwdin class) + 5 substrate drifts (the regen
  loop, gated); U14e's 11 refusals (the raw-AO windows' re-solve
  does not reconstruct; the 7s SA cells past e_noise; the lane legs
  on real reconstruction failures).

- **The U14e prior optimizer follow-on — EXECUTED 2026-11-18**
  (bb39a70 + the close commit; the prior-anchor-first scenario's
  dated claim re-tested by `--full`, both prior windows, 8 legs
  each): the escape class is UNIVERSAL — the free single-root
  CASSCF descends ~190 mEh below its CI anchor BY BREAKING ⟨L²⟩
  (65.0–65.2, [H,Lz] 5e+00–1.3e+01, cycle-capped 1e-5/200); SA
  holds the cure: −42.0 m at ⟨L²⟩ 72.0000 (4e-7), commutators
  1e-4, converged, on BOTH priors; the frozen pur class certifies
  machine-zero exact ⁵L under optimization (1e-13/1e-14, both
  objectives, both windows).  THE DEEP FINDING: the optimizer's
  landscape is window-INDEPENDENT — @r1 and @s1 free SA read
  −50.4238/−1.1418 identically (the escape floors coincide at
  −50.4298 too): the 5.6 mEh CI-level @s1 advantage was a
  span-capture artifact of the frozen-orbital read; the relaxation
  channel dominates.  The basis decision's evidence is complete:
  both priors optimization-stable under SA, neither prior's
  identity surviving to the optimizer level.  The twins' adoption
  decision carries this evidence (the same un-gating).
- **The swap-basin optimizer class re-run gauge-anchored** (the
  Stage-2 leftover): if the campaign wants in-window solutions on the
  twins, the optimizer legs re-run seeded -- rides the Stage-3
  adoption decision that un-gates `run_casscf`.
- **`_COMM_CAP` (cells) = 6000** skips the dense [H,L] meters past
  6000 dets (the U14e (6,13) sector's em-dashes): this box could
  afford the raise + an ~80-min targeted re-run of the 24 cached
  legs if the certification is ever wanted there.
- The U14e off-basin `moc pur` starts (the multi-solution class,
  honest meters in the payloads) + the `d14p sa moc_frozen pol`
  27.7-mEh reconstruction outlier (the frozen-moc gauge's
  analysis-ruler artifact class).
- The (4,19)/(4,21)-class lattice rows: the ndet² analysis ceiling —
  the same class as `_COMM_CAP`.
- **A data-derived note-class for the poverty rows** (the CASCI
  review's leftover, 2026-11-17): C/Si's D-window CASCI rows read
  +221…+615 mEh with exact closure certificates but an EMPTY note
  column — the poverty signal is invisible there.  A renderer
  note-class (magnitude-derived, keyed on CASCI + ΔE threshold) would
  surface it; the notes are currently identity-derived only.  Rides
  `build_methods_tables.py`'s next touch; the user's call.
- **The Stage-3 adoption decision** (the user's call, the measured
  basis in hand 2026-11-16; RE-SCOPED PROBE-FIRST 2026-11-18 — the
  U14e optimizer wave's evidence: the escape class is universal and
  the optimizer landscape is window-independent, so the twins' full
  51+51-row grids would mostly re-measure a characterized class):
  run a scoped probe batch first (@s2 pur, ~4–6 legs per twin:
  casscf/sa × maskless/frozen), then scale or decline the full grid
  ON THE EVIDENCE.  `subspace_rohf`'s certificate says the
  @s2 windows are MEAN-FIELD-SOUND everywhere (the confinement costs
  1.3–2.2 mEh at the ROHF level, every span, every filling, both
  twins), @s1 is filling-dependent (3.3–38 mEh), @r1 is
  head-filling-only (+1.4/+1.9 for its own head, up to +199 mEh for
  the d-richer members).  Un-gating `run_casscf` on the @s2 windows
  is certified mean-field ground; the queued CASCI review EXECUTED
  2026-11-17 (no analysis bugs found) — the execution waits on the
  probe batch's go.  The swap-basin re-run (below) bundles with this
  decision (the same un-gating).

### Watchlist

- The latent refusal on `--fresh`: the incumbent distilled FRAC
  cells predate the distillation declaration — their caches mask the
  construction refusal.
- The adjudication watch: the pinned solver's roots are proposals
  our exact H corrects (self-verification, not solver independence)
  — a shared-convention error is covered by the dense commutator
  certificates + the ion anchor; revisit if a block class appears
  where adjudication is unaffordable AND the solver skips.
- Revisit conditions: the lane-family extraction from cells (if the
  radial legs grow it); the m-replication/assembly duplication
  (windows vs basis — extract at the next legitimate change, bitwise
  A/B); the per-atom prepared-basis bundle (a second consumer or a
  fifth field); the prior_bands-defaults rule hoist from cells (when
  next touched); the projected-Gram alternative (if the per-radial
  spectrum shows truncated-ranking distortions); the D2h determinism
  `--verify` mode (a second regen cycle); the D2h mol shared builder
  (a third frozen reference); `run_frac_occ(dm0=)`'s disposition —
  kept as the reference family's symmetric explicit-seed form, its
  one consumer was the 2d measurement attempt, revisit if no
  frac-bearing atom appears by the next reference-layer wave.

### Settled (standing facts, one line each; git history the archive)

- **The U14e band-identity poverty is the compact radial, and the
  prior cures it** (tickets 39 + 40, 2026-11-17): the incumbent
  10f+13d window captures only 82-83% of the head's f-valence density
  (the full 9f-11f family 96-97%; 13d 92% vs 98.5%); the @r1 prior
  recovers 84% of the +923 mEh poverty at CASCI (+147) with the exact
  ⁵L on pur, and the @s1 series average is the BEST prior
  (+141.3/+143.7 — ~5.6 mEh below @r1, ~85% cured, exact ⁵L on pur).
  The census's member refusal was the WEIGHT HOOK's vocabulary, not
  the landscape: the character-only premium ties ~12 pure d/f valence
  virtuals + the tight tail at premium ≈ 2.000 (a combinatorial tie
  broken by column index); the self-seeded MOM lands both members
  (f⁴ −50.31571501 l3:4.000; f²d² −50.34788201 under two
  accelerators — plain newton is draw-sensitive there, diis+newton
  lands the anchor in three fresh draws; on f⁴ diis+newton RELAXES
  INTO THE HEAD).  The seed pick rule is ε-first within the majority
  floor, never weight-ranked (the tight-primitive trap).  @s2 is
  construction-refused: the ROHF 5f radial is configuration-inert —
  the series average carries no second radial (k=2 eigenvalue 5.9e-5
  relative), the d shell never asked (k uniform by the tag grammar).
  The basis carries NO near-dependence tail (minima 2.2e-6 s … 7.8e-3
  f; the conditioning cut NO-GO by the twins' precedent; the f block
  is the BEST-conditioned, κ 428).  `run_rohf` deletes a refused
  draw's chkfile (the gate's or the post-check's) — a refused
  draw never poisons the trusted cache.

- **The U14e campaign COMPLETE** (tickets 36A + 38): the free-mask
  SA legs certify the ground term END-TO-END — ⁵L, ⟨L²⟩ = 72.00
  exact, closure + commutators at machine zero, bookkeeping ~1e-3
  mEh; single-root legs break symmetry as declared; the input CI
  grounds are ⁵I/⁵F (the ~1 Eh gain IS the ⁵L relaxation); best
  energies (6,13) −50.456 < (4,15) −50.434 < (4,12) −50.431; the
  retired tier windows' −50.386 validates the promotion; 8
  unconverged legs ACCEPTED (the user, 2026-11-16 — the no-free-mask
  crawl; the physics line rests on the certified legs).
- **The honest prior table** (Sc/Ti/C/Si, the @r/@s windows): NO
  universal k=1 winner (atom-dependent, ~11–100 mEh spread) — the
  SA-over-ladder optimization is the only road to a dominant k=1;
  @s2's dominance is variational headroom at double budget (every
  comparison row carries the CAS shape); C's @s2 REFUSES (the
  rank-1 series average) while Si's stays (a real second radial);
  the committed probe logs are the flipped-metric era — the module's
  numbers are the honest successors.
- **The uncontracted-twins reference disease + the seeded repair**:
  fresh SCF paths diverge on the unc bases — the CHARACTER HOOK, not
  the SCF (metric-amplified virtuals carrying span weights ~1 +
  shell-counter labels selecting AO slivers); the seeded-MOM repair
  (the projected 6-31G source, source-vocabulary hooks, the
  character-of-occupied post-check) lands the declared anchors
  EXACTLY — the product path since Stage 1c; both twins' trusted
  records verified to machine epsilon.
- **Stage 2 completed the twins' prior axis** (2026-11-16): the
  members' recipe is `diis+newton` — Newton alone crosses into the
  deeper sibling basin (the d3→d2s1 fall, −128 mEh, reproduced
  exactly), plain DIIS holds the configuration but never converges;
  `projected_source` must thread the declaring ref's accelerator
  (the member sources are the DIIS-oscillation class). The @s legs
  ran 36/36 per twin, zero refusals; the contracted prior verdict
  holds at the basis limit (the series mean loses to the head's own
  block, @s2 buys the double-budget headroom). The UHF anchor row
  re-landed the censused basin through `uhf_seed_density` (the
  pre-repair record was 11.4 mEh shallow). The frac class
  MEASURED-REFUSED on the twins (both accelerators diverge seeded —
  the 43-AO block feeds the amplified subspace); the averaged core
  is their spherical field.
- **Stage 3 measured the mean-field certificate** (2026-11-16):
  `subspace_rohf` confined-valence solves on both twins, 12 rows
  each — the whole-space legs exact to ≤ 2e-6 mEh (the degenerate
  limit reproduces the anchors), @s2's confinement cost 1.3–2.2 mEh
  everywhere (mean-field-sound), @s1 3.3–38 (filling-dependent),
  @r1 1.4–199 (head-only — excellent for its own head, catastrophic
  for the d-richer members); within-@s monotone everywhere, the
  cross-family r1/s1 nesting measured NOT nested (different
  eigenproblems).
- **The subspace-ROHF's five measured defects** (the bench's
  evidence chain, 2026-11-16): the seed density must reach the
  kernel SPIN-SEPARATED (a summed dm is split (d/2, d/2) — the
  open-shell asymmetry lost, the filling off-basin from cycle one);
  the veff shim returns the TRUE per-spin pair through the substrate
  (a symmetric pair from the summed density carries the
  closed-shell exchange formula, Eh-scale wrong on opens); fills
  size COLUMNS not electrons; NO block-index pools (the subspace
  Fock's eigvectors interleave across the family/window layout — the
  4s ranks 9th globally, inside the family's index range); the MOM
  assignment runs on the CLEAN spin-summed Fock with the anchors in
  the clean gauge (the kernel's Roothaan effective Fock is
  MIRROR-DEGENERATE in the low spectrum — the pairs split at 1e-7
  and flip gauge every cycle; identical Fock, the picks flip from
  exact to garbage).
- **The twins' confined_core_orbitals floor**: the closed-core
  Roothaan on the 65-column family span floors at ~5e-10 past its
  own 100-cycle budget (the conditioning class) — the ion
  cross-check runs on the carbon-class span (exact), the twins'
  refusal recorded.
- **The folded-field records are draw-documented single draws** (no
  replay; the post-kernel converged-filling certificate guards
  production now); the +400 mEh folded offset is basis-stable
  construction physics (the missing core-valence self-consistency),
  not span bloat; the certified numbers stand as measurements.
- **The D2h substrate**: NEVER symmetrize a degenerate-gauge seed
  (`symmetrize_orb` amputates the open MO, `symmetrize_space`
  scrambles the pairing); the solver's irrep blocks erase the gauge;
  pinned-cell windows are D2h-born sources only (any future intake
  needs an irrep-purity check at the seam).
- **The basin zoo is census-resolved** (the twins + U14e): the
  deepest basins are the adopted anchors; plain aufbau can land 10+
  mEh high at the right filling — the seed, not the gate, is the
  cure when gate cycling diverges; the Newton path freezes
  occupations (the occ_ok certificate refuses converged wrong
  configurations).

## Examples (worked records in `elements/`)

- `elements/C/`, `elements/Si/` — the p² campaign.
- `elements/Sc/`, `elements/Ti/` — the transition-metal campaign
  (the radial axis, the d² lattice, the prior windows).
- `elements/Sc-unc/`, `elements/Ti-unc/` — the uncontracted twins
  (the folded/lane records, the prior axis, the seeded references).
- `elements/U14e/` — the ECP atom: the certified ⁵L record.
- `elements/U/` — the frozen pre-campaign ECP workdir;
  `elements/ecp/` — the ECP seam's certification workdir.

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
only); its one driver (`no_group_windows`) lazily imports the
assembly layers. No cycles, no upward edges. Division of labor over
windows: declaration NAMES (labels), windows BUILDS (the kinds +
substrate builders), cells BINDS (the start axis picks the builder).
