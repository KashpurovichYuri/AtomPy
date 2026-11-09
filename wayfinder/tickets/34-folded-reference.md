# Ticket 34 — the folded reference + the window characterization matrix (Ti, then Sc)

Opened 2026-11-08 (user-approved round; supersedes nothing -- new
work).  One session for the Ti leg; the Sc port is the closing step.

## Why (the problem, measured)

The unc lanes' totals sit +118.4 mEh (Tiunc) / +38.7 (Scunc) above the
free rows, and the deep dive showed why: the fold's confined core on a
wide span IS THE FREE ION (Ti4+ RHF reproduced to 0.016 mEh; the
confined core -845.190591 vs the neutral pur core -845.051227 =
-139.4 mEh), and the valence CI runs on that ion field with a window
cut from the NEUTRAL field's MOs -- the frozen window cannot respond
(the 2x2: fixed-wavefunction field cost +119.7 mEh, CI re-relaxation
recovery 0.02 mEh).  The contracted bases sit in the other regime
(their fold spans are the core's own contractions: +18.6 core / -17.1
valence = +1.5 net) -- one construction, span-width-driven regime
flip.  The inconsistency is prerequisite work for the term-ID round:
term splittings are correlation differences and a mismatched field
distorts them.  Elements/Scunc/logs/d2h_probe.log carries the probes.

## Step 1 -- the Tiunc basin census, then PROMOTE (user-decided)

Cycle ROHF init guesses on the unc basis at conv_tol 1e-12, spin 2,
ALL at the character-forced 3d2 4s2 filling where forcing applies
(the plain aufbau draws are basin-chaotic -- the zoo: -848.39639 plain
shallow / -848.40648875 the current trusted record / -848.40744343 the
deepest known, found by the D2h-seeded draw, open pair degenerate at
-0.090219, reverse-seed-stable on C1):

    for guess in (default, "atom", "huckel", "minao"):
        run ROHF(seed or guess) -> collect converged E
    plus the SEEDED starts: from elements/Ti/rungs/rohf.chk (contracted,
    cross-basis lstsq+Loewdin -- reproduces -848.40648875 exactly, the
    38393f8 probe) and from the D2h draw at -848.40744343.

ADOPT the deepest converged same-filling basin as the trusted
reference: re-store elements/Tiunc/rungs/rohf.chk, declare
ReferenceSpec(occ_closed=("1s","2s","3s","4s","2p","3p"),
occ_open=("3d",), occ_e_ref=<the adopted E>) on the case (through the
shell layer -- the forcing spans resolve via shell_realization),
re-run the Tiunc free cells (minutes at the compressed sizes),
regenerate the table.  If the census certifies -848.40744343 the
anchor and every dE shift ~1 mEh -- expected and fine.

## Step 2 -- the folded-reference utility (reference.py)

Frozen-core ROHF: the core density FIXED at the confined solution
(the lane's fold construction), the valence electrons (4 for Ti)
self-consistent in its field, the character forcing + basin
discipline as every runner.  One home; the substrate signature in,
(mol, mf) out.  Tests: (a) the frozen core block is BIT-UNCHANGED by
the run; (b) the folded total >= the full ROHF total (variational,
both same filling); (c) the recovery measurement: the folded valence
re-contracts on the ion field -- predicted to recover ~+120 mEh of
the lane-vs-free gap (a PREDICTION to measure, not assume; if it
does not reconcile, the decomposition machinery says why).

## Step 3 -- the characterization matrix (NO winner; MO default)

On ONE fold (3d4s@ar) and ONE field (the folded reference), run the
window variants and CERTIFY them side by side -- E, <L2>, ptp,
closure ([P,Lz]/[P,Lperp]), window purity, capture:

(i)   AO-compressed k=1, re-cut from the FOLDED reference's MOs (the
      math unchanged; the signal field is the change), report gains
      the PER-RADIAL CAPTURE SPECTRUM (see below);
(ii)  MO-PICKED m-complete -- THE DEFAULT (user-decided): the whole
      reference MOs whose set covers every m-member (5 d partners +
      the 4s), the span closed by set construction, window_closure
      verifying; the core projection owns the core tails; the
      valence-family tails reported by the purity meters;
(iii) NO-informed: one CI pass -> the term-averaged density ->
      natural.occupation_groups -> the correlation-weighted radials
      (the m-set-safe machinery exists).

The support continuum: (i) and (ii) are the ends of the
allowed-support dial (AO combinations restricted to the family = (i);
unrestricted AO combinations near MOs collapse to (ii); the middle =
wider-span radials via a declared support set -- only if the numbers
ask).  The k decision is data-driven: the per-radial spectrum + the
NO occupations (k=2 for the 4s stays one declaration away,
CAS(4,7)); the m-symmetric radial SCAN (scans.radial_scan,
optimizer-free) is the cheap oracle for whether the second radial
direction is energetically live before paying for it.

The MO default is DECLARED after its meters pass -- if a meter fails,
the default stays AO and this ticket records why.  The AO class
remains in the toolbox and the tables: it is the symmetry-preserving
class (band truncation is m-complete by construction; the MO echo's
<L2>=0.005 failure was the m-incomplete SET, not MOs per se).

## Step 4 -- the lane re-runs + the disclosure

The lane legs re-run on the default variant + the folded substrate
(the stored 3F/2D records stand as the constructions that validated
the multi-column fold; the field caveat documented).  The 3F/2D
EXACTNESS re-verified on the new field.  The totals reconciled; the
tables carry the ZERO-DISCLOSURE note (the lane's fold zero is the
confined core = the free ion on a wide span; the free rows' zero is
the neutral core).

## Step 5 -- the Sc port (T-A-prime)

Scunc needs no census (its trusted record is already the deeper of
its two known basins).  The folded reference + the matrix + the lane
re-run, same acceptance.

## Riding along (same commits)

- The ANISOTROPY COLUMN: render core_anisotropy in the tables (the
  lane rows show the machine floor -- itself the certificate).
- The PER-RADIAL CAPTURE SPECTRUM in window_compressed's report:
  p_ij = |e_j^dag S c_i|^2 against each S-orthonormal radial -- the
  k-decision sensor (bounded fractions; rank by G, report by
  projectors).
- The DOCSTRING REVIEW (review, not append -- user-directed): core.py
  (the m-scalar theorem -- the span Loewdin and the closed-shell Fock
  both block-diagonal with identical m-blocks -> replicated
  eigenvectors; the RHF-for-ion framing, brief); operators.py
  (spherical_core: the field defect is quadrupole-like = totally
  symmetric in D2h -- irreps cannot see it; only the SO(3) l=0
  projection isolates it); windows.py (the Gram is the SELECTION
  object -- member-summed ranking, S-metric necessity; the member
  chains); casscf.py (the state_average_mix_ sentence); the twins'
  TABLE preambles once the numbers are final.
- map.md: the projected-Gram DECLINE record (revisit condition: if
  the per-radial spectrum shows ranking distortions).

## Acceptance

The exact term certificates hold on the consistent construction;
the field-window consistency measured (the recovery number); the
captures/spectra/purity recorded; the contracted comparison DIAGNOSTIC
ONLY (evidence, never the anchor -- user-directed: the picking
procedure is judged by our own meters).

## Operational

One regen at a time (the U14e route batch holds its store until it
exits; the watcher completes that table automatically -- do not
interfere).  The fast tier is the per-commit gate.  NO /tmp paths (an
earlier slip): probe chks are in-repo *.tmp beside their targets,
removed after.  Tables presented pipe-checked.  U14e work (T-C) stays
GATED on this ticket + the Sc port.

## After this ticket

T-B the contraction-recipe review (provenance lookup for PySCF's
Ti/Sc ccpvdz; our radials vs the authors' contraction coefficients --
verified structure: Ti (20s16p8d2f) general contractions, s = 5
columns sharing 19 primitives + 1 free, nao 43); T-C U14e (gated);
T-D docs /code-review (anti-inflation); T-E the term-ID round (NIST
as the acceptance test).  Stop after T-E's review -- no physics waves
before it.

## The first round, scoped (2026-11-09, user)

**Ti and Sc fully first (contracted verification + the unc re-anchor),
user review, then U14e.**  The U14e route batch (the CAS(4,31)
lane's SA pur leg) was stopped by decision -- CAS(4,31) is too big
for now; the pol legs (casci/casscf/sa-casscf) are banked in the
store; the future U14e round reviews the DISTILLED 5f6d
(k=1 per shell = CAS(4,12), the m-complete minimal window) and the
promoted ladder (10f13d7s at (5,1) = CAS(6,19)) -- pick then.

The cells:
- Ti/Sc CONTRACTED: the spine/diagnostics CASCI legs re-run with
  `shell_radials = {"3d": 1, "4s": 1}` declared (the distilled
  construction) -- VERIFICATION legs: the 1:1 identity predicts
  span-identical energies vs the certified stores; a discrepancy
  stops the round.
- Tiunc/Scunc: the re-anchor -- the same structure as their
  contracted siblings (3d, 3d4s at the case sectors, casci,
  ao/mo/moc starts, raw/avg/pur folds) + the `3d4s@ar` lane leg on
  the NEUTRAL field (the offset structural ~+124 mEh, disclosed;
  the value = the independent machine-exact term certification).
  SKIPPED: the raw full-family AO windows (the monster spans), the
  optimizer class (gated this wave), the frac fold (unmeasured).
  The mo start STAYS (the m-delocalized control arm).
- C/Si: untouched (the 1:1 identity covers them; their certified
  tables stand).
- Consistency readout: one construction both sides -- term
  identities, splittings, energies (unc <= contracted), captures.

## The folded review (2026-11-09, step 3): the corrected interpretation,
## second revision -- the offset is basis-STABLE construction physics

The contracted control experiment (elements/Ti/logs/folded_review.log)
REFUTES the narrow-span prediction and revises the +400 story AGAIN:

- Contracted Ti: folded offset +334.4 mEh; contracted Sc: +136.9 mEh
  -- vs the unc twins' +399.7/+150.9.  The offset is BASIS-STABLE
  (the drift ≈ +65/+14 mEh is the wide-span confined core's
  over-relaxation, a minor term).
- The confined-vs-SCF cores on CONTRACTED bases are IDENTICAL
  (singular values 0.9999+, core dE +16.9/+10.3 mEh) -- the confined-
  core mismatch is NOT the cause.  The solve replays cleanly on both
  (conv, correct fillings) -- the machinery certificate delivered.
- The offset is therefore the folded CONSTRUCTION's own frozen-core
  SCF penalty: E_folded = E_confined + min{valence in the complement,
  frozen field} misses the core-valence self-consistency, and that
  cost is ~0.15-0.4 Eh for an [Ar]-core TM atom at the SCF level --
  the textbook frozen-core scale, partition-dependent, not
  basis-dependent.  The user's convergence intuition is CONFIRMED by
  measurement: richer basis, same offset.

**Consequences:**
- The folded zero is structurally ~0.3 Eh above the true minimum --
  the demotion (draw-documented study; tables on the TRUE anchors)
  stands, now with the measured reason.
- The first-round IMPLICATION for the lane: the unc lane's
  neutral-field offset (~+124 mEh) is ALSO structural (the lane's
  fold stays full-family, so its confined core ≠ the SCF core on the
  wide span regardless of window distillation).  The lane legs run on
  the neutral field with the offset DISCLOSED -- the lane's value is
  the independent machine-exact term certification, not the totals.

## The folded-field draw-chaos finding (2026-11-09)

The tail's window legs exposed it (the probe diverged where the
certified round converged): the certified folded fields do NOT
replay -- same chks (mtime predates the certified commits), the
exact old module re-imported, the same frozen env, both thread
counts.  Full replay matrix: elements/Tiunc/logs/folded_repro.log.

- The certified Ti/Sc folded numbers (the matrix logs, the lane
  recoveries 5.5e-13/0.085 mEh, the zero disclosures) stand as
  measurements of their draws -- not replayable artifacts.
- run_folded's "deterministic given its inputs" claim was falsified
  on the twins; the docstrings/AGENTS row state the honest form.
- Code forced by the finding (0df7491): the post-kernel converged-
  filling certificate (the certified-era code silently returned
  wrong-fill fields on replay -- the certificate catches all); the
  pick rule DECLARED per atom-gauge (picks='gates' | 'three-level' --
  each form certified on one twin, each fails the other by a
  recorded mechanism); the gates' dominance is a declared parameter.
- The pick-rule history's "each failure was a physics fact" record
  gains: the rules are ALSO draw-fragile -- on today's replays both
  forms abandon the projected basin on both twins.
- BLOCKED on this: the MO-subspace leg, the scan-oracle leg, the
  nested-basis Ti/Sc integration legs, the folded-field rows' entry
  into any store (their addresses exist -- the field tag -- but no
  replayable payload can be minted today).
- NOT blocked, delivered: the schema commit (e6123fe), the
  nested-basis carbon keeper (05af7b9 -- the contracted basis
  certified as its own split's realization, convention stack exact).
- FIX DIRECTIONS (user-gated): (a) the folded reference as a STORED
  artifact (a folded-field chk beside the ROHF chk -- the draw
  becomes the trusted input, the reproducibility unit moves to the
  store); (b) a folded-basin census (multi-start + the certificates,
  the ROHF census pattern); (c) re-scope the folded-field claims to
  draw-documented results.

## The wave review (2026-11-09) -- the two-axis code-review triage

The review of `8bb62e9...HEAD` (the code-review skill, Standards +
Spec sub-agents).  Verdict: the physics core clean (steps 1/5, tests
(a)/(b), the ³F/²D re-verification delivered as specified; the
ReferenceSpec non-declaration a documented measured outcome).

FIXED (158de72): the AGENTS.md lockstep row (rule 6); the run_folded
docstring generalized (rule 1 amended -- the measured numbers moved
to map/logs); the _occ fail-loud gate-shortfall guard (the silent
short-fill finding); the scf_summary honesty (e1 + e2 = e_val); the
redundant local imports.  map.md tense word.

HELD (accepted, no action): the watcher script + its absolute path
(the operational record); the private `_closed_field` test import
(defensible micro-pin); the docstring-length norm-stretch; the
recovery-prediction deviation (disclosed).

ADDED TO THE TAIL (this ticket's remaining work):
- the projector_weights helper extraction (the einsum idiom, ~8
  call sites across modules -- mechanical, lockstep with tests);
- the radial-scan oracle leg: the k=1 verdict currently rests on
  captures + NO groups alone -- run the scan alongside the
  MO-subspace probe below.

DECLINED (recorded): no full rename migration of the key grammar
(the addresses are store identities; the semantics get documented
instead).

## The decisions (2026-11-09, user)

- **The Sc free-route folded cut: dominance 0.4 APPROVED** (the C1
  gauge's 0.4628 signal declared; measured honest -- L2 = 6 exact).
- **The re-anchor DECLINED**: the tables keep the TRUE ROHF (and
  UHF) anchors -- "we should compare to true ROHF and UHF".  The
  folded-field rows, when they enter the campaign, carry the field
  tag + the offset disclosure (their zero = the confined core = the
  free ion on a wide span) alongside the true-anchored rows.
- **The m-complete MO-subspace leg APPROVED**: the window = the
  top-5-dimensional d-subspace (the member-summed Gram of the
  projected significant-d columns) + the 4s column -- rotation-closed
  by the pur theorem; the meters decide against the AO default at
  the same CAS size.  Attribution note: per-column labels are gauge
  (the weight vector + the natural occupations are the honest
  per-column identity); subspace projector traces are the physical
  certificates.
- **The nested-basis integration test APPROVED** (the user's design):
  the contracted AO space is a 43-dim subspace of the unc 122-dim
  space (identical primitives), so (i) the unc machinery restricted
  to the contracted subspace on the neutral field reproduces the
  certified contracted record -- the strongest machinery anchor;
  (ii) the same restriction on the folded field quantifies the
  contraction cost in the folded regime; (iii) the carbon pinned
  version at unit-tier sizes (both mols already in the tier) is the
  permanent keeper.  Feeds T-B (the authors' contraction = the
  k=1 analog, quantified against ours).
- **T-E sharpened**: the term splittings, folded field vs neutral
  field vs NIST -- the infrastructure now supports both fields
  consistently.

## The amended program (2026-11-09, post-finding; supersedes "The
## remaining execution order")

The user's narrowing after the draw-chaos finding and the +400 mEh
review.  The corrected interpretation first:

**The +400 mEh is NOT expected physics and NOT basis convergence.**
It is the cost of freezing the BARE confined core (the span's own
Roothaan, solved in the ABSENCE of the valence), which over-relaxes
on a rich span: measured decomposition ≈ −139.4 mEh (the core
deeper -- no valence repulsion shapes it) + ≈ +539 mEh (the valence's
strain in the mismatched field).  The PHYSICAL frozen-core cost is
the +1.5 mEh scale, basis-stable -- the contracted basis shows it
because its stiffness keeps both constructions' cores coincident.
Variational truth: E_folded ≥ E_neutral ALWAYS (the folded
configuration is accessible to the neutral minimization); freeze the
SCF's own core and the gap collapses to the mEh relaxation scale on
ANY basis.  The old +124.3 lane-vs-free gap is the same
confined-core-mismatch disease in a milder form.  Both die with the
bloated span: on distilled (narrow) windows the confined core ≈ the
SCF core, and the lane joins the true minimum's zero by construction.

**Locked decisions:**
- The m-complete MO-subspace windows are ENOUGH for the start window
  and the first round (C/Si/Sc/Ti/U14e).  NOs, the basis exporter,
  and the auto-basis-optimization scheme are FUTURE stages.
- The compression is DELETED (user: "bad and unclear" -- top-k of one
  reference's occupied space is not how bases are built).  Sequencing:
  step 1 certifies the replacement, deletion follows immediately.
- The folded field: demoted to a draw-documented study.  The lane
  runs on the neutral field.  Stage 2/3 of the basis pipeline are
  specced, not built: the CONTRACT is "a distillation delivers,
  per chemical shell, k radial vectors over that shell's primitives,
  m-complete by replication" -- sources pluggable (SCF-MO subspace,
  NO groups, the authors' coefficients); stage 2 = the exporter
  (inverse make_uncontracted, certified by the nested-basis seam);
  stage 3 = the on-basis road (own zero).

**Execution order:**
1. The m-complete MO-subspace window (neutral field): the Ti probe --
   meters (window_closure), the ³F readout vs the AO-compressed
   default at the same CAS size, the radial-scan oracle alongside;
   unit-tier tests (carbon fixture).  Home in windows.py only if the
   deletion test admits it (probe first).
2. Delete the compression: window_compressed, compressed_ncas,
   _TM_UNC_RADIALS, the capture plumbing, their tests; re-anchor the
   twins; amend AGENTS.md fact 26 + the windows.py layout row; logs
   and committed records stay as history.
3. The contracted-vs-unc folded review (the control experiment):
   run_folded on CONTRACTED Ti/Sc ONLY -- expected mEh-scale offset,
   replay-clean (the unc side stays draw-documented; it does not
   replay).  Optional diagnostic leg: the SCF-core-frozen variant on
   unc (predicted mEh -- proves the +400 was the confined-core
   choice).
4. Finish the first round: C/Si/Sc/Ti/U14e on the m-complete MO
   windows; the OWED verification -- the lane legs on the distilled
   windows vs the free rows, residual at mEh scale expected.
5. Deferred: NOs, the exporter, stage 3, the auto-optimization
   scheme -- discussed after the round.
6. The tail/queue: projector_weights extraction, the anisotropy
   column, the docstring review (the per-radial spectrum died with
   the compression), T-C gated on the U14e batch, T-B, T-D, T-E
   (splittings vs NIST).
