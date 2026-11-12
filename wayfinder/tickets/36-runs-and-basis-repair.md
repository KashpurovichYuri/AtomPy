# Ticket 36 — the runs wave + the basis-repair tracks

## The approved batches (B -> C -> D -> E -> A; U14e LAST, the user's
## ordering 2026-11-12) -- APPROVED, NOT YET RUN

- **B (Sc)**: the prior treatments -- `regen Sc --group planned
  --full` (safe: Sc's priors all constructible; no same-l in the
  group).  ~24 cells: @r1/@s1/@s2 x CASSCF/frozen/SA/SA-frozen at
  both spans (+ the NO store-only certificates).
- **C (Ti)**: the same treatments, PER-WINDOW `--full` x6 (+no
  labels) -- NEVER blanket: Ti's planned group also carries the
  same-l windows and materializing those hits the distillation
  refusal.
- **D (twins)**: `--group grid --full` each (the primary optimizer
  legs, the _unc_gates anchor-first re-test), then per-window
  `--full` for 3d4s@r1 / 3d@r1.  The @s windows NEVER queued (the
  fresh-SCF disease -- their repair is the deferred track below).
- **E (C/Si priors)**: DECLARED (cases.py, censused 2026-11-12; the
  tables carry the pending rows, 73 each).  The CI anchors
  (raw/avg/pur + kinds, seconds each), then `--group planned --full`
  for the treatments.
- **A (U14e, LAST)**: the (6,13) CI anchors + BOTH lane families
  (Rn*7s*lane, Rn*lane) + `--full` optimizer legs across ALL THREE
  CAS sections.  The heaviest (17-root SA, pinned solvers/irrep).

Mechanics: one regen at a time, OMP 16, background chains, per-cell
stores; re-render + commit per batch.  Expected physics: local-
minimum flags on free single-root legs where they sit above their
span's SA ground root; SA anchors the member gauge on the twins.

## The deferred tracks (the user's design, converged 2026-11-12)
## -- BENCH NOT STARTED ("later and more accurately")

**Track P -- projection.**  The on-demand 6-31G source (forced-SCF
there is tame -- no dependence on any cached solution); the
projector C~_B = S_B^-1 X C_A (X = intor_cross int1e_ovlp) with the
CONDITIONING TRUNCATION INSIDE (S_B^-1 amplifies the null space of
a near-singular B -- the robust projector solves in the truncated
eigenbasis; P needs C internally whether or not C ever cuts the
working basis).  Uses: (1) the seed (project + re-relax, MOM-gated;
6-31G coarseness absorbed by the re-relax), (2) the density-direct
@s source (the prior needs only k=1-2 radials -- coarse densities
carry that robustly).  Battery: B=A identity -> superset-exactness
-> H/He/C -> the known-answer transfer (the contracted heads ->
the twins' trusted heads) -> the near-dependent stress.

**Track C -- conditioning.**  The spectrum measurement first (S per
l-block on both twins; the design rests on a clean gap); the cut =
canonical orthogonalization expressed as generalized contractions
(a standalone basis object); the synthetic-C bug-and-fix; the
head certificates.  Guaranteed convergence, the truncation floor
measured.  Fallback if seeded P still falls into the null space.

**The dependency the user named**: cut/projection development => the
@s sources exist => the twins' @s legs run (CI anchors, then the
treatments -- the same classes as B/C/D) => _unc_s_block retires.
The bench is DEFERRED until the user calls it; P first, C fallback,
they compose (conditioned basis + projected seeds).

## Then: the ROHF basis exploration

The subspace-ROHF (the mean-field certificate): valence chains
contracted (@r1 first -- the runs' data anointed it), the CORE
EXACT (the cusp lesson), the full-basis Fock with rotations
confined to the span (the run_folded mirror; core.py's span-
Roothaan the closed-shell precedent); DeltaE vs the anchor per
contraction; Sc first.  Its outcome un-gates the two pending
optimizer classes (the prior treatments of B/C/D -- if the basis
decision redirects them -- and the twins' primary legs).
