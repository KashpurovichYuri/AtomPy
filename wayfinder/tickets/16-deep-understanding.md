# Ticket 16: Deep understanding of C/Si/Sc orbital effects (Phases 0-2)

Status: RESOLVED (Phase 0: 13a1b15; Phase 1: b1c3c52; Phase 2: fact 19 —
`window_closure` / `avas_window` / `span_overlap`, StageFacts + ladder
carry the one-particle columns, `scripts/sa_decomposition.py` artifacts
in `results/*_sa_decomposition.md`). The m-scalar-constraint rung is
DECIDED AGAINST (fact 19e: the member gauge is exactly flat on p² and
sa_radial already realizes the m-symmetric optimum reproducibly).

## Question

Tickets 14/15 produced findings (facts 16/17) that we accept on
measurement but want to *understand* and *show* step by step, with
built-in/simple methods first: why radial mixing carries the whole
CASSCF gain on p2 atoms, why the m-breaking is energetically flat, what
exactly each SA deviation is made of, and whether an independent window
constructor (AVAS) agrees with ours. No Ti until this lands (user
decision 2026-10-15).

## Resolution (measured)

Phase 0 (13a1b15): RUNGS.md + GLOSSARY.md plain-language rewrite.

Phase 1 (b1c3c52, fact 18): frozen scans -- C/Si gain is a
one-parameter m-symmetric radial rotation at exact symmetry (theta* =
0.61/0.85 rad, -0.2508/-0.1267 Eh); single-member tilts cost EXACTLY 0
(M_L degeneracy dodge); [H,Lz] blind spots (m=0 member, equal +-pairs);
Sc inverts (symmetric gains nothing, member-asymmetric descends).

Phase 2 (fact 19): (a) one-particle closure meter ||[P,Lz]||/||[P,Lp]||
-- complete certificate, closes the blind spots; carried by every
ladder row; exposes C's radial-rung member tilt ([P,Lp] ~ 1) invisible
to sector commutators. (b) SA-deviation decomposition: core-fold witness
<= 1e-13 everywhere (deviations are 100% window); free SA HEALS below
the input baseline on Si (4e-6 vs 3.4e-3), only constrained sa_radial
sits at it; Sc sa = swap basin (sv 0.00) / sa_radial = baseline /
sa_pyscf = collapse; canonical-start baselines huge (0.24-3.7). (c)
span overlap != closure (AVAS-Sc: 0.9855 overlap, 0.5 defects). (d)
AVAS = mean-field-inherited relaxation: C/Si exactly the theta* optimum
(-250.8/-126.7 mEh, near m-symmetric); Sc +80.8 mEh WORSE with
masquerade <L^2> 5.90 -- raw shells strictly dominate. (e) Sc
tilt-follow to delta=0.35: descent to -2.7/-3.6 mEh with L drifting
2.00 -> 1.84 (the 2D->2P relabeling anatomy); an Lz-BLIND descent
family caught only by [H,L^2].

Next: ticket 13 (= Phase 3, the closed core span) on the strength of
fact 17 + the measured input baselines.

## Resolution (the original spec, superseded by "measured" above)

Phase 0 -- plain-language docs:
- RUNGS.md: every ladder rung explained in plain physics (what is held
  fixed, what is released, which question it answers) + the report
  columns and the divergence guard.
- GLOSSARY.md rewritten in plain language (window, rung, member, gauge,
  flat direction, escape basin, input-window baseline, fold).

Phase 1 -- the radial mechanism, step by step (C first, then Si, light
Sc), NO optimizer anywhere:
1. scans.py: scan windows cos(th)*A_m + sin(th)*B_m over two radial
   families (S-normalized columns; m-members orthogonal exactly).
2. Frozen m-SYMMETRIC scan: common th for all members; E(th) for the
   ground state split into one-body and two-body parts; compare E(0)
   (= casci_avg, must match to machine precision), the optimum vs the
   CASSCF gains (-251 mEh C, -127 mEh Si, -12 mEh Sc vs relax).
3. Frozen PER-MEMBER scan around the optimum: measure the flatness of
   the member direction directly (expected: plateau to ~1e-5 Eh with a
   quadrupole-scale ripple) -- the step-by-step proof of
   commutator = O(eps), energy = O(eps^2).
4. Cross-check the mask rungs against the scan-predicted composition.

Phase 2 -- SA-deviation decomposition + AVAS comparison:
- Decomposition table per atom: window-identity change / member gauge /
  input-window baseline / fold convention, including the measured row
  "free SA in the complete {3d,4s,4p} window does NOT escape" and the
  L~0 masquerade caveat (small <L^2> in any window -- only commutators
  measure closure).
- AVAS probe (C + Sc): window from AVAS orbitals vs raw shell AOs, same
  CASCI fold, compare energy + commutators (expected: span-close, not
  l-pure; verify, don't assume). Review avas.py + arXiv:1701.07862.
- Optional m-scalar-constraint probe -- decide after Phase 1 quantifies
  what the per-member gauge buys (expected: nothing).

## Acceptance (met)

RUNGS.md/GLOSSARY readable without prior context ✓; scans reproduce
E(casci_avg) at th=0 exactly and bracket the CASSCF gains ✓ (fact 18);
the flatness is measured, not inferred ✓ (exactly 0 to 1e-12); the
decomposition table covers every SA rung on C/Si/Sc ✓ (artifacts in
results/*_sa_decomposition.md); AVAS comparison done with numbers ✓
(fact 19d, verified rather than assumed -- the same-span guess was
refuted); facts/tests land with each phase ✓; suite green throughout
(78 tests).
