# Ticket 16: Deep understanding of C/Si/Sc orbital effects (Phases 0-2)

Status: open (plan: /home/kashpur/.opencode/plan/deep-understanding-pre-ti.md)

## Question

Tickets 14/15 produced findings (facts 16/17) that we accept on
measurement but want to *understand* and *show* step by step, with
built-in/simple methods first: why radial mixing carries the whole
CASSCF gain on p2 atoms, why the m-breaking is energetically flat, what
exactly each SA deviation is made of, and whether an independent window
constructor (AVAS) agrees with ours. No Ti until this lands (user
decision 2026-10-15).

## Resolution (to implement)

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

## Acceptance

RUNGS.md/GLOSSARY readable without prior context; scans reproduce
E(casci_avg) at th=0 exactly and bracket the CASSCF gains; the flatness
is measured, not inferred; the decomposition table covers every SA rung
on C/Si/Sc; AVAS comparison done with numbers; facts/tests land with
each phase; 70-test suite stays green.
