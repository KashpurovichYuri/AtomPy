# Ticket 13: generalized core purification (the closed core span)

Status: open (PROMOTED from optional by fact 17 / ticket 14: an exactly
closed core span removes the ~1e-3 input-window commutator baseline on
Si/Sc, making optimized rungs exactly closed. Plan: Phase 3 of
`~/.opencode/plan/deep-understanding-pre-ti.md`; runs after ticket 16.)

## Question

`purified_s_core` handles s-type cores. How to purify p-core (Si 2p⁶) or
future d-core (Cu 3d¹⁰) orbitals -- and can the purified core span make
the CASSCF input window exactly rotation-closed?

## Resolution (to implement -- refined & approved 2026-10-15)

`purified_core_orbitals(mol, mf, n_core) -> (C, contamination)`:

1. For each canonical core MO, identify its dominant (n, l) shell = the
   AO shell with the largest S-METRIC PROJECTOR weight (the fact-15c
   lesson: projector weights, never bare quadratic forms).
2. Project each MO onto ITS OWN shell's AO subspace (s MOs → s AOs, p MOs →
   p AOs, ...): drops all foreign-l AND foreign-n same-l tails at once.
3. SHARED RADIAL COMES FREE for closed-shell cores (C 1s2s, Si [Ne], Sc
   [Ar]): the number of core MOs of each closed shell equals its AO
   count, so the projected set spans each FULL (n,l) AO subspace -- and
   a full AO subspace is rotation-closed with shared radial vectors by
   construction. PROVE, don't assume: a rank/m-completeness check per
   closed shell is a pinned test.
4. Orthonormalize with ONE global Loewdin on the projected set: S is
   exactly zero between different l, so l-purity survives and the span
   (direct sum of full AO subspaces) is closure-invariant.
5. Return per-MO contamination `1 − shell weight`.

Verification ladder, each step gated before the next:
(a) purified core ~= SCF core: density change ~1e-9 (fact 9 precedent);
    if bigger, the reference changed -- stop.
(b) THE HEADLINE: build_full_mo / radial_basis windows orthogonalized
    against the purified span, measured with window_closure -- predicted
    Si 3.4e-3 -> ~1e-13, Sc 5.7e-3/2.1e-2 -> ~1e-13 (the projector onto
    a full AO subspace commutes with L, so the window projection leaks
    nothing). If (b) fails the construction is wrong -- do not optimize.
(c) only then re-run relax / casscf_radial / sa_radial (+ sa) on the
    purified substrate (thread c_core= through; build_full_mo /
    constrained_casscf already accept core overrides): predicted
    machine-precision closure on C/Si; Sc drops to the ~1e-13 baseline
    with only the genuine JT gauge residual (fact 18d/19e) left. sa on
    Sc may stay bimodal (optimizer attractor, not baseline artifact) --
    keep branching tests on window weights.
(d) energy effect of the purification itself: rung energies purified vs
    canonical substrate (predicted ~1e-9 Eh) + the decomposition
    before/after table (sa_decomposition.py gains a purified pass).
    max|dE| columns may shift class (fold interacts with the ~1e-9
    different core): expected, document, don't chase.

Phase 4 (post-hoc closure repair) stays CONDITIONAL on (c) leaving a
residual on p2 atoms -- after fact 19b (free SA heals below baseline on
Si) it is expected NOT to be needed.

Validation: same three-tier chain as `purified_s_core` (purity via
`L² c = l(l+1) c`, orthonormality, commutators) + the per-shell rank
check + (b) as a pinned invariant. One session; commit per milestone;
suite green.
