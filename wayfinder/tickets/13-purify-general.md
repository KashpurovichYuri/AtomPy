# Ticket 13: generalized core purification (the closed core span)

Status: open (PROMOTED from optional by fact 17 / ticket 14: an exactly
closed core span removes the ~1e-3 input-window commutator baseline on
Si/Sc, making optimized rungs exactly closed. Plan: Phase 3 of
`~/.opencode/plan/deep-understanding-pre-ti.md`; runs after ticket 16.)

## Question

`purified_s_core` handles s-type cores. How to purify p-core (Si 2p⁶) or
future d-core (Cu 3d¹⁰) orbitals -- and can the purified core span make
the CASSCF input window exactly rotation-closed?

## Resolution (to implement)

`purified_core_orbitals(mol, mf, n_core) -> (C, contamination)`:

1. For each canonical core MO, identify its dominant (n, l) shell = the
   AO shell with the largest metric weight (c† S c over that shell's AOs).
2. Project each MO onto ITS OWN shell's AO subspace (s MOs → s AOs, p MOs →
   p AOs, d MOs → d AOs): drops all foreign-l numerical tails at once.
3. Löwdin-orthonormalize per l-block (shells of different l are exactly
   orthogonal, so per-block Löwdin suffices and preserves l-purity).
4. Returns per-MO contamination `1 − shell weight`.

SHARED-RADIAL refinement (the fact-17 killer): closure of the core span
additionally requires the core's radial vectors to be shared across the
m members of each l (see GLOSSARY "m-completeness"): m-average the
radial content within each closed l-shell before step 3. This is what
build_full_mo/radial_basis then project the window against.

Invariant to assert: m-completeness — for a closed p-shell core the three
projected p MOs must span the whole (n, l) AO subspace (projection is
near-identity since contamination is small; verify rank).

Step-by-step verification BEFORE any re-optimization:
(a) purified core ~= SCF core (density change ~1e-9; fact 9 precedent);
(b) input window orthogonalized against the purified span: commutators
    predicted to drop from ~1e-3 to ~1e-13 (Sc/Si);
(c) re-run relax / casscf_radial / sa_radial on the purified substrate:
    predicted machine-precision closure on C/Si, gauge-only residual on
    Sc.

Validation: same three-tier chain as `purified_s_core` (purity via
`L² c = l(l+1) c`, orthonormality, commutators) — plus for a closed-shell
core atom (e.g. Si²⁺ or Mg) the full-space fixed-core test.
