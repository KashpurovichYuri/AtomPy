# Ticket 10 (optional): generalized core purification

Status: open (not blocking 07–09; needed only for core-explicit TM tests)

## Question

`purified_s_core` handles s-type cores. How to purify p-core (Si 2p⁶) or
future d-core (Cu 3d¹⁰) orbitals for core-explicit second-quantized tests?

## Resolution (to implement)

`purified_core_orbitals(mol, mf, n_core) -> (C, contamination)`:

1. For each canonical core MO, identify its dominant (n, l) shell = the
   AO shell with the largest metric weight (c† S c over that shell's AOs).
2. Project each MO onto ITS OWN shell's AO subspace (s MOs → s AOs, p MOs →
   p AOs, d MOs → d AOs): drops all foreign-l numerical tails at once.
3. Löwdin-orthonormalize per l-block (shells of different l are exactly
   orthogonal, so per-block Löwdin suffices and preserves l-purity).
4. Returns per-MO contamination `1 − shell weight`.

Invariant to assert: m-completeness — for a closed p-shell core the three
projected p MOs must span the whole (n, l) AO subspace (projection is
near-identity since contamination is small; verify rank).

Validation: same three-tier chain as `purified_s_core` (purity via
`L² c = l(l+1) c`, orthonormality, commutators) — plus for a closed-shell
core atom (e.g. Si²⁺ or Mg) the full-space fixed-core test.
