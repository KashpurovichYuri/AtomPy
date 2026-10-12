# Operator construction: first quantization

## Question

How should the many-electron operators (H, Lz, L², S², N) be represented?

## Resolution (closed — user decision)

**First quantization only**: full matrices in the determinant basis, built by
explicit Slater–Condon rules from one-body matrices and two-electron integral
tensors. No second quantization: no creation/annihilation algebra, no
transition-RDM machinery.

- `detbasis.py` (numpy-only) implements the rules; generic over complex,
  non-symmetric two-body kernels (fully antisymmetrized 4-term form).
- `L² = Σᵢ l²(i) + 2 Σ_{i<j} lⱼ·lᵢ` with composite integrals
  `(pq|rs)^LL = Σ_a l_a[p,r]·l_a[q,s]`; the same composite primitive builds S².
- Validation: eigenvalues vs `fci.direct_spin1.kernel` (independent solver,
  same integrals) and vs the hand check `E(³P det) = 2 h_2p + J − K`;
  textbook 2p² spectrum {0×1, 2×9, 6×5} in the Ms=0 sector.
