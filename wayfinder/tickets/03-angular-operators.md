# Angular momentum operators in the AO basis

## Question

How to build exact one-electron Lx/Ly/Lz/L² matrices for an atom?

## Resolution (closed)

Two independent constructions, tested to agree:

1. **libcint (primary)**: `int1e_irp` (9 components) decodes as
   `T[3a+b] = <mu| r_a d_b |nu>`; `Lz = −i(T1−T3)`, `Lx = −i(T5−T7)`,
   `Ly = −i(T6−T2)`. Verified: Hermitian, pure imaginary, `[Lx,Ly] = iLz`,
   p-shell eigenvalues {−1,0,1}, Condon–Shortley consistent.
2. **Analytic (cross-check)**: per-shell real↔complex harmonic transform U
   from exact tesseral algebra (PySCF p-ordering `(px,py,pz)` for l=1,
   m-ascending otherwise); `L^a` block `= U† gen_a U`; matrix elements from
   overlaps within the shell. `L² = Σ_a l_a l_a` computed in the
   symmetrically orthonormalized AO basis (matrix products are exact only in
   an orthonormal basis) and transformed back with `S^{±1/2}`.

Limitation (documented): single atom at the origin, spherical basis.
