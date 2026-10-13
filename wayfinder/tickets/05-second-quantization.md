# Second-quantization approach

## Question

Should the many-electron operators also exist in a second-quantized
construction, and how is it validated?

## Resolution (closed)

Yes. `secondquant.py` implements correct fermionic operator application:
`annihilate`/`create` primitives with exact sign tracking, one-body
`Σ h a+ a` and two-body `½ Σ G[p,q,r,s] a+ a+ a_s a_r` strings applied to
occupation vectors (column-wise: M[:, J] = O|J>). The two-body index
arrangement and the ½ prefactor are derived in the module docstring and
pinned elementwise against the first-quantized Slater–Condon matrices.

Validation on the carbon CAS(2,3) tests, three tiers:

1. elementwise identity of both constructions for H, Lz, L², S², N in the
   Ms=1 and Ms=0 sectors (1e-12);
2. core-explicit full space CAS(6,5) (raw integrals, no folding, 50/100
   determinants): eigenvalues vs PySCF `CASCI(ncore=0)` to 1e-8; full-space
   `[H, Lz] = [H, L²] ≈ 0`; N = 6; ground state ³P with ⟨L²⟩ = ⟨S²⟩ = 2;
3. the doubly-occupied-core block of the full-space H equals the folded
   active-space Hamiltonian + explicit core energy (1e-10 against the
   first-quant folded path, 1e-8 against PySCF).

Implementation notes: core columns for the full space are the canonical core
MOs **purified to exact l = 0** (`drivers.purified_s_core`: project onto the
s-AO subspace, drop the SCF p/d tail — measured 1s 1.9e-7 / 2s 6.3e-6 — then
Löwdin-reorthonormalize in the AO metric). Raw s AOs are an equivalent
fallback; unpurified canonical MOs break orthonormality and show up in the
commutators at 1e-4. `build_operators(..., method="slater_condon"|"secondquant")`
dispatches the many-body evaluation; `full_space_integrals` provides
unfolded integrals.
