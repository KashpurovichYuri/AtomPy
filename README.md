# AtomTerm

Atomic term-symbol (L-multiplicity) analysis on top of [PySCF](https://pyscf.org).

For a carbon atom (ROHF/cc-pVDZ, triplet) the pipeline runs ROHF → CASCI →
CASSCF → SA-CASSCF, then classifies every obtained state by **exact** angular
momentum operators `Lz`, `L²`, `S²` built as plain matrices in the
Slater-determinant basis (first quantization — no second quantization), and
quantifies rotational-symmetry preservation by `‖[Ĥ, L̂z]‖_F` and `‖[Ĥ, L̂²]‖_F`.

## Quick start

```bash
conda activate pyscf          # PySCF 2.13 environment
pytest                        # 21 tests, a few seconds
python scripts/run_carbon.py  # writes results/ (chkfiles + report.md/.json)
```

## Example output (carbon, CAS(2,3) = 2p², Ms=0 sector)

```
 #           E (Eh)      <L^2>  L      <S^2>  S  term
 0      -1.10060875   2.000000  1   2.000000  1    3P
 1      -1.10060875   2.000000  1   2.000000  1    3P
 2      -1.10060875   2.000000  1   2.000000  1    3P
 3      -1.01251492   6.000000  2   0.000000  0    1D      (x5)
 8      -0.88037417   0.000000  0   0.000000  0    1S
Lz eigenvalues (M_L content): [-2, -1, -1, 0, 0, 0, 1, 1, 2]
||[H, Lz]||_F  = 1.5e-15      ||[H, L^2]||_F = 3.2e-15
```

Textbook 2p² multiplet (Hund: ³P < ¹D < ¹S) recovered *exactly*; commutators
vanish to machine precision because the active space (the 2p shell AOs) is
exactly rotation-closed. Single-root CASSCF measurably breaks the symmetry
(⟨L²⟩ ≈ 1.67); state-averaging over the three M_L components restores it.

## Design

First-quantized operator algebra: one-body matrices + two-electron integral
tensors in, full `(ndet × ndet)` matrices out (explicit Slater–Condon rules).
`L² = Σᵢ l²(i) + 2Σ_{i<j} lⱼ·lᵢ` with composite integrals
`(pq|rs)^{LL} = Σ_a l_a[p,r]·l_a[q,s]`; same primitive builds `S²`.

See `AGENTS.md` for module map, conventions, and hard-won implementation facts.
