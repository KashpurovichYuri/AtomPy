# AtomTerm

Atomic term-symbol (L-multiplicity) analysis on top of [PySCF](https://pyscf.org).

The pipeline runs ROHF → CASCI → CASSCF → SA-CASSCF (all checkpointed),
then classifies every obtained state by **exact** angular momentum
operators `Lz`, `L²`, `S²` built as plain matrices in the
Slater-determinant basis (first quantization — no second quantization),
and quantifies rotational-symmetry preservation by `‖[Ĥ, L̂z]‖_F` and
`‖[Ĥ, L̂²]‖_F`.

Atoms covered: **C** (2p²), **Si** (3p²), **Sc** (3d¹4s² — first doublet,
first transition metal).

## Quick start

```bash
conda activate pyscf          # PySCF 2.13 environment
pytest                        # 57 tests, ~2 min (chkfile caching)
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
exactly rotation-closed.

## The symmetry story (what the commutators measure)

- **Closed shell window (raw shell AOs)**: commutators vanish to ~1e-13
  and the full free-ion term table comes out exactly (Sc CAS(3,6) over
  {3d, 4s}: all 90 states of 3d¹4s² + 3d²4s¹ + 3d³, ground ²D).
- **Single-root CASSCF**: measurably breaks symmetry. On C/Si the window
  stays l-pure (m-incompleteness); on Sc it escapes into p and collapses
  to a dominantly-L=0 basin — orthonormal and variational, a real
  multi-solution of the TM landscape, not a bug.
- **SA-CASSCF over the term's M_L components**: restores ⟨L²⟩ on C/Si to
  machine precision. On Sc the shell-start SA restores ⟨L²⟩ but not
  [H, Lz]: the optimizer swaps the window's 4s for a core 3p (m-incomplete
  window) — which motivates the constrained-optimization axis
  (same-(l,m) radial-only rotations) in `wayfinder/tickets/`.

## Design

First-quantized operator algebra: one-body matrices + two-electron integral
tensors in, full `(ndet × ndet)` matrices out (explicit Slater–Condon rules;
an independent second-quantized construction is pinned elementwise-equal by
tests). `L² = Σᵢ l²(i) + 2Σ_{i<j} lⱼ·lᵢ` with composite integrals
`(pq|rs)^{LL} = Σ_a l_a[p,r]·l_a[q,s]`; same primitive builds `S²`.

Two independent L constructions cross-check each other (libcint `int1e_irp`
decode vs exact tesseral algebra). Half-integer S is handled on the
half-grid (doublets/quartets label correctly).

See `AGENTS.md` for module map, conventions, and hard-won implementation
facts; `wayfinder/map.md` for the effort map and live tickets.
