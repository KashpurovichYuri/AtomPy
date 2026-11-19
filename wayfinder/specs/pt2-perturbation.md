# SPEC — The perturbation analyzation (`perturb.py`), PySCF-mrpt-primary

The wave's contract, REWRITTEN 2026-11-19 (the replan): the original
premise — "no PT2 implementation exists in-env" — was FALSE. PySCF
2.13.1 ships `pyscf/mrpt` natively (`NEVPT` = strongly contracted
NEVPT2; `dfnevpt2`; **no `caspt2` in 2.13.1** — master-only). The
home-grown ladder (FOIS/EN rungs 0–1, the Dyall rung, the density
layer) is DELETED — git history is the archive; the Dyall rung's
committed numbers were invalid anyway (its virtual-range `eigh` mixed
α/β virtuals through degenerate eigenspaces; the coupling norm
collapsed 0.016 → 5.4e-11 on the carbon fixture). The engine is PySCF's;
we own the analyzation.

## Problem Statement

The campaign's tables certify *static* correlation only: a CAS span,
its exact symmetry, its commutators. But the deeper question is the
**non-saturation of CASSCF itself**: enlarging the active window keeps
buying large, non-converging amounts of energy (the poverty class:
+250.8 mEh on C, +221.5 on Si between the compact and @r1 windows, at
identical determinant counts and exact certificates). Two readings
compete: the window is missing genuinely STATIC correlation
(near-degenerate radial directions that must be variational), or the
CAS is absorbing DYNAMIC correlation — badly, at growing cost, never
saturating. Only perturbation theory can tell: PT2 measures what the
external space supplies on top of each window. If the PT2-corrected
ladder saturates where the CAS ladder does not, the "missing" energy
was dynamic correlation in fact; if the total keeps moving, the
missing part is static and belongs in the window.

## Architecture

**The seam**: a certified CAS state in its orbital context → PT2
variant rows. The adapter module keeps the address `perturb.py` (the
analyzation keeps its home; fresh contract, PySCF-backed — the module
is NOT pure-numpy anymore; its docstring owns the new boundary).
Variants are a registry at the seam: `nevpt2-sc` now; CASPT2 joins
when an env ships it (a NEW conda env alongside the pinned 2.13.1 —
pre-approved); each variant produces the same row shape.

**Interface sketch** (one deep function family): consume the substrate
signature `(mol, mf, mo_full, ncore, ncas, nelecas)` + explicit
`roots`; internally: rebuild the multi-root CASCI on the stored MO
(the readout-twin recipe — `casscf.run_irrep_casci`'s pattern), assert
the rebuilt roots match the stored `ci_stack` BIT-EXACTLY (the
ticket-43 precedent), then `mrpt.NEVPT(mc, root=k).kernel()` per root,
extracting the 8 Koopmans-class `(norm, energy)` pairs (public module
functions `Sr/Si/Sijrs/Sijr/Srsi/Srs/Sij/Sir`). SA solutions enter
through the readout pass (the kernel refuses `StateAverageFCISolver`).
`density_fit` stays effectively inert (plain CASCI → exact 4-index
ERIS).

**The meters** (replacing EN's floor/w_ref, which died with EN):
per-class norms/energies; the effective per-class denominator
(energy/norm — the compression meter, SC's softened intruder read);
and the SYMMETRY PROBES (README · Analysis):

- **gauge redraw** — random unitary inside an exactly degenerate
  block; E_CAS/⟨L²⟩ as controls; E(2) spread = the disclosure. Per-L
  expectations: L=0 no gauge; L=1 spread ≈ 0 (transitivity — any
  spread is a covariance DEFECT, machinery certificate); L≥2 genuine
  spread possible (formulation property — canonical gauge primary).
- **component identity** — E(2) across pure M_L components equal by
  covariance; deviation = a window-closure defect, never gauge.
- **the canonical gauge**: pure M_L projections via our own Lz/L²
  analysis (the CI-order bridge) — the escalation if gauge spread
  reaches interpretation margins.
- PySCF's internal `canonicalize(cas_natorb=True)` is bypassable via
  `NEVPT(..., canonicalized=True)` — the designated hook for
  symmetry-controlled gauges, later if the probes demand it.

**The reference-class axis**: the same adapter runs NEVPT2 on
single-root CASSCF and CASCI cells (any solution class feeds the
seam); E(2) reads against each reference's symmetry-breaking
amplitude (⟨L²⟩ deviation, commutator norms — in the payloads). The
question: does PT2 compensate, penalize, or wash out a broken
reference.

## Delivery chain

`elements/<Atom>/perturbation.json` (payload: rows per cell × root ×
variant, per-class pairs, meters, full provenance) →
`elements/<Atom>/PERTURBATION_TABLE.md` + CSV mirror (ONE table
between markers; the generator refuses stale provenance — the
stale-render precedent) → the gate `scripts/validate_pt2.py` (the
certificate chain, green output committed to `logs/`) → README's
Analysis chapter is the reading discipline (landed 1adb82a). Commit +
push per batch; handoff at ticket boundaries.

## Steps (the gradual ladder)

0. **Smoke** — DONE 2026-11-19: examples 01 + 41 verbatim in-env; O2
   anchor agreement 3.2e-7 Eh (vs the example's retained reference);
   C2 SA pattern 0.3 s. The engine works in-env.
1. **Deletion** — DONE 2026-11-19: `perturb.py` + `test_perturb.py` +
   `validate_perturb.py` + the committed log removed; counts 352/258;
   spec/tickets/map/AGENTS lockstep (this commit).
2. **The adapter** — the new `perturb.py` + unit tests: O2/H2 anchors,
   per-class sum identity (Σe_class == e_corr), stored-roots
   bit-exactness, the SA guard exercised, class extraction pinned.
3. **Simplest real case** — one C SA cell, per-root rows; the two
   probes (the symmetry answer gets measured here).
3b. **Reference-class comparison** — NEVPT2 on stored CASSCF/CASCI
   cells (C, Si) vs SA, against symmetry amplitude.
4. **The ladder read + render** — C then Si across the window ladder;
   PERTURBATION_TABLE.md + reading key (R, gauge spread, reference
   class, basis thinness).
5. **Gated later** — Ti/Sc + the unc twins (U14e beyond); CASPT2 in a
   new env; wavefunction-level analyzation (⟨L²⟩/⟨S²⟩ of
   Ψ⁽⁰⁾+λΨ⁽¹⁾ — explicit first-order reconstructable from git
   history if the probes demand it).

## Verification

- (a) **Anchors**: O2/6-31G CASCI(6,8) → −0.169785157128 (the
  example's retained value; agreement band ~1e-6 Eh), H14/STO-3G →
  −0.094164359938.
- (b) **Identities**: per-class sum == e_corr (the kernel's own
  arithmetic, read twice); bit-exact stored-roots assertion before
  any PT2; NEVPT's SA guard exercised (raises on an SA solver).
- (c) **The probes**: gauge redraw controls flat (E, ⟨L²⟩), the
  component identity — the L=1 invariance is itself a machinery
  certificate.
- Renderer pins: byte-determinism, declared-roots policy. The gate is
  a run instrument (the `run.py` precedent), not a pytest tier; the
  per-commit gate stays fast-tier + targeted files.
- Test tier: unit tests on the frozen carbon fixture (NO SCF) for the
  pure parts (row shaping, class extraction, sum identity on
  synthetic MC objects); the campaign-size evidence priced at doors.

## Out of Scope (wave 1)

- CASPT2 construction (joins on the env upgrade; the reading caution —
  its level shift masks the intruder signal — lives in the reading
  key then).
- Multi-state / quasi-degenerate PT2 (per-state corrections only).
- Atoms beyond C and Si (Step 5's gated list).
- Campaign integration: no CellKey axis, no regen grammar change —
  post-processing on the stores.
