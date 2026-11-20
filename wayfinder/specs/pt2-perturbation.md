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
the rebuilt roots match the stored `ci_stack` (the ticket-43
precedent, in its DECLARED form: energies strict at 1e-11 absolute;
CI amplitudes per-root at 1e-8 with the exactly degenerate block's
lawful gauge freedom honored through the joint projector — the store
keeps flat vectors, the reader (n,1) columns, amplitudes are the
assertion), then `mrpt.NEVPT(mc, root=k).kernel()` per root,
extracting the 8 Koopmans-class `(norm, energy)` pairs (public module
functions `Sr/Si/Sijrs/Sijr/Srsi/Srs/Sij/Sir`, plus the kernel's
canonicalization prelude mirrored through the module's internal ERIS
/ contraction helpers — the plumbing-vs-kernel twin pins it), SA
solutions enter through the readout pass (the kernel refuses
`StateAverageFCISolver`).
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
2. **The adapter** — DONE 2026-11-19: `perturb.py` rebuilt at the
   seam; 9 unit tests on the frozen carbon fixture (NO SCF): the
   exact ³P block on the pur core, the raw-core split pinned 4.4e-5
   as the unit-scale detector, the gauge-tolerant stored-roots
   refusal, the plumbing-vs-kernel twin at 1e-10 on the small CAS;
   4 anchors (O2/H14 at 1e-6, two-root H2 root selection, the
   kernelled SA guard); counts 365/271.  DISCOVERY: PySCF's
   SC-NEVPT2 wobbles ~3e-7 Eh run-to-run (thread-parallel reductions
   in the big contractions; the kernel's own three-run record
   −0.16978547/−0.16978542/−0.16978512 on O2) — the twin pins sit at
   that floor; the class plumbing shares it.
3. **Simplest real case** — DONE 2026-11-19 (`scripts/probe_pt2_block.py`,
   the certified pur cell `c_2p_2-0_sa-casscf_ao_-_pur_avg`): E(2) =
   −58.1340 mEh (**−1.5819 eV**) IDENTICAL across the ³P block;
   probe A (8 redraws × 3 roots): controls 1.2e-14, spread
   **7.3e-13 eV** — the engine wobble floor; probe B: pure M_L
   states exact (−1/0/+1), the real ML=0 component reads the same
   constant.  SYMMETRY KEPT, per-root rows stand as stored.  The
   guard's first real bite: store-flat vs reader-column CI shapes —
   the assertion compares amplitudes (fixed, pinned).  Note: E(2) is
   ~7.5× the EN-era preview (−7.7 mEh) — the H₀-choice axis (Dyall
   denominators vs EN's exact-configuration diagonals) is a big
   effect; the ladder read quantifies it per window.
3b. **Reference-class comparison** — DONE 2026-11-19: NEVPT2 on
   stored CASSCF/CASCI cells (C, Si) vs SA, against symmetry
   amplitude; E(2) rotation-insensitive at every measured amplitude.
4. **The ladder read + render** — DONE 2026-11-19/20 (C and Si, then
   the whole campaign): PERTURBATION_TABLE.md + reading key (R, the
   gauge spread, the reference class, the basis thinness); the
   optimized ladder SATURATES on C/Si (the optimizer repairs the
   poverty class; the thin spans fuel, static remainder 43.2/~41
   mEh).
5. **Gated later** — DONE 2026-11-20 for the ladder: Ti/Sc + the
   unc twins + U14e collected/rendered/gated (the combined gate ALL
   GREEN, 1259 roots); the CAS-size flagship read: the corrected
   ladder PREFERS (4,12) -- the (4,15) promotion double-counts
   dynamic correlation.  Still open here: CASPT2 in a new env (a
   disclosed variant on the env upgrade); wavefunction-level
   analyzation (declined -- the resurrection path recorded).  The
   lane seam (the follow-up) DONE 2026-11-20: cells.readout_lane +
   nevpt2_rows_on; the lane rows priced on C/Si/Sc/Ti/U14e (the
   declared-out legs recorded in ticket 45).

## Verification

- (a) **Anchors**: O2/6-31G CASCI(6,8) → −0.169785157128 (the
  example's retained value; agreement band ~1e-6 Eh), H14/STO-3G →
  −0.094164359938.
- (b) **Identities**: per-class sum == e_corr (the kernel's own
  arithmetic, read twice); the declared stored-roots assertion before
  any PT2 (strict energies + gauge-aware amplitudes, per
  Architecture); NEVPT's SA guard exercised (raises on an SA solver).
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
- Atoms beyond C and Si (Step 5's gated list) — SUPERSEDED
  2026-11-20: the whole campaign joined (the user's directive); the
  lane legs' fresh constructions needed the regen grammar's `lane`
  group (the store-recovery loop, not a grammar change) -- the
  post-processing principle stands for the GRID stores.
- Campaign integration: no CellKey axis, no regen grammar change —
  post-processing on the stores.
