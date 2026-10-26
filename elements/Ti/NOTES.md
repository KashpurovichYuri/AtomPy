# Ti — NOTES (factual, dated; details behind TABLES.md)

- **Gated reference (ticket 12/21a).** Default ROHF = aufbau runaway:
  doubly-occupied 3d in the core, +123.07 mEh above 3d²4s² ³F; every
  energy-ordered start lands there. Fix: character-forced occupations
  (S-projector weights). Second layer: a metastable +11.2 mEh basin
  with the SAME occupations (polarized core smears later constructions
  at 1.6e-2 vs 3e-8 window defect); fix: the energy gate with cycled
  guesses, |E − E_ref| ≤ 1e-4. E_ref self-calibrated by diagnostics.
- **Escape basin.** Free single-root: −31.3 mEh, L ≈ 0 window
  ({s,p,d}-mixed), not a term; deepest minimum, identity-destroying.
- **Swap basin.** Free SA: the exact 4s↔3p trade (window p = 1/6),
  totals −10.1 mEh, ⟨L²⟩ 12.003, commutator O(1), term table scrambles
  (³F 38 vs 35, ³G 25 vs 18) — a different partition; totals only.
- **JT class.** d² member tilts descend (mechanism: README · Motivation
  deformations block); sa_radial_pur stops above its floor
  (1.4e-4/5.3e-4) — flat gauge, no restoring force.
- **Raw-core caveat.** Ground ⟨L²⟩ UNSHIFTED at percent-class
  commutators (4.9e-2/2.6e-1): ⟨L²⟩ alone is not a raw-core detector.
- **Pur substrate.** Machine-exact closures (6.2e-14) on every pur
  rung; sa_pur = partial swap (s-slot ⅔ s + ⅓ p, −3.8 mEh).
- **⟨L²⟩ = 0 MO windows (open).** Same class as the U ⁵S° reference;
  AOC-smearing hypothesis; see elements/U/NOTES.md.
- **MOs (pick) — 2026-10-26 (b).** The MOs start = the stock default:
  canonical ROHF MOs, active slice `[ncore : ncore+ncas]` by energy
  (`run_casscf_pyscf_start` passes no start; PySCF's own pick). Here
  that slice is 3d(5) + 4s(1) — it coincides with the s/d-class
  content for this window. Distinct from the SCF *occupations*, where
  energy order runs away and the reference forces by character
  (fact-class 15c/21a). Other picks (character-class, energy-window)
  are admissible future slots; the axis is stated per slot, not
  hardcoded.
- **⟨L̂z²⟩ column — 2026-10-26 (b).** Ground-root ⟨L̂z²⟩ in the rung's
  own fold (`scripts/lz2_analysis.py`; v|L²|v matches the stored
  ground_l2 to 1e-14 on all 17 rungs). Raw-fold detector: 4.0000
  exactly = an |M_L| = 2 root (the Λ branching). Exact-fold values are
  degenerate-cluster gauges (starred). OPEN CHECK: the ladder
  validation (joint H ∧ L̂z resolution, expected 0,1,1,4,4,9,9 + 0,1,1
  on CAS(2,5)) rides with the real→complex transform task, and settles
  the 4.000-vs-1/3/3/2 root-group tension.
