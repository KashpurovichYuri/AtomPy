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
- **JT class.** d² member tilts descend (mechanism: README §2.4
  deformations block); sa_radial_pur stops above its floor
  (1.4e-4/5.3e-4) — flat gauge, no restoring force.
- **Raw-core caveat.** Ground ⟨L²⟩ UNSHIFTED at percent-class
  commutators (4.9e-2/2.6e-1): ⟨L²⟩ alone is not a raw-core detector.
- **Pur substrate.** Machine-exact closures (6.2e-14) on every pur
  rung; sa_pur = partial swap (s-slot ⅔ s + ⅓ p, −3.8 mEh).
- **⟨L²⟩ = 0 MO windows (open).** Same class as the U ⁵S° reference;
  AOC-smearing hypothesis; see elements/U/NOTES.md.
