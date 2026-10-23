# Sc — NOTES (factual, dated; details behind TABLES.md)

- **Escape basin.** Free single-root CASSCF: −31.5 mEh below ROHF, an
  ²S-like window ({s,p,d}-mixed, ⟨L²⟩ 0.05) — variational and
  orthonormal but not a term; the deepest minimum. Not radially
  reachable (the radial ceiling sits ~+25 mEh above it — a window
  statement: the lattice crosses below via 4f, −3.06 mEh).
- **Swap basin.** Free SA: the 4s↔3p trade (window p = 1/6 exactly);
  totals −13.3 mEh; ⟨L²⟩ 6.002 readable; commutator O(1). A different
  partition — compare in totals only.
- **Bimodal SA (fact 15g).** The shell-start SA lands in the swap or
  the healthy basin across SCF chaos; the swap-basin commutator
  fresh-drew 0.046 (2026-10-23) against the pinned O(0.1–1) band —
  band widened; the p-weight stays the discriminator.
- **Frac-occ reference.** +1.7 Eh vs ROHF; its ²D terms are
  quasi-degenerate at the bottom and their ordering flips across SCF
  chaos — never pinned in tests.
- **Radial rung chaos.** Single-root radial rung: ²D→²P ground-term
  relabeling across draws; sa_radial reproducible to µEh with exact
  ⟨L²⟩ (the ensemble anchors the gauge).
- **Pur substrate.** relax_pur / sa_relax_pur at 1e-13 with unbroken
  Lz/L²; the swap attractor survives purification (branch on window
  weights); the flat-gauge floor stays (~1e-6..1e-4 window tilt).
- **NO windows (ticket 22).** Groups 4s(1)/3d(5)/4p(3)/[3d+4d](5); NO
  beats the shell lattice at every equal-det count measured and the
  ROHF-MO column at the top pairs; every group window at ~1e-13.
- **AVAS cross-check (fact 19).** The ecosystem window: right size,
  mean-field-inherited span — 0.9855 overlap with the shell span yet
  O(0.5) closure defects: overlap says nothing about closure.
