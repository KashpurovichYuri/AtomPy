# Si — NOTES (factual, dated; details behind TABLES.md)

- **Ticket 09 (the multi-solution that wasn't).** Si's CASSCF/SA once
  looked like it escaped to an ³S-like basin; the cause was a
  NON-orthonormal `build_full_mo` (raw 3p AOs overlap the canonical
  2p-core MOs by ~0.3; the complement projector became inexact; the
  optimizer diverged below the variational bound). With the block-safe
  construction Si behaves exactly like carbon. The divergence guard
  (`ortho_err` / DIVERGED) exists for this class.
- **Gauge draws.** Radial rung == free rung to ≤3e-7 Eh with exact
  shell purity — the entire gain is radial; which member gauge a
  single-root optimizer lands in is chaotic (Si drew 1.2e-1 member
  asymmetry where carbon drew 5e-11).
- **Fact-17 baseline.** Every Si optimized rung inherits ~1e-3 from its
  STARTING window (core-orthogonalization against a non-scalar core
  span); the pur substrate removes it (rungs close machine-exact).
