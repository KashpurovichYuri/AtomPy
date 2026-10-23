# C — NOTES (factual, dated; details behind TABLES.md)

- **Confinement finding.** CASCI over the ROHF's own 2p-like orbitals
  (canonical window) reproduces the ROHF total exactly (ΔE = 0.00):
  the ³P pair carries no net correlation the CAS can express; the
  +250.8 mEh "CASSCF gain" is the recovery of the raw-shell
  confinement cost. The optimized window (~2/3 2p + 1/3 3p at
  θ\* = 0.61) rebuilds the relaxed ROHF-like orbital.
- **Flat-gauge census.** Five fresh draws of the single-root radial
  rung: one common energy (spread < 1e-6 mEh) while ⟨L²⟩ ranges
  1.69–2.00; one draw landed the symmetric window itself. Endpoint =
  noise on a flat direction.
- **Free-CASSCF anatomy.** The tilted window keeps ‖[H,Lz]‖ ≈ 5e-7
  (single-irrep sector — the commutator is vacuous there) while
  ‖[P,L⊥]‖ ≈ 1.1 exposes it; the window-projected L̂² eigenvalues
  read 1.723/1.862/1.862 against exactly {2, 2, 2} on the raw shell.
  Each tilted member is still an l = 1 m-eigenfunction; the ladder
  cross-terms carry the defect.
- **Pur substrate.** The pur core closes both input windows at machine
  precision; the pur family carries a window-definition offset vs the
  raw-fold CASCI rows (bound, not pinned).
