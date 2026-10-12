# Grilling: scope & configuration

## Question

Which system, basis, active space, state-averaging, and validation strategy
for the first AtomTerm slice?

## Resolution (closed)

- **System**: carbon atom at the origin, triplet (spin=2), cc-pVDZ (spherical).
- **Active space**: CAS(2,3) — the three 2p orbitals only; textbook ¹S/¹D/³P.
- **SA-CASSCF**: 3 equal-weight triplet roots (Ms=1) — the ³P M_L components.
- **Operator validation**: use the libcint integral machinery inside PySCF
  (user's call) rather than a hand-rolled grid quadrature; an analytic
  overlap-based construction serves as the independent cross-check.
- Deliverables: .chkfile caching at every stage; L-multiplicity analysis;
  [H, Lz] and [H, L²] commutator tests; AGENTS.md.
