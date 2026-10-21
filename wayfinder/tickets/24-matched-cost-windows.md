# Ticket 24: matched-cost window comparison — AO-exact vs raw-ROHF-broken (U + Sc + Ti)

Status: OPEN (execute first among the new wave; light items only)

Origin: user direction 2026-10-21 (post-4f2d discussion). Frame: pure
AO windows make "below ROHF" expectable — enough span trivially does
it; the INFORMATIVE comparison is exact-symmetry vs broken-symmetry AT
MATCHED COST (same determinant count), plus where the two converge as
the window grows.

## The physics questions

1. **Does exactness cost energy at matched cost?** Window A (AO
   primitives, m-complete, ⟨L²⟩ exact) vs window B (raw ROHF columns,
   m-incomplete, ⟨L²⟩ broken) at the SAME det count: energies and
   symmetry fingerprints side by side. Expectation: A ≥ B in energy
   (B's extra freedom = member/asymmetric radials) — measure by how
   much, per atom class.
2. **Where does the fact-16b plateau melt?** (user's full-space
   intuition: exact and broken converge as the window → full space —
   the +24.4/+25.1 mEh radial-vs-escape gap on Sc MUST close once the
   lattice includes the escape's directions m-completely; the escape
   was a window-identity change, and identity is available exactly in a
   wide-enough m-complete window). Measure E_lattice vs the escape
   energy / radial ceiling on Sc (Ti as the second TM point).

## The runs (ALL short; ≤ 30 min/sector, nothing bigger — user rule)

- **U methodology-B scan**: window = the 4 occupied open MOs + up to 4
  character-selected virtuals (S-projector metric — canonical order is
  junk on this basis), admitted one at a time: CAS(4,4) → CAS(4,12)
  class. Per step: E, ⟨L²⟩ off-grid magnitude, ‖[H,Lz]‖, ‖[H,L²]‖ —
  the broken-symmetry-size-vs-extension curve. Minutes per step.
- **Sc/Ti lattices, methodology A**, both basis kinds side by side:
  contracted (shell granularity: {3d}, {4d}, {4s}, … — one converged
  radial per shell; select_valence_ao route) and decontracted
  (primitive granularity; primitive_chain route). Per-block ΔE tables
  (references EXPLICIT: Δf = E(k f) − E(k−1 f) at fixed d; Δd
  analogously) + ⟨L²⟩ columns; read against the escape/ceiling numbers
  of facts 15d/16b.
- **Sc/Ti methodology B at matched det counts**: the question-1
  verdict.

## Deliverables

- Comparison tables under results/ (per atom × methodology × basis
  kind), per-block ΔE matrices with defined references.
- **README §7.x section** (draft reviewed by the user IN-CONVERSATION
  before landing — the standing docs rule): the U window-lattice table
  + the brief per-block ΔE matrix + the representability comment (a
  CASCI bounds the ROHF only when the window span contains the
  occupied MOs; window nesting ≠ determinant nesting; the crossing
  measures how much radial span the mean field uses). Brief.
- §3 entries in ticket 20 for any new friction.

## Out of scope / postponed (user rule: no computationally intensive)

CAS(6,13) on U (even the short (5,1) sector — 16,731 dets — is
postponed, not cancelled); any sector > 30 min; 4f+2d-class runs; the
3rd-d-primitive U point; NO compression (ticket 22); SA wave (23).
