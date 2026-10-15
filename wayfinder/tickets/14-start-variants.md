# Ticket 14: Constrained-optimization & start-variant axis (window/provenance/constraint study)

Status: RESOLVED (pre-item refactor at a024d26; experiments + facts 16/17
in the follow-up commit; 70 tests green over two consecutive cold suite
runs)

## Question

The active window's PROVENANCE is a free design choice, and so is what the
ORBITAL OPTIMIZER may do to it. Sc (ticket 11) showed the free choices
break rotational closure profitably: free CASSCF escapes to an L≈0 basin,
shell-start SA restores <L^2> but swaps the window's 4s for a core 3p
(m-incomplete window, [H, Lz] = 0.83). Can we make the symmetry-vs-energy
trade-off MEASURABLE: impose closure on the optimizer and compare against
the broken-symmetry variational minimum? And can we prove the divergence
guard works by deliberately breaking orthonormality?

## Resolution (to implement)

1. **Constraint rungs** (the closure-preserving optimizer): rotations
   allowed only between functions of the SAME (l, m) -- for AO-ordered
   windows this is pure RADIAL mixing (3d_m may borrow 4d_m/5d_m, the 4s
   any s; m-structure and l-content untouchable). NOTE: "same-l only" is
   NOT sufficient -- a general rotation inside the d block mixes m's and
   breaks closure exactly like l-mixing. New rungs on Sc (+ C/Si controls):
   - `casscf_radial`: single-root, same-(l,m) constraint.
   - `sa_radial`: SA(5) over the 2D M_L set, same-(l,m) constraint --
     "imposed-symmetry SA"; compare vs free SA (which restores <L^2> but
     not [H, Lz]) and check whether the constraint kills the sa/sa_pyscf
     basin split.
2. **Enlarged complete window**: CAS(3,11) over {3d, 4s, 4p} (605 dets) --
   free optimization in a window that is m-complete for every l it
   contains. Question: was the Sc escape the optimizer crudely buying 4p
   correlation? If the big window reaches lower E at exact closure, yes.
3. **Start-variant axis**: `start=shell|canonical|shell+virtual|avas` --
   PySCF AVAS with our shell labels as independent cross-validation of
   the window span (ecosystem precedent: examples/mcscf/43-avas.py).
4. **Negative control**: deliberately feed the pre-fix-style
   non-orthonormal matrix; assert the guard fires (ortho_err > threshold,
   DIVERGED flag; C survives, Sc/Si-class divergence).
5. **Comparison view**: per-rung symmetry-vs-energy table (dE vs a frozen
   reference, ||[H,Lz]||, ||[H,L2]||, ground <L^2>, window l-weights, core
   anisotropy) as a first-class artifact, so new experiments are readable.

Pre-items (architecture, from the ticket-11 post-mortem): rung registry
(Rung dataclass; run_ladder = build refs -> run registry -> report),
declarative constraint spec (one generic constrained-CASSCF module over
PySCF's variable freezing), deduplicate the fold-and-analyze logic shared
by pipeline.stage_analysis and ladder._rung.

## Acceptance

- Constraint-rung table for Sc (+ C/Si control); measured answer to "how
  much of the escape gain is radial vs l/m-breaking" and "does the big
  window beat the escape at exact closure"; AVAS span agreement asserted;
  guard-fire test green; new facts in AGENTS.md.

## Resolution

Machinery: `drivers.lm_sectors` / `radial_basis` / `radial_mask` (the
(l, member)-labeled substrate; construction lessons in fact 16e — fused
shells via bas_nctr, member-preserving Gram-Schmidt because the symmetric
Loewdin gauge is free within degenerate member sets, strong-vs-needle
directions split at 1e-4), two new registry rungs (casscf_radial,
sa_radial), `rung_tags=` subset selection, comparison columns in _report.

Answers to the five questions (facts 16/17):

1. "How much of the escape gain is radial?" — on p² atoms ALL of it: the
   radial rung equals the free minimum to <=3e-7 Eh at exact shell purity,
   so the m-breaking is an energetically flat member GAUGE. On Sc NONE of
   it below the ceiling: the escape basin is ~25 mEh BELOW the best radial
   solution (stable: +24.4/+25.1 mEh) — it needs l/window-identity change.
2. "Does the complete {3d,4s,4p} window beat the escape at exact
   closure?" — CAS(3,9) (not 11: 5+1+3=9 orbitals, 324 dets) gains ~25 mEh
   of static 4p correlation at exact closure at CASCI, but free CASSCF
   STILL escapes (window grabs even f): the escape is an optimizer
   attractor, not correlation hunger; complete start != complete
   converged.
3. "Does sa_radial restore BOTH commutators to ~1e-13?" — C yes (8e-11);
   Sc no: L² = 5.9997 exact but commutators land at ~1e-3, which turned
   out to be the INPUT window's core-orthogonalization baseline (fact 17)
   — SA suppresses the member asymmetry (4e-2 single-root -> ~2e-3) but
   cannot go below the baseline. Removing it = ticket 13.
4. "Does the constraint kill the sa/sa_pyscf basin split?" — with window
   identity frozen there is no basin to split (sa_radial reproducible to
   µEh across cold runs). Bonus finding: the FREE shell-start sa is itself
   BIMODAL across chaos (4s<->3p swap or healthy window) — fact 15g.
5. "Attribute the relax-rung comm (~1e-3)" — RESOLVED as fact 17: it is
   the input window's defect (measured pre-optimization: 6.76e-4 == the
   relax rung's converged value), not core polarization (optimized core +
   raw window refold commutes at 1e-13) and not the optimizer (sa_radial
   with frozen canonical core shows the same 1e-3).

Plus: AVAS cross-validation (ncas/nelecas exact, span overlaps >= 0.985
Sc / 0.917 C); guard negative control (non-orthonormal input -> DIVERGED).
Test discipline learned: single-root radial commutators/L and the TM SA
basins are chaotic — pin structural invariants (window identity, ortho,
energy gaps) and branch on window weights, never the gauge.
