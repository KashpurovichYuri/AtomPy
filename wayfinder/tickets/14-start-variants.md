# Ticket 14: Constrained-optimization & start-variant axis (window/provenance/constraint study)

Status: open

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
