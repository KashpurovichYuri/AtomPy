# Ticket 14: Start-variant axis for the ladder (window provenance study)

Status: open

## Question

The active window's PROVENANCE is a free design choice (raw shell AOs,
canonical MOs, shell + Löwdin-virtual picks, AVAS projection). The ecosystem
already validates shell-projected starts (PySCF examples 43-avas.py: AVAS
builds the active space from `ao_labels = ['Fe 3d', 'C 2pz']` and feeds the
projected orbitals to CASSCF; also 43-dmet_cas, 34-init_guess_localization,
42-compare_cas_space). Can the ladder treat window provenance as a controlled
axis — and can we prove our divergence guard works by deliberately breaking
orthonormality?

## Resolution (to implement)

1. `run_ladder(..., start=...)`: "shell" (default, exact closure),
   "canonical" (already exists as casscf_pyscf/sa_pyscf rungs -- unify),
   "shell+virtual" (shell m-set + one complement/Löwdin virtual pick --
   NOT m-complete: approximate-analysis rung), "avas" (delegate to
   pyscf.mcscf.avas with our shell labels; on C/Si should reproduce our
   window span -- independent cross-validation).
2. Negative control: a rung deliberately fed the pre-fix-style
   non-orthonormal matrix (raw stack, no block safety). Assert the guard
   fires: ortho_err > threshold, DIVERGED flag; quantify the sensitivity
   gradient (C survives, Si diverges). Test-only, not a default rung.
3. Physics caveat kept explicit: only shell windows are exactly closed;
   alternate starts are measured-and-labeled, never gated.

## Acceptance

- Start-variant table for C/Si (and Sc if ticket 11 has landed);
  AVAS-vs-shell span agreement asserted; guard-fire test green.
