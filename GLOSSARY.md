# Glossary

Domain vocabulary for AtomTerm. Names here are load-bearing: code, tests,
tickets, and AGENTS.md use them exactly.

## The calculation

- **atom config** (`AtomConfig`): symbol, basis, spin (2S), charge — the
  single input object; validated to one atom at the origin, spherical basis.
- **reference** (`mf`): an SCF solution of the atom. The ROHF reference is
  symmetry-broken (localized open-shell member); the fractional-occupation
  reference is spherically symmetric by construction.
- **core**: the inactive orbitals `[0:ncore]`. Folded into the active h1
  (mcscf convention); the fold may use the raw core density or its
  **spherical average** (the l=0 component — the physical reference).
- **window** (active space): the `ncas` orbitals at `[ncore:ncore+ncas]`.
  Canonical choice: the **raw shell AOs** (one or more complete (n,l)
  shells) — exactly rotation-closed, which is what makes the L analysis
  exact. Closure claims: l-purity is necessary, NOT sufficient (see
  m-completeness).
- **sector**: a fixed (na, nb) determinant space; Ms = (na−nb)/2. The
  primary sector follows cfg.spin; extra Ms sectors may exist.

## The analysis

- **fold**: active-space integrals with the core field included
  (`active_space_integrals`) — spherical_core, c_core are its conventions.
- **stage**: one finished calculation (CASCI, CASSCF, SA, ...) returning a
  `CASResult`.
- **stage facts** (`StageFacts`, `stage_facts`): the ONE analysis pass over
  a stage — fold with the stage's own core, term analysis, energies vs
  PySCF, window purity, ortho guard. Presentation-free; pipeline and
  ladder only render it.
- **commutator norms** ‖[H, Lz]‖, ‖[H, L²]‖: THE closure measurement. Zero
  = the window is rotation-closed for this H.
- **window purity** (`shell_weights`): S-metric l-weights of the window
  columns. Necessary-not-sufficient for closure (m-completeness below).
- **m-completeness**: a window carries the full m-set with SHARED radial
  vectors per l — the actual closure condition. An l-pure window of radial
  mixtures can still be m-incomplete (the free-CASSCF breaking, fact 8).
- **ortho_err**: max |C†SC − I| of a stage's MO matrix; > 1e-6 = DIVERGED,
  energies meaningless (ticket-09 guard).

## The experiment axis

- **rung**: one orbital-treatment variant of the same window (a `Rung`:
  stage-to-run + how to fold on it). The ladder = the rung registry +
  all-vs-all comparison.
- **rotation mask**: an `uniq_var_indices`-shaped callable declaratively
  constraining which orbital pairs the CASSCF optimizer may rotate
  (core↔virtual only, no-core↔window, same-(l,m) radial ...). The seam of
  `constrained_casscf`.
- **radial constraint** (same-(l,m)): rotations only between functions of
  the same l and member index — pure radial mixing; keeps every slot's
  (l, m) character. NOTE: exact closure additionally requires the radial
  rotation to be applied IDENTICALLY to all m partners (shared radial
  vectors); independent per-member radial mixing preserves each slot's
  purity but can still break [H, Lz] at the level of the radial-vector
  differences. Measured by the radial rungs.
- **escape**: free CASSCF leaving the intended window basin (Sc: dominantly
  L=0 basin of the escaped window — optimizer artifact, not a physical
  term).
- **max|dE|**: our reconstructed energies vs PySCF's own — a fold/bookkeeping
  check with known physical offsets (polarized-core detector).
