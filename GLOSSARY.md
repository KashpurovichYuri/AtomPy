# Glossary

Plain-language vocabulary for AtomTerm. Every term the reports and code
use is defined here; the precise statement follows the plain one.

## The objects

- **atom config** — the input card: which atom, which basis, which spin
  (2S), which charge. Everything else is derived.
- **reference** — a solved mean-field (an SCF). The ROHF reference is
  *lopsided*: it occupied one specific open-shell orbital (say px), so
  its density remembers that choice. The fractional-occupation reference
  puts equal fractions of an electron in every open-shell orbital and is
  round (spherical) *by construction*.
- **core** — the orbitals we do NOT correlate (everything below the
  window). The core enters the active problem as an electric field: the
  "fold".
- **fold** — building the effective Hamiltonian for the active electrons
  with the core's field included (one-body part) instead of treating
  core electrons explicitly.
- **window** — the active space: the few orbitals whose electrons we
  correlate. The natural choice for an atom: the raw functions of one
  valence shell (e.g. the three 2p). A window made of *complete shells*
  is rotation-closed, which is what makes the angular-momentum analysis
  exact.
- **sector** — the determinant space for fixed numbers of alpha and beta
  electrons; the label is Ms = (na − nb)/2. Atoms need several sectors
  to hold all terms of a configuration.

## The analysis

- **member** — the m-component of a shell: px, py, pz are the three
  members of p. Under rotation they mix into each other.
- **m-completeness** — the real closure condition: a window is safe under
  rotations when, for each l it contains, it holds ALL members and those
  members share the SAME radial mixture (all built from the same
  combination of 2p/3p/…). l-purity alone is not enough.
- **commutator norms** ‖[H, Lz]‖, ‖[H, L²]‖ — the symmetry meters. Exactly
  zero when core+window is rotation-closed; their size is a measured
  number for every calculation ("rung") we do.
- **window purity** — how much of each window orbital lies in its
  nominal l (1.000 = no leakage). Necessary, not sufficient: a pure-p
  window can still be m-incomplete (see above).
- **stage / stage facts** — a finished calculation (CASCI, CASSCF, …) is
  a *stage*; `stage_facts` performs the single analysis pass every stage
  gets: fold, classify terms, compare energies with PySCF, check purity
  and orthonormality.
- **ortho_err / DIVERGED** — max |C†SC − I| of a rung's orbitals. Above
  1e-6 the optimizer lost orthonormality and every number of that rung
  is meaningless; the report says DIVERGED.
- **max|dE|** — our energies vs PySCF's own for the same calculation.
  Zero means identical bookkeeping; known physical offsets (e.g. we
  average a polarized core, PySCF uses it raw) show up here on purpose.

## The experiment axis

- **ladder** — the comparison series: the same atom and window analyzed
  under every orbital treatment, from frozen to fully optimized.
- **rung** — one step of the ladder (see RUNGS.md for what each rung
  holds fixed and releases).
- **rotation mask** — a rule declaring which pairs of orbitals the
  optimizer may rotate. Masks are how we impose physics questions on
  PySCF's optimizer ("only radial mixing", "no core exchange", …).
- **radial constraint** — the mask "same (l, member) only": 2px may mix
  with 3px, never with py or any s/d. Pure radial breathing.
- **gauge** — a choice of description that does not change the physics
  (like a phase). For us: *which* linear combination of the members each
  window orbital is — the span is the same, the energy is the same, only
  the description differs.
- **flat direction / flat member gauge** — a deformation of the window
  along which the energy does not change (measured < 1e-6 Eh on C/Si).
  Giving each member a *different* radial mixture is such a direction:
  it breaks the commutators (first order) but not the energy (second
  order at a stationary point). Optimizers wander along it for free —
  which is why free CASSCF "breaks symmetry" without gaining energy.
- **escape / escape basin** — on Sc the *free* optimizer rebuilds the
  window from ~half p-character (borrowed from virtual 4p and core 3p):
  the active space is no longer the 3d/4s shell at all. The "basin" is
  that solution's region of attraction; it is a real local minimum
  (variational, orthonormal), just not the window we intended. Not
  reachable by radial mixing (the radial ceiling sits ~25 mEh above it).
- **input-window baseline** — the ~1e-3 commutator level that every Si/Sc
  optimized rung inherits from its STARTING window: to make it
  orthonormal we must orthogonalize the shell AOs against the core, and
  the canonical core carries small asymmetric tails, so the projected
  window leaks under rotations before any optimization runs (fact 17).
  Not core polarization; not the optimizer. Fix = the closed core span
  (ticket 13).
- **closure repair** — (planned) post-processing that re-shares the
  radial mixture equally among members, restoring exact closure at
  (predicted) ~µEh energy cost, because that direction is flat.
