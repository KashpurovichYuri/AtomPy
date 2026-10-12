# Active space choice & orbital handling

## Question

Which orbitals represent the CAS(2,3) active space, and how are they fed to
PySCF's CASCI/CASSCF?

## Resolution (closed)

- **Active orbitals = the raw 2p shell AO columns** (columns of the identity),
  not canonical MOs. Canonical valence p MOs mix 2p/3p with *different radial
  profiles per m-axis*, so their span is not rotation-closed; the shell AOs
  (one contracted shell = same radial function, self-orthonormal solid
  harmonics) are exactly closed under L.
- Full MO matrix = [canonical core MOs | shell AOs | orthonormal complement],
  active window at `[ncore : ncore+ncas]` as PySCF requires.
- Core folded into h1 with the **spherically averaged** (l=0) core density —
  removes ~1e-5 SCF noise that would split ³P by ~4e-5 Eh.
- CASSCF/SA-CASSCF stages are analyzed on their *optimized* active orbitals;
  the resulting non-exact ⟨L²⟩ and finite ‖[H, Lz]‖ are physical results
  (orbital-optimization symmetry breaking), reported as such.
