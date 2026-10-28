"""AtomTerm: atomic term symbols (^{2S+1}L) on top of PySCF.

ROHF -> CASCI / CASSCF / SA-CASSCF for an atom, then classification
of the states by exact Lz / L^2 / Sz / S^2 operator matrices built
directly in the Slater-determinant basis (a second-quantized twin
exists only as an elementwise cross-check).

No root API -- import module-qualified:
``from atomterm.windows import build_full_mo``.
"""
