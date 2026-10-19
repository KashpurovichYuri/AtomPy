"""AtomTerm: atomic term-symbol (L-multiplicity) analysis on top of PySCF.

Pipeline: ROHF -> CASCI -> CASSCF -> SA-CASSCF for an atom, then classification
of the obtained states by exact Lz / L^2 / Sz / S^2 operator matrices built in the
Slater-determinant basis (first quantization, no second quantization anywhere).
"""

from atomterm.system import AtomConfig, build_mol
from atomterm.detbasis import (
    determinant_basis,
    slater_condon_matrix,
    commutator,
    frobenius_norm,
)
from atomterm.angular import l_matrices, l2_matrix, l_matrices_analytic

__all__ = [
    "AtomConfig",
    "build_mol",
    "determinant_basis",
    "slater_condon_matrix",
    "commutator",
    "frobenius_norm",
    "l_matrices",
    "l2_matrix",
    "l_matrices_analytic",
]
