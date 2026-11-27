"""The analysis meters: <<L^2>> and <<S^2>> of CI roots in one orbital
set.

atomterm's angular operators through the matrix-free application;
a CI vector enters in PySCF's order and leaves through the CI-order
bridge.  The one piece of atomterm analysis the routes keep -- the
construction is double-certified in ``angular.py``; the meters are
the terms' identity, not part of what is under test.
"""

import numpy as np

from atomterm.detbasis import determinant_basis, operator_apply
from atomterm.operators import (ci_in_detbasis_order, l2_operator_parts,
                                s2_operator_parts)


class Meters:
    """The meters on ONE active block: <<L^2>> reads the block's AO
    columns, <<S^2>> is orbital-free.

    Parameters
    ----------
    mol : the molecule (the L^2 AO matrices read its basis).
    active : (nao, ncas) the block's S-orthonormal columns -- the
        orbitals the CI vectors were solved in.
    sector : (na, nb) the active electrons.
    """

    def __init__(self, mol, active, sector):
        norb = active.shape[1]
        dets = determinant_basis(norb, *sector)
        self._apply = {
            "l2": operator_apply(dets, *l2_operator_parts(mol, active)),
            "s2": operator_apply(dets, *s2_operator_parts(norb)),
        }
        self._order = (norb, *sector)

    def of(self, ci) -> tuple[float, float]:
        """(<<L^2>>, <<S^2>>) of one CI root (PySCF order, any shape)."""
        w = ci_in_detbasis_order(np.asarray(ci).ravel(), *self._order)
        return (float(np.real(w @ self._apply["l2"](w))),
                float(np.real(w @ self._apply["s2"](w))))
