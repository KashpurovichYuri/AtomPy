"""The analysis meters: atomterm's angular operators through the
matrix-free application; the many-body meters read per root through
the CI-order bridge."""

import numpy as np  # noqa: E402

from atomterm.detbasis import determinant_basis, operator_apply  # noqa: E402
from atomterm.operators import (  # noqa: E402
    ci_in_detbasis_order, l2_operator_parts, s2_operator_parts)

# the analysis meters --------------------------------------------------

def meters(mol, c_act, na: int, nb: int):
    """The analysis meters as their ACTION (atomterm's angular
    operators through the matrix-free application; the many-body
    meters read per root through the CI-order bridge).  The one piece
    of atomterm analysis the legs keep -- the construction is
    double-certified in angular.py; the meters are the terms'
    identity, not part of the derivation under test."""
    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    apply_l2 = operator_apply(dets, *l2_operator_parts(mol, c_act))
    apply_s2 = operator_apply(dets, *s2_operator_parts(norb))
    apply = {"l2": apply_l2, "s2": apply_s2}

    def expect(kind: str, ci_pyscf) -> float:
        w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(),
                                 norb, na, nb)
        return float(np.real(w @ apply[kind](w)))

    return expect


def meters_all(blocks, expect):
    """The meters over EVERY block component (in branch order) --
    independent of any PT2 root cut."""
    l2s, s2s = [], []
    for _, mc in blocks:
        for k in range(len(np.atleast_1d(mc.e_cas))):
            ci = np.asarray(mc.ci[k]).ravel()
            l2s.append(expect("l2", ci))
            s2s.append(expect("s2", ci))
    return l2s, s2s
