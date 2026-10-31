"""The SO(3)->D2h dictionary pins: branching, the XOR product
structure, and the determinant map -- each anchored against PySCF's
own symmetry machinery."""

import numpy as np
import pytest
from pyscf import ao2mo, fci, gto
from pyscf.symm import param as symm_param

from atomterm.irreps import (
    IRREP_IDS,
    IRREP_NAMES,
    _CHARACTERS,
    _SIGNS,
    branch_d2h,
    det_irreps,
    orb_irreps,
)
from atomterm.windows import shell_slots

pytestmark = pytest.mark.fast

# _CHARACTERS/_SIGNS micro-pins: the flip encoding IS the module's
# derivation core -- the public functions cannot reach it directly.


def _d2h_mol(basis) -> "gto.Mole":
    """A D2h atom over an uncontracted one-primitive-per-l basis (m
    labels unambiguous, whole-l spans exact)."""
    return gto.M(atom="C 0 0 0", basis=basis, symmetry="D2h", verbose=0)


_SHELLS_BASIS = gto.parse("""
C S
  1.0  1.0
C P
  1.0  1.0
C D
  1.0  1.0
C F
  1.0  1.0
""")


def test_branch_multiplicities_sum_to_2l1():
    counts = {(L, p): sum(branch_d2h(L, p).values())
              for L in range(9) for p in ("g", "u")}
    assert all(v == 2 * L + 1 for (L, _), v in counts.items())


def test_branch_hand_known_patterns():
    # the one-electron shell patterns (textbook tessereal splitting)
    assert branch_d2h(0, "g") == {"Ag": 1}
    assert branch_d2h(0, "u") == {"Au": 1}
    assert branch_d2h(1, "u") == {"B1u": 1, "B2u": 1, "B3u": 1}
    assert branch_d2h(2, "g") == {"Ag": 2, "B1g": 1, "B2g": 1, "B3g": 1}
    # the f shell: Au holds xyz, each B_u holds two components
    assert branch_d2h(3, "u") == {"Au": 1, "B1u": 2, "B2u": 2, "B3u": 2}


def test_branch_matches_real_shell_labels():
    mol = _d2h_mol(_SHELLS_BASIS)
    ids = orb_irreps(mol, np.eye(mol.nao_nr()))
    for _, l, _, rows in shell_slots(mol):
        parity = "g" if l % 2 == 0 else "u"
        expected = {IRREP_IDS[name]: n
                    for name, n in branch_d2h(l, parity).items()}
        actual = dict(zip(*np.unique(ids[rows], return_counts=True)))
        assert {int(k): int(v) for k, v in actual.items()} == expected


def test_character_table_equals_pyscfs():
    # PySCF's table rows are in id order, columns in the same op
    # order -- the flip derivation must reproduce it elementwise.
    pyscf_table = np.array([row[1:] for row in
                            symm_param.CHARACTER_TABLE["D2h"]])
    assert np.array_equal(_CHARACTERS, pyscf_table)


def test_irrep_product_is_xor():
    for a, sa in _SIGNS.items():
        for b, sb in _SIGNS.items():
            product = tuple(s * t for s, t in zip(sa, sb))
            assert _SIGNS[a ^ b] == product


def test_det_irreps_doubly_occupied_is_identity():
    # (1,1) sector on two orbitals: same-orbital pairs cancel to Ag,
    # mixed pairs carry the product B1u ^ B2u = B3g
    assert det_irreps([5, 6], 1, 1).tolist() == [0, 3, 3, 0]


def test_det_irreps_support_equals_wfnsym_solver():
    mol = _d2h_mol("sto-3g")
    c = np.eye(mol.nao_nr())[:, 2:5]  # the 2p shell: tesseral AOs
    orbsyms = orb_irreps(mol, c)
    ids = det_irreps(orbsyms, 2, 0)
    h1 = c.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ c
    eri = ao2mo.kernel(mol, c)
    for ir in np.unique(ids):
        solver = fci.direct_spin1_symm.FCI(mol)
        solver.orbsym = orbsyms
        solver.wfnsym = IRREP_NAMES[ir]
        solver.conv_tol = 1e-12
        _, ci = solver.kernel(h1, eri, 3, (2, 0))
        theirs = set(np.flatnonzero(
            np.abs(np.asarray(ci).ravel()) > 1e-10).tolist())
        ours = set(np.flatnonzero(ids == ir).tolist())
        assert theirs == ours, ir


def test_invalid_inputs_fail_fast():
    with pytest.raises(ValueError):
        branch_d2h(-1, "g")
    with pytest.raises(ValueError):
        branch_d2h(1, "x")
    with pytest.raises(ValueError):
        det_irreps([8], 1, 0)
