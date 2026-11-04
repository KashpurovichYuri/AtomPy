"""Unit pins of atomterm.reference: the ReferenceSpec contract, the
character-forced occupation hooks (ROHF and the per-spin UHF twin),
and the UHF contamination meter (numpy-only given a Mole -- no SCF
runs here; the SCF-bearing reference tests live in the integration
tier)."""

from types import SimpleNamespace

import numpy as np
import pytest

from atomterm.reference import (ReferenceSpec, _forced_rohf_occ,
                                _forced_uhf_occ, uhf_contamination)
from atomterm.windows import lm_sectors
from helpers import uncontracted_carbon_mole

pytestmark = pytest.mark.fast


def test_reference_spec_validation():
    """The occ-pattern contract: the two pattern lists come as a pair,
    and the basin gate needs them."""
    ReferenceSpec()  # plain reference: all-None is the default run
    ReferenceSpec(occ_closed=("4s",), occ_open=("3d",))
    ReferenceSpec(occ_closed=("4s",), occ_open=("3d",),
                  occ_e_ref_uhf=-759.73789123)
    with pytest.raises(ValueError, match="together"):
        ReferenceSpec(occ_closed=("4s",))
    with pytest.raises(ValueError, match="needs the"):
        ReferenceSpec(occ_e_ref=-848.40648875)
    with pytest.raises(ValueError, match="needs the"):
        ReferenceSpec(occ_e_ref_uhf=-759.73789123)


def test_forced_occ_l_specifiers_select_whole_subspaces():
    """'l0'..'l9' force occupations onto whole angular-momentum
    subspaces (the honest vocabulary on uncontracted bases)."""
    mol = uncontracted_carbon_mole()  # nelec = 6: closed = 2, open = 2
    rows_l = {l: sorted(r for (ll, m), rws in lm_sectors(mol).items()
                        if ll == l for r in rws) for l in range(4)}
    hook = _forced_rohf_occ(mol, ["l0", "l1"], ["l2"])
    occ = hook(np.zeros(mol.nao_nr()), np.eye(mol.nao_nr()))
    closed = set(np.flatnonzero(occ == 2.0))
    assert closed <= set(rows_l[0] + rows_l[1])       # s/p-pure block
    with pytest.raises(ValueError, match="matched no AO labels"):
        _forced_rohf_occ(mol, ["4z"], ["l2"])


def test_forced_occ_open_block_carries_a_strict_premium():
    """The open block is ranked by 2 w_open + w_closed: open-shell
    character beats closed leftovers even where they tie."""
    mol = uncontracted_carbon_mole()
    rows_l = {l: sorted(r for (ll, m), rws in lm_sectors(mol).items()
                        if ll == l for r in rws) for l in range(4)}
    hook = _forced_rohf_occ(mol, ["l0", "l1"], ["l2"])
    occ = hook(np.zeros(mol.nao_nr()), np.eye(mol.nao_nr()))
    open_ = set(np.flatnonzero(occ == 1.0))
    # on identity MOs the unchosen s/p rows (w_closed=1, w_open=0) tie
    # the d rows (w_open=1, w_closed=0) under a 1:1 sum; the 2:1
    # premium breaks it toward d. Measured on U: a 1:1 sum let DIIS
    # close s/p garbage (E > +1e3); the closed term still catches
    # early-cycle open character hiding in core-like MOs (Sc/Ti).
    assert open_ <= set(rows_l[2])


# --------------------------------------------------------- the UHF class

def _synthetic_uhf(a_cols, b_cols, nao, c_a=None):
    """A UHF-shaped namespace over identity (S-orthonormal) AOs:
    alpha occupies the columns ``a_cols`` of ``c_a`` (default the
    identity), beta the columns ``b_cols`` of the identity."""
    eye = np.eye(nao)
    occ_a = np.zeros(nao)
    occ_a[list(a_cols)] = 1.0
    occ_b = np.zeros(nao)
    occ_b[list(b_cols)] = 1.0
    return SimpleNamespace(mo_coeff=(eye if c_a is None else c_a, eye),
                           mo_occ=(occ_a, occ_b),
                           get_ovlp=lambda: eye)


def test_forced_uhf_occ_blocks_per_spin():
    """The per-spin hook: the closed pairs on the closed span for BOTH
    spins, the unpaired electrons alpha-only (the maximal-Ms counts),
    beta's open block empty."""
    mol = uncontracted_carbon_mole()  # nelec = 6, spin = 2: 2 pairs
    rows_l = {l: set(r for (ll, m), rws in lm_sectors(mol).items()
                     if ll == l for r in rws) for l in range(4)}
    hook = _forced_uhf_occ(mol, ["l0", "l1"], ["l2"])
    eye = np.eye(mol.nao_nr())
    occ = hook((np.zeros(mol.nao_nr()),) * 2, (eye, eye))
    assert occ.shape == (2, mol.nao_nr())
    beta = set(np.flatnonzero(occ[1] == 1.0))
    assert len(beta) == 2 and beta <= rows_l[0] | rows_l[1]
    alpha = set(np.flatnonzero(occ[0] == 1.0))
    assert len(alpha) == 4
    assert len(alpha & (rows_l[0] | rows_l[1])) == 2
    assert len(alpha & rows_l[2]) == 2
    with pytest.raises(ValueError, match="need the MOs"):
        hook((np.zeros(mol.nao_nr()),) * 2, None)


def test_uhf_contamination_paired_reference_is_the_floor():
    """The paired (ROHF-like) determinant: sigma exactly 1, <S^2>
    exactly S(S+1) -- the overlap formula's fixed point."""
    mf = _synthetic_uhf((0, 1, 2), (0, 1), 4)
    s2, sigma = uhf_contamination(mf)
    assert s2 == 0.75                      # the 2S+1=2 doublet
    assert np.allclose(sigma, 1.0)


@pytest.mark.parametrize("a_set,b_set", [
    ((0, 1, 2), (0, 1)),   # paired (the floor, exact S(S+1))
    ((0, 1, 2), (0, 3)),   # a displaced beta open: classic contamination
    ((0, 2, 3), (1, 2)),   # permutation overlap: one shared, one disjoint
])
def test_uhf_contamination_matches_exact_s2(a_set, b_set):
    """The overlap formula == the diagonal element of the EXACT
    determinant-space S^2 on the same (possibly broken) determinant
    -- the independent-construction certification, no SCF."""
    from atomterm.detbasis import determinant_basis
    from atomterm.operators import s2_operator_parts
    from atomterm.secondquant import secondquant_matrix

    nao = 4
    mf = _synthetic_uhf(a_set, b_set, nao)
    s2, sigma = uhf_contamination(mf)
    dets = determinant_basis(nao, len(a_set), len(b_set))
    s2o, s2t = s2_operator_parts(nao)
    mat = secondquant_matrix(dets, one_body=s2o, two_body=s2t)
    det = tuple(sorted(tuple(a_set) + tuple(nao + b for b in b_set)))
    i = dets.index(det)
    exact = float(np.real(mat[i, i]))
    assert abs(s2 - exact) < 1e-12
    assert np.all(sigma <= 1 + 1e-12)


def test_uhf_contamination_rotated_pair():
    """A continuously broken pair (angle theta): the formula == the
    exact <psi(theta)|S^2|psi(theta)> with psi the one-body-rotated
    paired determinant, and sigma reads cos(theta)."""
    from scipy.linalg import expm

    from atomterm.detbasis import determinant_basis
    from atomterm.operators import s2_operator_parts
    from atomterm.secondquant import one_body_matrix, secondquant_matrix

    nao = 4
    a_set, b_set = (0, 1, 2), (0, 1)
    dets = determinant_basis(nao, 3, 2)
    s2o, s2t = s2_operator_parts(nao)
    mat = secondquant_matrix(dets, one_body=s2o, two_body=s2t)
    paired = dets.index(tuple(sorted(a_set + tuple(nao + b for b in b_set))))
    psi0 = np.zeros(len(dets), dtype=complex)
    psi0[paired] = 1.0
    # rotates the alpha-1 spin-orbital in the (1, 3) plane
    # (anti-Hermitian generator -> a unitary on the sector)
    g = np.zeros((2 * nao, 2 * nao))
    g[3, 1] = 1.0
    g[1, 3] = -1.0
    rot = one_body_matrix(dets, g)
    for theta in (0.0, 0.6, 1.3):
        psi = expm(theta * rot) @ psi0
        exact = float(np.real(psi.conj() @ mat @ psi))
        c_a = np.eye(nao)
        c_a[:, 1] = (np.cos(theta) * c_a[:, 1]
                     + np.sin(theta) * c_a[:, 3])
        s2, sigma = uhf_contamination(_synthetic_uhf(
            a_set, b_set, nao, c_a=c_a))
        assert abs(s2 - exact) < 1e-12, theta
        assert abs(sigma[0] - 1.0) < 1e-12
        assert abs(sigma[1] - np.cos(theta)) < 1e-12
