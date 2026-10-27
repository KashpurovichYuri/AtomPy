"""The cell layer: key addressing, the store contract, and the CI
ordering bridge.

The bridge is the load-bearing piece: PySCF FCI vectors and
``detbasis`` enumerations order determinants differently for norb >= 4
with >= 2 electrons in a block (they coincide on carbon's norb=3
window -- why the difference went unnoticed). The certificate test
verifies the mapping through the eigen-equation against an
independently built Slater-Condon H: generic integrals admit no
permutation automorphism, so a wrong mapping cannot satisfy it.
"""

import numpy as np
import pytest

from atomterm.cells import (
    CellData,
    CellKey,
    cell_chk,
    ci_in_detbasis_order,
    load_cell,
    pyscf_ci_permutation,
)
from atomterm.detbasis import determinant_basis, slater_condon_matrix

pytestmark = pytest.mark.fast


# --------------------------------------------------------------- the key

def test_address_round_trip():
    keys = [
        CellKey("C", "2p", (1, 1)),
        CellKey("C", "2p", (1, 1), fold="raw"),
        CellKey("c", "2p", (2, 0), method="sa", mask=None, substrate="pur"),
        CellKey("Ti", "3d4s", (3, 1), method="sa", start="mo",
                mask="radial", substrate="pur", fold="avg"),
        CellKey("Ti", "3d4s", (2, 2), method="casscf", start="ao",
                mask="virtmix"),
        CellKey("Sc", "3d4s", (1, 1), method="casscf", mask="frozen",
                fold="frac"),
        CellKey("U", "3f2d", (4, 2), method="sa", start="ao"),
    ]
    for k in keys:
        assert CellKey.from_address(k.address) == k
        assert k.address.count("_") == 7
    assert keys[0].address == "c_2p_1-1_casci_ao_-_pol_avg"
    assert keys[3].address == "ti_3d4s_3-1_sa_mo_radial_pur_avg"


def test_key_validation_and_canonicalization():
    # the pur-fold degeneracy rule: raw on pur canonicalizes to avg
    k = CellKey("C", "2p", (1, 1), substrate="pur", fold="raw")
    assert k.fold == "avg"
    # frac is a different REFERENCE -- meaningless on the purified span
    with pytest.raises(ValueError, match="frac"):
        CellKey("C", "2p", (1, 1), substrate="pur", fold="frac")
    # a mask needs an optimizer method
    with pytest.raises(ValueError, match="mask"):
        CellKey("C", "2p", (1, 1), method="casci", mask="radial")
    for bad in [dict(method="fci"), dict(start="aos"), dict(mask="relax"),
                dict(substrate="polarized"), dict(fold="averaged")]:
        with pytest.raises(ValueError, match="not in"):
            CellKey("C", "2p", (1, 1), **bad)
    with pytest.raises(ValueError, match="sector"):
        CellKey("C", "2p", (0, 1))
    with pytest.raises(ValueError, match="shell label"):
        CellKey("C", "3d_4s", (1, 1))
    # the window label parses to shell selectors
    assert CellKey("Ti", "3d4s", (3, 1)).shell_patterns == ["3d", "4s"]
    assert CellKey("C", "2p", (1, 1)).shell_patterns == ["2p"]


# ------------------------------------------------------------- the bridge

SECTORS = [(3, 1, 1), (4, 2, 1), (5, 3, 1), (6, 3, 2), (6, 4, 2)]


def test_ci_permutation_is_a_permutation():
    from math import comb

    for norb, na, nb in SECTORS:
        perm = pyscf_ci_permutation(norb, na, nb)
        assert len(perm) == comb(norb, na) * comb(norb, nb)
        assert np.array_equal(np.sort(perm), np.arange(len(perm)))


def test_orders_coincide_only_in_trivial_sectors():
    """The empirical fact that kept the bug invisible: carbon's norb=3
    window forces pyscf's integer-ascending string order to coincide
    with detbasis' lexicographic order; norb >= 4 breaks it."""
    assert np.array_equal(pyscf_ci_permutation(3, 1, 1), np.arange(9))
    assert not np.array_equal(pyscf_ci_permutation(4, 2, 1), np.arange(24))


def _random_eri(norb, rng):
    """A generic chemist-notation eri with the full 8-fold symmetry
    (average over the orbit under the pair transpositions) -- PySCF's
    contraction assumes it."""
    a = rng.standard_normal((norb,) * 4)
    seen, fr = [], [a]
    while fr:
        x = fr.pop()
        if any(np.array_equal(x, s) for s in seen):
            continue
        seen.append(x)
        for p in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
            fr.append(np.ascontiguousarray(x.transpose(p)))
    eri = sum(seen) / len(seen)
    for p in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
        assert np.allclose(eri, eri.transpose(p))
    return eri


@pytest.mark.parametrize("norb,na,nb", SECTORS)
def test_bridge_certified_by_eigen_equation(norb, na, nb):
    """THE bridge certificate: a PySCF-solved CI vector, mapped into
    detbasis order, satisfies the eigen-equation of the independently
    built Slater-Condon Hamiltonian. No ordering assumption anywhere --
    generic integrals admit no permutation automorphism."""
    from pyscf import fci
    from atomterm.operators import spin_block_one_body, spin_block_two_body

    rng = np.random.default_rng(20261027 + 10 * norb + na)
    h1 = rng.standard_normal((norb, norb))
    h1 = (h1 + h1.T) / 2
    eri = _random_eri(norb, rng)

    e_fci, ci = fci.direct_spin1.kernel(h1, eri, norb, (na, nb),
                                        conv_tol=1e-13)
    dets = determinant_basis(norb, na, nb)
    # the owners' spin-orbital expansion; detbasis eats the PHYSICIST
    # arrangement of the same eri (fact 3), pyscf the chemist one
    ham = slater_condon_matrix(
        dets, one_body=spin_block_one_body(h1),
        two_body=spin_block_two_body(eri.transpose(0, 2, 1, 3)))

    w = ci_in_detbasis_order(ci, norb, na, nb)
    e_mine = float(np.real(w.conj() @ ham @ w))
    assert abs(e_mine - e_fci) < 1e-9
    resid = ham @ w - e_fci * w
    assert np.linalg.norm(resid) / np.linalg.norm(ham) < 1e-10


def test_configuration_summary_uses_the_bridge(tmp_path):
    """The consumer fix: a pyscf-ordered unit amplitude must label the
    determinant at ITS string address, not dets[k]."""
    from pyscf import gto
    from pyscf.fci import cistring
    from atomterm.analysis import configuration_summary

    mol = gto.M(atom="C 0 0 0", basis="cc-pvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    # a mixed-l window (2s + 2p) so column identity matters
    cols = [1, 2, 3, 4]
    c_win = np.eye(mol.nao_nr())[:, cols]
    na, nb = 2, 1
    # pyscf address 5 of the (2,1)/(4) sector -> its determinant
    a_strs = cistring.make_strings(range(4), na)
    b_strs = cistring.make_strings(range(4), nb)
    ia, ib = divmod(5, len(b_strs))
    sa, sb = a_strs[ia], b_strs[ib]
    occ_a = [i for i in range(4) if sa >> i & 1]
    occ_b = [i for i in range(4) if sb >> i & 1]
    # expected label: the l-mass of exactly those columns
    from atomterm.drivers import shell_weights
    names, wl = shell_weights(mol, c_win, s)
    occ_cols = sorted(occ_a + occ_b)
    mass = {n: sum(wl[c, k] for c in occ_cols)
            for k, n in enumerate(names)}
    expected = "".join(f"{n}^{mass[n]:.2f}" for n in names if mass[n] > 0.01
                       ) + " (|c|=1.00)"
    v = np.zeros(len(a_strs) * len(b_strs))
    v[5] = 1.0
    [label] = configuration_summary(mol, s, [v], c_win, na, nb, top=1)
    assert label == expected


# -------------------------------------------------------- the store layer

def test_load_cell_synthetic_roundtrip(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1))  # any atom label; synthetic store
    chk = cell_chk(key, tmp_path)
    chk.parent.mkdir(parents=True)
    rng = np.random.default_rng(7)
    mo = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
    energies = np.array([-1.0, -0.5, -0.25])
    ci = [rng.standard_normal(9), rng.standard_normal(9)]
    for name, val in (("energies", energies),
                      ("ci", np.asarray(ci, dtype=object)),
                      ("core_energy", 5.25), ("mo_coeff", mo)):
        lib.chkfile.save(str(chk), f"atomterm/{key.address}/{name}", val)

    data = load_cell(key, tmp_path)
    assert isinstance(data, CellData)
    assert data.key == key
    assert np.array_equal(data.mo_coeff, mo)
    assert np.array_equal(data.energies, energies)
    assert data.core_energy == 5.25
    assert len(data.ci) == 2
    assert np.allclose(data.ci[1], ci[1])


def test_load_cell_prefers_sa_roots_and_reports_missing(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1), method="sa")
    chk = cell_chk(key, tmp_path)
    chk.parent.mkdir(parents=True)
    for name in (key.address, f"{key.address}_roots"):
        lib.chkfile.save(str(chk), f"atomterm/{name}/energies",
                         np.array([-2.0 if "roots" in name else -1.0]))
        lib.chkfile.save(str(chk), f"atomterm/{name}/ci",
                         np.asarray([np.zeros(9)], dtype=object))
        lib.chkfile.save(str(chk), f"atomterm/{name}/core_energy", 1.0)
    lib.chkfile.save(str(chk), f"atomterm/{key.address}/mo_coeff",
                     np.eye(6))
    assert load_cell(key, tmp_path).energies[0] == -2.0

    with pytest.raises(FileNotFoundError, match="run it first"):
        load_cell(CellKey("Xe", "5p", (2, 0)), tmp_path)
