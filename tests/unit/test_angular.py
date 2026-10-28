"""Tests of the one-electron angular momentum operator construction."""

import numpy as np
import pytest

from pyscf import gto

from atomterm.angular import (
    _complex_generators,
    l2_matrix,
    l_matrices,
    l_matrices_analytic,
    real_to_complex,
)
from atomterm.detbasis import commutator, frobenius_norm

pytestmark = pytest.mark.fast


# a fused-shell basis: general contraction (two columns per shell)
_FUSED = gto.parse("""
C s
 8.0 1.0 0.0
 1.0 0.0 1.0
C p
 1.0 1.0 0.0
 0.5 0.0 1.0
""")


def _fused_mol():
    return gto.M(atom="C 0 0 0", basis={"C": _FUSED}, spin=2, verbose=0)


def test_libcint_matrices_hermitian_imaginary(p_shell_mol):
    l_ops = l_matrices(p_shell_mol)
    for a, mat in l_ops.items():
        assert frobenius_norm(mat.real) < 1e-12
        assert frobenius_norm(mat - mat.conj().T) < 1e-12


def test_libcint_ladder_algebra(p_shell_mol):
    l_ops = l_matrices(p_shell_mol)
    assert frobenius_norm(commutator(l_ops["x"], l_ops["y"]) - 1j * l_ops["z"]) < 1e-12
    assert frobenius_norm(commutator(l_ops["y"], l_ops["z"]) - 1j * l_ops["x"]) < 1e-12
    assert frobenius_norm(commutator(l_ops["z"], l_ops["x"]) - 1j * l_ops["y"]) < 1e-12


def test_lz_p_shell_quantum_numbers(p_shell_mol):
    """Pure p shell: L^2 = 2 * S; Lz eigenvalues = {-1, 0, +1} (+ s block 0)."""
    mol = p_shell_mol
    s = mol.intor("int1e_ovlp")
    p = slice(1, 4)  # normalized p shell
    l2 = l2_matrix(mol)
    assert np.allclose(l2[p, p], 2.0 * s[p, p], atol=1e-10)

    lz = l_matrices(mol)["z"]
    p_blk = lz[p, p]
    vals = np.sort(np.real(np.linalg.eigvals(p_blk)))
    assert np.allclose(vals, [-1.0, 0.0, 1.0], atol=1e-10)


def test_l2_shell_eigenvalues_carbon(carbon):
    """For a single-center spherical basis L^2 chi_nu = l_nu(l_nu+1) chi_nu."""
    import re

    mol = carbon[1]
    s = mol.intor("int1e_ovlp")
    letters = "spdfghi"
    labels = [line.split()[-1] for line in mol.ao_labels()]
    weight = np.array(
        [
            letters.index(re.match(r"\d+([spdfghi])", lbl).group(1))
            * (letters.index(re.match(r"\d+([spdfghi])", lbl).group(1)) + 1)
            for lbl in labels
        ]
    )
    l2 = l2_matrix(mol)
    expected = s * weight[None, :]
    scale = max(1.0, frobenius_norm(expected))
    assert frobenius_norm(l2 - expected) < 1e-9 * scale


def test_analytic_matches_libcint(p_shell_mol, carbon):
    """Two independent constructions of Lx/Ly/Lz must agree."""
    for mol in [p_shell_mol, carbon[1]]:
        ref = l_matrices(mol)
        ana = l_matrices_analytic(mol)
        for a in "xyz":
            assert frobenius_norm(ref[a] - ana[a]) < 1e-9, a


def test_analytic_matches_libcint_on_fused_shells():
    """General contraction (nctr > 1): the analytic path works per
    slot, so fused shells agree with libcint exactly like segmented
    ones (this class crashed the whole-shell-block form)."""
    mol = _fused_mol()
    ref, ana = l_matrices(mol), l_matrices_analytic(mol)
    for a in "xyz":
        assert frobenius_norm(ref[a] - ana[a]) < 1e-9, a


def test_analytic_matches_libcint_sc():
    """Sc/cc-pVDZ carries a fused p shell (nctr = 4): the campaign
    atom the old whole-shell-block form never ran on."""
    mol = gto.M(atom="Sc 0 0 0", basis="ccpvdz", spin=1, verbose=0)
    ref, ana = l_matrices(mol), l_matrices_analytic(mol)
    for a in "xyz":
        assert frobenius_norm(ref[a] - ana[a]) < 1e-12, a  # measured 2.4e-15


# ----------------------------------------------- the tesseral machinery

def test_real_to_complex_unitary():
    for l in range(5):
        u = real_to_complex(l)
        assert np.abs(u.conj().T @ u - np.eye(2 * l + 1)).max() < 1e-12


def test_real_basis_lz_from_transform(synthetic_p_operators):
    """U^(1)-transported Lz in the real (px, py, pz) basis equals the
    Condon-Shortley matrix Lz|px> = i|py>."""
    lz_c = _complex_generators(1)["z"]
    lz_real = real_to_complex(1).conj().T @ lz_c @ real_to_complex(1)
    assert np.abs(lz_real - synthetic_p_operators["z"]).max() < 1e-12


# ------------------------------------------------------ degenerate anchors

def test_l_matrices_vanish_on_s_only_basis():
    """An l=0-only space carries no angular momentum: L and L^2 are
    exactly zero."""
    mol = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    for m in l_matrices(mol).values():
        assert np.abs(m).max() == 0.0
    assert np.abs(l2_matrix(mol)).max() == 0.0
