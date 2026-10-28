"""Unit pins of atomterm.windows: build_full_mo's two representations
(AO rows, column vectors) on a synthetic substrate, plus the STO-3G
degenerate substrates (null core, whole space as window)."""

import numpy as np
import pytest

from pyscf import gto

from atomterm.windows import (Window, build_full_mo, complete_orthonormal,
                              lm_sectors, lowdin, select_valence_ao)

pytestmark = pytest.mark.fast


class _FakeMf:
    """The synthetic reference: an SCF-free substrate provider (mol +
    S-orthonormal "canonical" MOs + their occupations)."""

    def __init__(self, mol, mo, ncore):
        self.mol, self.mo_coeff = mol, mo
        self.mo_occ = np.zeros(mo.shape[1])
        self.mo_occ[:ncore] = 2.0

    def get_ovlp(self):
        return self.mol.intor_symmetric("int1e_ovlp")


def _substrate():
    """C/cc-pvdz + a synthetic reference: per-l S-orthonormal
    "canonical" MOs and a hand-built exactly-closed pur core (one p
    radial replicated on every member), theta-rotated against each
    other so canonical columns carry cross-span content at a planted
    scale (sin 0.3) -- the physical situation the projections exist
    for."""
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    n = mol.nao_nr()
    sec = lm_sectors(mol)
    # NOTE: two independent l-pure blocks would make the projected
    # window Gram STRUCTURALLY singular -- the theta-rotation against
    # the complement is what keeps the substrate non-degenerate
    rows = {l: sorted(r for (ll, _), rs in sec.items() if ll == l
                      for r in rs)
            for l in sorted({ll for ll, _ in sec})}
    blocks = {l: lowdin(np.eye(n)[:, rows[l]], s) for l in rows}

    # the pur core: 2 s radials (s has one member) + one p radial
    # replicated on every member -- E (x) m-set, exactly closed (each
    # member's chain carries the same radial: its first element)
    p_pur = lowdin(np.column_stack(
        [np.eye(n)[sec[(1, m)][0]] for m in range(3)]), s)
    c_pur = np.hstack([blocks[0][:, :2], p_pur])
    ncore = c_pur.shape[1]

    # canonical MOs: [pur | complement] with the p radials theta-rotated
    # against the first complement directions (an S-orthonormal pair
    # stays S-orthonormal under a 2x2 rotation)
    full = complete_orthonormal(c_pur, s)
    theta = 0.3
    c, w = np.cos(theta), np.sin(theta)
    core_can = full[:, :ncore].copy()
    win_can = full[:, ncore:].copy()
    for i in range(3):
        core_can[:, 2 + i], win_can[:, i] = (
            c * full[:, 2 + i] + w * full[:, ncore + i],
            -w * full[:, 2 + i] + c * full[:, ncore + i])
    mo = np.hstack([core_can, win_can])
    return _FakeMf(mol, mo, ncore), mol, s, c_pur, rows


def test_build_full_mo_ao_rows_equal_ao_columns():
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    act = select_valence_ao(mol, ["2p"])
    a = build_full_mo(mf, list(range(ncore)), act)
    b = build_full_mo(mf, list(range(ncore)), np.eye(mol.nao_nr())[:, act])
    assert np.array_equal(a, b)


def test_build_full_mo_window_columns_on_pur_core():
    """Canonical columns on the pur core: orthonormal assembly, the
    pur block intact, the planted cross-span leak removed."""
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    mo = build_full_mo(mf, list(range(ncore)), mf.mo_coeff[:, ncore:],
                       c_core=c_pur)
    assert mo.shape == (mol.nao_nr(), mol.nao_nr())
    err = np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max()
    assert err < 1e-12
    # the core block is the pur core (Loewdin of an S-orthonormal block
    # is itself, to machine precision)
    assert np.allclose(mo[:, :ncore], c_pur, atol=1e-12)
    # the window columns are exactly S-orthogonal to the core span
    w = mo[:, ncore:ncore + 3]
    assert np.abs(w.conj().T @ s @ c_pur).max() < 1e-12
    # and the projection did real work: the canonical window columns
    # carried pur-core-span content at the planted scale
    w_can = mf.mo_coeff[:, ncore:ncore + 3]
    leak = np.linalg.norm(w_can.conj().T @ s @ c_pur)
    assert leak > 1e-3, leak      # the fact-17-class component existed
    assert np.linalg.norm(w.conj().T @ s @ c_pur) < 1e-12  # ...and is gone


def test_window_columns_assemble_through_build_full_mo():
    """Window.columns(c_win) == build_full_mo directly; columns already
    orthogonal to the core pass through unchanged."""
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    # an S-orthonormal window orthogonal to the pur core: the d block
    c_win = lowdin(np.eye(mol.nao_nr())[:, rows[2]], s)[:, :3]
    assert np.abs(c_win.conj().T @ s @ c_pur).max() < 1e-12
    win = Window("synthetic no", None, c_win, False)
    mo = win.columns(mf, ncore, c_core=c_pur)
    closed = list(range(ncore))
    assert np.array_equal(mo, build_full_mo(mf, closed, c_win, c_core=c_pur))
    err = np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max()
    assert err < 1e-12
    assert np.allclose(mo[:, ncore:ncore + 3], c_win, atol=1e-12)


# ---------------------------------------------------------------------
# degenerate substrates (STO-3G): null core, the whole AO space as the
# window -- rotation-closed by construction, so the certificates are
# exact; any deviation is a bug.

def test_h_sto3g_null_core_full_window():
    """Empty core + full window: the zero-column paths and an
    exactly-zero closure certificate."""
    from pyscf import scf
    from atomterm.operators import window_closure

    mol = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    mf = scf.RHF(mol).run(verbose=0)
    n = mol.nao_nr()
    s = mol.intor_symmetric("int1e_ovlp")
    mo = build_full_mo(mf, [], list(range(n)))
    assert mo.shape == (n, n)
    assert np.abs(mo.conj().T @ s @ mo - np.eye(n)).max() < 1e-12
    cl = window_closure(mol, mo)
    assert cl["lz"] < 1e-14 and cl["lperp"] < 1e-14


def test_c_sto3g_null_core_full_window():
    """Null core, whole AO space as window (CAS(6,5), sector (4,2)):
    certificates exactly zero, [H, Lz]/[H, L2] at machine precision,
    3P ground."""
    from pyscf import scf
    from atomterm.operators import build_operators, window_closure

    mol = gto.M(atom="C 0 0 0", basis="sto-3g", spin=2, verbose=0)
    mf = scf.ROHF(mol).run(verbose=0)
    n = mol.nao_nr()
    mo = build_full_mo(mf, [], list(range(n)))
    cl = window_closure(mol, mo)
    assert cl["lz"] < 1e-13 and cl["lperp"] < 1e-13

    ops = build_operators(mol, mf, mo, 0, 4, 2)
    h, lz, l2, s2 = ops["H"], ops["Lz"], ops["L2"], ops["S2"]
    assert np.linalg.norm(h @ lz - lz @ h) < 1e-11
    assert np.linalg.norm(h @ l2 - l2 @ h) < 1e-11
    ev, evec = np.linalg.eigh(h)
    v = evec[:, 0]
    assert abs(float(np.real(v.conj() @ l2 @ v)) - 2.0) < 1e-8  # 3P ground
    assert abs(float(np.real(v.conj() @ s2 @ v)) - 2.0) < 1e-8
