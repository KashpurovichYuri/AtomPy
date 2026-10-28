"""Fast (SCF-free) pins of atomterm.windows: the generalized
``build_full_mo`` -- window COLUMNS as a first-class input (ticket 29
T1a) -- and the ``Window.columns`` collapse (both representations
assemble through the ONE project-renormalize-assemble treatment; the
c_win sibling builder is gone).

The substrate is synthetic: a gto molecule + per-l S-orthonormal
"canonical" blocks + a hand-built exactly-closed pur core (the fact-20
construction in miniature: one radial replicated on every member). No
SCF anywhere -- sub-second.
"""

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
    """C/cc-pvdz + the synthetic reference. The pur core is the fact-20
    construction in miniature (2 s radials + one p radial replicated on
    every member); the canonical MOs are the pur core theta-rotated
    against its S-orthogonal complement -- canonical core columns carry
    window content and canonical window columns carry pur-core content
    at sin(theta), exactly the physical situation the mo-pur projection
    exists for. (A naive 'two independent l-pure blocks' substrate is
    degenerate: the projected window Gram has a STRUCTURAL null where
    the two 3-dim p subspaces of the 6-dim p space meet.)"""
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    n = mol.nao_nr()
    sec = lm_sectors(mol)
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
    """THE adapter pin: AO-row indices and the explicit AO columns are
    the same window -- byte-equal assemblies (the strangler guarantee
    for every existing AO-index caller)."""
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    act = select_valence_ao(mol, ["2p"])
    a = build_full_mo(mf, list(range(ncore)), act)
    b = build_full_mo(mf, list(range(ncore)), np.eye(mol.nao_nr())[:, act])
    assert np.array_equal(a, b)


def test_build_full_mo_window_columns_on_pur_core():
    """The mo-pur path (ticket 29 T1a): canonical window+virtual columns
    re-orthogonalized against the pur core span. The assembly is
    S-orthonormal, the core block IS the pur core, the window is exactly
    S-orthogonal to it, and the projection has real work to do (the
    planted rotation left pur-core content in the canonical columns --
    the fact-17 class, removed here at the source)."""
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
    """The sibling-death pin: Window.columns routes BOTH representations
    through build_full_mo's one treatment -- for c_win columns
    constructed S-orthogonal to the core (the NO windows) the assembly
    is the no-op it must be, and identical to calling build_full_mo
    directly."""
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
