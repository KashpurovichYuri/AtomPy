"""The PySCF-mrpt anchors at example scale: the retained O2/H14
SC-NEVPT2 values (examples/mrpt/01) as the layer-(b) anchor, the
multi-root root-selection read (stretched H2, the near-degenerate
pair), and the SA discipline (example 41's guard, kernelled)."""

import numpy as np
import pytest
from pyscf import gto, scf

from atomterm.perturb import nevpt2_rows, readout_mc

pytestmark = pytest.mark.fast

O2_REF = -0.169785157128082   # examples/mrpt/01's retained value
H14_REF = -0.094164359938171


def _o2():
    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None
    mol.atom = [["O", (0.0, 0.0, 0.0)], ["O", (0.0, 0.0, 1.207)]]
    mol.basis = "6-31g"
    mol.spin = 2
    mol.build()
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-20
    mf.scf()
    return mol, mf


def test_o2_anchor():
    """O2/6-31G CASCI(6,8): our class plumbing sums to the example's
    retained NEVPT2 value (fresh agreement 3.2e-7 across versions) and
    to the kernel itself (the plumbing-vs-kernel twin)."""
    from pyscf.mrpt import NEVPT

    mol, mf = _o2()
    mc = readout_mc(mol, mf, mf.mo_coeff, 6, (5, 3), 1)
    rows = nevpt2_rows(mol, mf, mf.mo_coeff, 6, (5, 3), 1)
    assert len(rows) == 1
    row = rows[0]
    assert row.e_corr == pytest.approx(O2_REF, abs=1e-6)
    # the engine's own run-to-run wobble: PySCF's kernel (and every
    # construction through its big contractions) wobbles ~3e-7 Eh
    # across runs (fresh observed: kernel -0.16978547 / -0.16978542 /
    # -0.16978512; plumbing within the same class) -- the twin pin
    # sits at that floor; the small-CAS fixture case stays 1e-10
    assert row.e_corr == pytest.approx(
        float(NEVPT(mc, root=0, density_fit=False).kernel()), abs=1e-6)
    assert len(row.classes) == 8
    assert row.e_total == pytest.approx(row.e_ref + row.e_corr,
                                        abs=1e-14)


def _h2():
    """Stretched H2/STO-3G (1.4 A): the shared small substrate."""
    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None
    mol.atom = [["H", (0.0, 0.0, 0.0)], ["H", (0.0, 0.0, 1.4)]]
    mol.basis = "sto-3g"
    mol.build()
    mf = scf.RHF(mol)
    mf.scf()
    return mol, mf

def test_h14_anchor():
    """H14/STO-3G CASCI(8,10): the second retained anchor, a different
    system class (closed-shell chain, larger active space)."""
    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None
    mol.atom = [["H", (0.0, 0.0, 2.0 * k + d)] for k in range(7)
                for d in (0.0, 0.8)]
    mol.basis = "sto-3g"
    mol.build()
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-20
    mf.scf()
    rows = nevpt2_rows(mol, mf, mf.mo_coeff, 8, (5, 5), 1)
    assert rows[0].e_corr == pytest.approx(H14_REF, abs=1e-6)


def test_two_root_h2_root_selection():
    """Stretched H2/STO-3G CASCI(2,2), two roots: root 1's row carries
    its own CASCI total, and the class plumbing matches the kernel for
    the SECOND root (the multi-root path, not just the ground)."""
    from pyscf.mrpt import NEVPT

    mol, mf = _h2()
    mc = readout_mc(mol, mf, mf.mo_coeff, 2, (1, 1), 2)
    rows = nevpt2_rows(mol, mf, mf.mo_coeff, 2, (1, 1), 2, roots=[0, 1])
    assert [r.root for r in rows] == [0, 1]
    assert rows[1].e_ref > rows[0].e_ref
    assert rows[1].e_corr == pytest.approx(
        float(NEVPT(mc, root=1, density_fit=False).kernel()), abs=1e-10)


def test_sa_kernelled_object_refuses():
    """The example-41 discipline, kernelled: after a state-average
    CASSCF, the SA solver object itself refuses at NEVPT -- the
    multi-root CASCI readout is the required form."""
    from pyscf import mcscf
    from pyscf.mrpt import NEVPT

    mol, mf = _h2()
    mc = mcscf.CASSCF(mf, 2, 2).state_average_(weights=(0.5, 0.5))
    mc.kernel()
    with pytest.raises(RuntimeError):
        NEVPT(mc, root=0).kernel()
