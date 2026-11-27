"""The reference substrate of the V&V library (scripts/u14e/vv/
reference.py) on the frozen carbon reference: the D2h-tagged wrapper
IS the C1 molecule (overlap and one-electron integrals certified
where it is built, the two-electron ones compared here), a
declaration that does not describe the reference is refused, the
wrapper is built once, a rung without an SCF record is refused, and
the free routes refuse the symmetry-tagged substrate -- there PySCF
pins the FCI to one irrep (no SCF; tests/carbon_reference.py)."""

import numpy as np
import pytest
from pyscf import ao2mo, gto
from pyscf.lib import chkfile

from atomterm.system import AtomConfig

from vv.cases import visible_window
from vv.reference import (Substrate, wrapper_certificate,
                          wrapper_substrate)
from vv.routes import free_casci, free_casscf, free_sa

pytestmark = pytest.mark.fast

SECTOR = (2, 0)


def test_the_wrapper_is_the_c1_molecule(vv_substrate, carbon):
    mol, mf = vv_substrate.mol, vv_substrate.mf
    mol_d2h, mf_d2h = wrapper_substrate(mol, mf, carbon[0])
    assert mol_d2h.groupname == "D2h" and mol.groupname == "C1"
    cert = wrapper_certificate(mol, mf, mol_d2h, mf_d2h)
    assert cert["identical"] and cert["dev_s"] == 0.0 and cert["dev_h"] == 0.0
    # the record rides along: the wrapper's orbitals ARE the C1 ones
    assert mf_d2h.mo_coeff is mf.mo_coeff
    assert mf_d2h.e_tot == mf.e_tot and mf_d2h.converged


def test_a_declaration_that_is_not_the_reference_is_refused(vv_substrate):
    cfg = AtomConfig(basis="sto-3g")        # the reference is cc-pVDZ
    with pytest.raises(ValueError, match="not the C1 molecule"):
        wrapper_substrate(vv_substrate.mol, vv_substrate.mf, cfg)


def test_the_certificate_reads_the_difference(vv_substrate, carbon):
    mol2 = gto.M(atom="C 0 0 0", basis="sto-3g", spin=2, symmetry="D2h",
                 verbose=0)
    from pyscf import scf
    cert = wrapper_certificate(vv_substrate.mol, vv_substrate.mf, mol2,
                               scf.ROHF(mol2))
    assert not cert["identical"] and cert["dev_s"] == float("inf")


def test_the_substrate_builds_its_wrapper_once_and_needs_a_declaration(
        vv_substrate):
    d1 = vv_substrate.d2h
    assert vv_substrate.d2h is d1                       # cached
    assert d1[0].groupname == "D2h"
    with pytest.raises(ValueError, match="no case declaration"):
        _ = Substrate(vv_substrate.mol, vv_substrate.mf).d2h
    ready = Substrate(vv_substrate.mol, vv_substrate.mf, d2h=d1)
    assert ready.d2h is d1                              # supplied as is


def test_a_rung_without_an_scf_record_is_refused(tmp_path):
    from vv.reference import load_rung
    chk = str(tmp_path / "empty.chk")
    chkfile.save_mol(gto.M(atom="H", basis="sto-3g", spin=1, verbose=0), chk)
    with pytest.raises(RuntimeError, match="no converged SCF record"):
        load_rung(chk)


def test_the_wrapper_carries_the_c1_two_electron_integrals(vv_substrate):
    # the overlap and the one-electron integrals are certified where
    # the wrapper is built; the two-electron ones over the core +
    # window block are compared directly (observed 0.0)
    space = visible_window(vv_substrate, (("2p", 3),), SECTOR).space
    c = space.mo[:, :space.ncore + space.ncas]
    mol_d2h, _mf = vv_substrate.d2h
    a = ao2mo.restore(1, ao2mo.kernel(vv_substrate.mol, c), c.shape[1])
    b = ao2mo.restore(1, ao2mo.kernel(mol_d2h, c), c.shape[1])
    assert np.abs(a - b).max() == 0.0


def test_a_free_route_refuses_a_symmetry_tagged_substrate(vv_substrate):
    # on the tagged substrate PySCF's symmetry-adapted SCF makes a
    # plain CASCI use its symmetry FCI solver, pinned to one irrep --
    # another calculation (2e-02 Eh above the sector ground on a U14e
    # window; on carbon the 3P ground happens to sit in that block)
    mol, mf = vv_substrate.d2h
    twin = Substrate(mol, mf, vv_substrate.spec, d2h=(mol, mf))
    space = visible_window(vv_substrate, (("2p", 3),), SECTOR).space
    assert free_casci(vv_substrate, space).certificate["residual"] < 1e-12
    for route in (free_casci, free_casscf,
                  lambda s, w: free_sa(s, w, nroots=3)):
        with pytest.raises(ValueError, match="symmetry-tagged substrate"):
            route(twin, space)
