"""The reference substrate of the V&V library (scripts/u14e/vv/
reference.py) on the frozen carbon reference: the D2h-tagged wrapper
IS the C1 molecule (certified where it is built), a declaration that
does not describe the reference is refused, the wrapper is built
once, a rung without an SCF record is refused, and a frozen CASCI
does not depend on the symmetry tag (no SCF; tests/carbon_reference.py)."""

import pytest
from pyscf import gto
from pyscf.lib import chkfile

from atomterm.system import AtomConfig

from vv.cases import visible_window
from vv.reference import (Substrate, wrapper_certificate,
                          wrapper_substrate)
from vv.routes import free_casci

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


def test_a_frozen_casci_does_not_depend_on_the_symmetry_tag(vv_substrate):
    # the leg's wrapper check: the same window through the C1
    # substrate and through the D2h-tagged one (observed 1.1e-13)
    space = visible_window(vv_substrate, (("2p", 3),), SECTOR).space
    a = free_casci(vv_substrate, space).states[0].e_ref
    mol, mf = vv_substrate.d2h
    twin = Substrate(mol, mf, vv_substrate.spec, d2h=(mol, mf))
    b = free_casci(twin, space).states[0].e_ref
    assert abs(a - b) < 1e-10
