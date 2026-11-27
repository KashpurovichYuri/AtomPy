"""The active spaces of the V&V library (scripts/u14e/vv/windows.py)
on the frozen carbon reference: the visible window IS the production
build (core + window elementwise, the complement by span) and closes
under rotations, the whole non-core space is an exact-by-construction
window (no complement at all), a wrong declared width is refused, the
matrix a space carries cannot be written to, a broken matrix or a
partition that misses electrons is refused, and the diagnostic
readings (survival norms, dominant labels, the s column, the core
pair count) are exact on constructed inputs (no SCF;
tests/carbon_reference.py)."""


import numpy as np
import pytest

from atomterm.core import purified_core_orbitals
from atomterm.windows import build_full_mo, select_shell_ao

from vv.cases import core_pairs, visible_window
from vv.windows import (ActiveSpace, dominant_labels, pair_mo, s_column,
                        survival_norms, window_from_shells)

pytestmark = pytest.mark.fast

SECTOR = (2, 0)


@pytest.fixture(scope="module")
def core(vv_substrate):
    """The purified core of the frozen reference."""
    return purified_core_orbitals(vv_substrate.mol, vv_substrate.mf, 2)[0]


def test_the_visible_window_is_the_production_build(vv_substrate, core):
    built = window_from_shells(vv_substrate.mol, vv_substrate.mf, (("2p", 3),),
                               core, SECTOR)
    space = built.space
    mol, mf = vv_substrate.mol, vv_substrate.mf
    assert (space.ncore, space.ncas, space.sector) == (2, 3, SECTOR)
    assert "visible construction == production build_full_mo" in built.report[0]
    # core + window elementwise against the builder, independently here
    prod = build_full_mo(mf, [0, 1], select_shell_ao(mol, "2p"), c_core=core)
    assert np.abs(space.mo[:, :5] - prod[:, :5]).max() < 1e-12
    assert space.certify(mf) < 1e-13
    cl = space.closure(mol)
    assert cl["lz"] < 1e-13 and cl["lperp"] < 1e-13
    assert np.abs(space.core - core).max() == 0.0
    assert space.active.shape == (mol.nao, 3) and space.virtual.shape == (mol.nao, 9)


def test_the_whole_non_core_space_is_a_window(vv_substrate, core):
    # exact by construction: everything outside the core active (the
    # s radial the core leaves, the whole p and d shells) -- no
    # complement, every closure meter at roundoff
    mol, mf = vv_substrate.mol, vv_substrate.mf
    shells = (("3s", 1), ("2p", 3), ("3p", 3), ("3d", 5))
    space = window_from_shells(mol, mf, shells, core, SECTOR).space
    assert space.ncas == mol.nao - 2
    assert space.virtual.shape == (mol.nao, 0)
    assert space.certify(mf) < 1e-12
    cl = space.closure(mol)
    assert max(cl["lz"], cl["lperp"]) < 1e-12


def test_a_wrong_declared_width_is_refused(vv_substrate, core):
    with pytest.raises(ValueError, match="found 3 AOs, expected 4"):
        window_from_shells(vv_substrate.mol, vv_substrate.mf, (("2p", 4),),
                           core, SECTOR)


def test_the_matrix_a_space_carries_cannot_be_written(vv_substrate, core):
    space = window_from_shells(vv_substrate.mol, vv_substrate.mf, (("2p", 3),),
                               core, SECTOR).space
    with pytest.raises(ValueError, match="read-only"):
        space.mo[0, 0] = 1.0


def test_a_broken_matrix_or_partition_is_refused(vv_substrate, core):
    space = window_from_shells(vv_substrate.mol, vv_substrate.mf, (("2p", 3),),
                               core, SECTOR).space
    mo = np.array(space.mo)
    mo[:, 0] *= 1.1
    with pytest.raises(ValueError, match="not S-orthonormal"):
        ActiveSpace(mo, 2, 3, SECTOR).certify(vv_substrate.mf)
    with pytest.raises(ValueError, match="accounts for .* electrons"):
        ActiveSpace(np.array(space.mo), 3, 3, SECTOR).certify(vv_substrate.mf)


def test_on_keeps_the_partition_on_another_orbital_set(vv_substrate, core):
    space = window_from_shells(vv_substrate.mol, vv_substrate.mf, (("2p", 3),),
                               core, SECTOR).space
    other = space.on(space.mo * np.array([1.0, -1.0] + [1.0] * 12))
    assert (other.ncore, other.ncas, other.sector) == (2, 3, SECTOR)
    assert other.certify(vv_substrate.mf) < 1e-13
    with pytest.raises(ValueError, match="read-only"):
        other.mo[0, 0] = 1.0


# ---------------------------------------------------------- the diagnostics

def test_survival_norms_are_exact_on_constructed_columns(vv_substrate, core):
    s = vv_substrate.mf.get_ovlp()
    p = np.eye(vv_substrate.mol.nao)[:, select_shell_ao(vv_substrate.mol, "2p")]
    inside = core[:, 0] + 0.5 * core[:, 1]
    act = np.column_stack([p[:, 0], inside])
    norms = survival_norms(act, core, s)
    assert abs(norms[0] - 1.0) < 1e-12          # core-orthogonal by symmetry
    # a remainder norm is sqrt(diag - projection): its floor is the
    # square root of the cancellation's roundoff (observed 2.1e-08)
    assert norms[1] < 1e-7                      # inside the core's span


def test_dominant_labels_and_the_s_column(vv_substrate, core):
    space = window_from_shells(vv_substrate.mol, vv_substrate.mf, (("2p", 3),),
                               core, SECTOR).space
    assert dominant_labels(space.active, vv_substrate.mol) == ["2px", "2py", "2pz"]
    with pytest.raises(ValueError, match="s sector is 0 columns"):
        s_column(space.active, vv_substrate.mol)


def test_pair_mo_finds_the_closed_s_pair_outside_the_core(vv_substrate):
    # the frozen carbon reference: 1s, 2s closed; outside the 1-pair
    # core the only s-dominated occupied column is the 2s
    j, col = pair_mo(vv_substrate.mol, vv_substrate.mf, 1)
    assert j == 1
    assert np.abs(col - np.asarray(vv_substrate.mf.mo_coeff)[:, 1]).max() == 0.0
    with pytest.raises(ValueError, match="closed block has 2"):
        pair_mo(vv_substrate.mol, vv_substrate.mf, 3)


# ------------------------------------------------------- the case-level pieces

def test_core_pairs_follow_the_electron_count(vv_substrate):
    assert core_pairs(vv_substrate, (2, 0)) == 2
    assert core_pairs(vv_substrate, (1, 1)) == 2
    assert core_pairs(vv_substrate, (4, 0)) == 1
    with pytest.raises(ValueError, match="not a positive even count"):
        core_pairs(vv_substrate, (3, 0))
    with pytest.raises(ValueError, match="not a positive even count"):
        core_pairs(vv_substrate, (4, 2))


def test_visible_window_builds_on_the_purified_core(vv_substrate):
    built = visible_window(vv_substrate, (("2p", 3),), SECTOR)
    assert built.space.ncore == 2 and len(built.report) == 1
    assert built.space.certify(vv_substrate.mf) < 1e-13
