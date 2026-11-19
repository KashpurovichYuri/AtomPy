"""perturb.py: the PySCF-mrpt adapter -- the readout, the stored-roots
refusal (drift + degenerate-block gauge), the class breakdown against
the kernel twin, and the row guards.  All on the FROZEN carbon fixture
(NO SCF): the C 2p2 Sz=+1 CAS(3o,(2,0)) -- three determinants, the
exactly degenerate 3P block on the pur core."""

import numpy as np
import pytest

from atomterm.perturb import (
    ClassRow,
    PT2Row,
    assert_stored_roots,
    nevpt2_rows,
    readout_mc,
)

pytestmark = pytest.mark.fast

# the C/cc-pVDZ partition: core = 1s2s (the electron-count derivation),
# active = the 2p shell (AO-row selectors), Sz=+1 sector (na, nb) = (2, 0)
NCAS, NELECAS, NROOTS = 3, (2, 0), 3


@pytest.fixture(scope="module")
def block(carbon):
    """The 3-root readout on the campaign's certified construction:
    the pur (rotation-closed) core + the raw AO 2p window, assembled
    S-orthonormal by ``build_full_mo``."""
    _, mol, mf, act_idx = carbon
    from atomterm.core import purified_core_orbitals
    from atomterm.windows import build_full_mo

    c_core, _ = purified_core_orbitals(mol, mf, 2)
    mo_full = build_full_mo(mf, [0, 1], act_idx, c_core=c_core)
    mc = readout_mc(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    return mol, mf, mc, mo_full


def test_readout_reproduces_the_degenerate_p2_block(block):
    """The Sz=+1 p2 sector holds ONLY the 3P states (ML = +1, 0, -1):
    three determinants, one exact energy -- the spread certifies the
    pur core (the raw core's polarization splits it; the next test)."""
    mol, mf, mc, mo_full = block
    e = np.atleast_1d(np.asarray(mc.e_cas, dtype=float))
    assert e.shape == (3,)
    assert np.max(e) - np.min(e) < 1e-9


def test_raw_core_field_splits_the_block(carbon):
    """The raw-core anisotropy detector at unit scale: the same readout
    on the RAW canonical core polarizes the 3P block (fresh observed
    4.4e-5 Eh -- the declared class; the pur core above cures it)."""
    _, mol, mf, act_idx = carbon
    from atomterm.windows import build_full_mo

    mo_raw = build_full_mo(mf, [0, 1], act_idx)
    mc = readout_mc(mol, mf, mo_raw, NCAS, NELECAS, NROOTS)
    e = np.atleast_1d(np.asarray(mc.e_cas, dtype=float))
    assert 1e-6 < np.max(e) - np.min(e) < 1e-3


def test_root_energies_follow_the_store_convention(block):
    """Root totals = e_cas + the common core part; root 0 reproduces
    the kernel's own total."""
    mol, mf, mc, mo_full = block
    from atomterm.perturb import _root_energies

    e = _root_energies(mc)
    e_tot0 = float(np.atleast_1d(np.asarray(mc.e_tot)).ravel()[0])
    assert abs(e[0] - e_tot0) < 1e-11


def test_stored_roots_assertion_passes_on_truth_and_refuses_drift(block):
    """The assertion accepts the solution it read out and refuses any
    drift: energy tamper and vector tamper both raise."""
    mol, mf, mc, mo_full = block
    from atomterm.perturb import _root_energies

    e_ref = _root_energies(mc).copy()
    ci_ref = [np.asarray(c).copy() for c in mc.ci]
    assert_stored_roots(mc, ref_energies=e_ref, ref_ci=ci_ref)  # truth
    # the container shape is not the assertion: the store keeps flat
    # vectors, the reader returns (n, 1) columns (observed on the
    # campaign store)
    assert_stored_roots(mc, ref_ci=[c.ravel() for c in ci_ref])

    e_bad = e_ref.copy()
    e_bad[1] += 1e-9
    with pytest.raises(ValueError):
        assert_stored_roots(mc, ref_energies=e_bad)

    ci_bad = [c.copy() for c in ci_ref]
    ci_bad[2] = ci_bad[2] * 0.9 + 1e-3
    ci_bad[2][0] += 0.5
    with pytest.raises(ValueError):
        assert_stored_roots(mc, ref_ci=ci_bad)


def test_stored_roots_assertion_is_gauge_tolerant_in_the_block(block):
    """A lawful gauge redraw of the exactly degenerate block (a random
    unitary inside it) passes -- the joint projector is the
    invariant; zeroing one root's vector does not."""
    mol, mf, mc, mo_full = block
    from atomterm.perturb import _root_energies

    e_ref = _root_energies(mc).copy()
    ci_ref = [np.asarray(c).copy() for c in mc.ci]
    rng = np.random.default_rng(3)
    q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    rotated = [sum(q[k, j] * ci_ref[j] for j in range(3)) for k in range(3)]
    assert_stored_roots(mc, ref_energies=e_ref, ref_ci=rotated)  # gauge

    zapped = [c.copy() for c in ci_ref]
    zapped[1] = np.zeros_like(zapped[1])
    with pytest.raises(ValueError):
        assert_stored_roots(mc, ref_energies=e_ref, ref_ci=zapped)


def test_rows_carry_the_full_class_breakdown_and_match_the_kernel(block):
    """The row: eight classes, nonnegative norms, the class sum as
    e_corr -- pinned against NEVPT(...).kernel() on a FRESH readout
    (the kernel's loop writes the canonicalized CI into the shared
    ``mc.ci``; the fixture must stay untouched)."""
    mol, mf, mc, mo_full = block
    rows = nevpt2_rows(mol, mf, mo_full, NCAS, NELECAS, NROOTS,
                       roots=[0])
    assert len(rows) == 1
    row = rows[0]
    assert isinstance(row, PT2Row)
    assert row.variant == "nevpt2-sc"
    assert len(row.classes) == 8
    assert [c.label for c in row.classes] == [
        "Sr", "Si", "Sijrs", "Sijr", "Srsi", "Srs", "Sij", "Sir"]
    assert all(isinstance(c, ClassRow) and c.norm >= -1e-12
               for c in row.classes)
    assert row.e_corr == pytest.approx(
        sum(c.energy for c in row.classes), abs=1e-14)
    assert row.e_total == pytest.approx(row.e_ref + row.e_corr,
                                        abs=1e-14)
    from pyscf.mrpt import NEVPT

    mc_twin = readout_mc(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    assert row.e_corr == pytest.approx(
        float(NEVPT(mc_twin, root=0, density_fit=False).kernel()),
        abs=1e-10)


def test_row_refs_exercise_the_assertion_path(block):
    """nevpt2_rows with refs re-reads the solution and passes; a
    tampered store refuses the whole call (fail-fast before any PT2)."""
    from atomterm.perturb import _root_energies

    mol, mf, mc, mo_full = block
    e_ref = _root_energies(mc).copy()
    ci_ref = [np.asarray(c).copy() for c in mc.ci]
    rows = nevpt2_rows(mol, mf, mo_full, NCAS, NELECAS, NROOTS,
                       roots=[1], ref_energies=e_ref, ref_ci=ci_ref)
    assert rows[0].root == 1
    e_bad = e_ref.copy()
    e_bad[0] += 1e-6
    with pytest.raises(ValueError):
        nevpt2_rows(mol, mf, mo_full, NCAS, NELECAS, NROOTS,
                    roots=[1], ref_energies=e_bad)


def test_variant_whitelist_refuses_unknown(block):
    mol, mf, _, mo_full = block
    with pytest.raises(ValueError):
        nevpt2_rows(mol, mf, mo_full, NCAS, NELECAS, NROOTS,
                    variant="caspt2")


def test_sa_solver_objects_refuse_at_the_kernel(carbon):
    """The engine's SA guard: a state-averaged solver object never
    enters NEVPT (the multi-root readout is the required form)."""
    from pyscf import mcscf
    from pyscf.mrpt import NEVPT

    _, mol, mf, _ = carbon
    mc = mcscf.CASSCF(mf, NCAS, NELECAS).state_average_(
        weights=(0.5, 0.5))
    with pytest.raises(RuntimeError):
        NEVPT(mc, root=0).kernel()
