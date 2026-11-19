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

# the fixture's fresh root-0 E(2) (Eh) -- the absolute anchor against
# the engine itself, observed at pin time (2026-11-19, pyscf 2.13.1).
# NOTE the coincidence that is not one: the fixture's construction
# (frozen carbon, pur core, pure-AO 2p readout) is the campaign's
# poverty-class cell class -- c_2p_2-0_casci_ao_-_pur_avg prices the
# same E(2), -265.7585 mEh, from an independent substrate rebuild
# (agreement ~1e-7 Eh, the engine wobble floor)
E2_FRESH = -0.2657585148


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
    """The Sz=+1 p2 sector holds only the 3P states: one exact energy."""
    mol, mf, mc, mo_full = block
    e = np.atleast_1d(np.asarray(mc.e_cas, dtype=float))
    assert e.shape == (3,)
    assert np.max(e) - np.min(e) < 1e-9


def test_raw_core_field_splits_the_block(carbon):
    """The raw-core field splits the block (the unit-scale detector,
    fresh observed 4.4e-5 Eh -- the band below carries the record)."""
    _, mol, mf, act_idx = carbon
    from atomterm.windows import build_full_mo

    mo_raw = build_full_mo(mf, [0, 1], act_idx)
    mc = readout_mc(mol, mf, mo_raw, NCAS, NELECAS, NROOTS)
    e = np.atleast_1d(np.asarray(mc.e_cas, dtype=float))
    assert 1e-6 < np.max(e) - np.min(e) < 1e-3


def test_root_energies_follow_the_store_convention(block):
    """Root totals = e_cas + the common core part (store convention)."""
    mol, mf, mc, mo_full = block
    from atomterm.perturb import _root_energies

    e = _root_energies(mc)
    e_tot0 = float(np.atleast_1d(np.asarray(mc.e_tot)).ravel()[0])
    assert abs(e[0] - e_tot0) < 1e-11


def test_stored_roots_assertion_passes_on_truth_and_refuses_drift(block):
    """The assertion passes on truth and refuses energy/vector drift."""
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
    """A degenerate-block gauge redraw passes; a zapped root refuses."""
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
    """Rows carry the full class breakdown and match the kernel
    (fresh-readout pin -- the kernel writes the canonicalized CI into
    the shared ``mc.ci``)."""
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
    """Refs re-read and pass; a tampered store refuses before PT2."""
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
    """A state-averaged solver object refuses at the kernel (SA guard)."""
    from pyscf import mcscf
    from pyscf.mrpt import NEVPT

    _, mol, mf, _ = carbon
    mc = mcscf.CASSCF(mf, NCAS, NELECAS).state_average_(
        weights=(0.5, 0.5))
    with pytest.raises(RuntimeError):
        NEVPT(mc, root=0).kernel()


def test_absolute_e2_pin_on_the_fixture(block):
    """The fixture's root-0 E(2) sits at its observed fresh value."""
    mol, mf, _, mo_full = block
    rows = nevpt2_rows(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    assert rows[0].e_corr == pytest.approx(E2_FRESH, abs=1e-6)


def test_refusal_bands_are_forceable(block):
    """The forced bands pin their own semantics: e_tol the energy
    band; o_tol the same-state angle band (sub-percent noise passes,
    a one-percent mixture refuses); c_tol governs exact amplitudes
    with o_tol forced off; a length-mismatched energy stack refuses
    by shape."""
    mol, mf, _, mo_full = block
    rows = nevpt2_rows(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    from atomterm.perturb import _root_energies, readout_mc

    mc = readout_mc(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    e = _root_energies(mc).copy()
    # (a) energy noise inside the default band passes; the tightened
    # band refuses the same stack
    assert_stored_roots(mc, ref_energies=e + 5e-13)
    with pytest.raises(ValueError):
        assert_stored_roots(mc, ref_energies=e + 5e-13, e_tol=1e-15)
    # (b) the o_tol band: a one-percent amplitude mixing is a
    # different state at the default ...
    mc1 = readout_mc(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    ci_mix = [np.asarray(c).ravel() + 0.01 * np.eye(c.size)[1]
              for c in mc1.ci]
    with pytest.raises(ValueError):
        assert_stored_roots(mc1, ref_ci=ci_mix)
    # ... and the widened o_tol admits it (forced-pin)
    mc1b = readout_mc(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    assert_stored_roots(mc1b, ref_ci=ci_mix, o_tol=0.05)
    # (c) the c_tol band governs exact amplitudes with o_tol off:
    # 5e-9 structured noise passes at the default ...
    mc2 = readout_mc(mol, mf, mo_full, NCAS, NELECAS, NROOTS)
    ci_pert = [np.asarray(c).ravel() + 5e-9 * np.eye(c.size)[1]
               for c in mc2.ci]
    assert_stored_roots(mc2, ref_ci=ci_pert, o_tol=0.0, c_tol=1e-8)
    # ... and the tightened c_tol refuses it
    with pytest.raises(ValueError):
        assert_stored_roots(mc2, ref_ci=ci_pert, o_tol=0.0, c_tol=1e-12)
    # (d) a wrong-length energy stack refuses by shape (the declared
    # message, not a numpy broadcast error)
    with pytest.raises(ValueError, match="shape mismatch"):
        assert_stored_roots(mc, ref_energies=np.append(e, 0.0))


def test_phase_and_davidson_noise_clauses(block):
    """The per-root clauses: a global-phase flip and a sub-percent
    davidson-noise mixing of one root pass; a physical mixture (ten
    percent) refuses."""
    mol, mf, _, mo_full = block
    from atomterm.perturb import readout_mc

    mc = readout_mc(mol, mf, mo_full, NCAS, NELECAS, 1)
    v = np.asarray(mc.ci).ravel()
    # a phase flip passes (single-root form: the array IS the root)
    mc2 = readout_mc(mol, mf, mo_full, NCAS, NELECAS, 1)
    mc2.ci = -v.reshape(np.asarray(mc.ci).shape)
    assert_stored_roots(mc2, ref_ci=[v])
    # a 1e-4 davidson-noise mixing passes (the declared class)
    partner = np.zeros_like(v)
    partner[1] = 1.0
    v_noise = v + 1e-4 * partner
    mc3 = readout_mc(mol, mf, mo_full, NCAS, NELECAS, 1)
    mc3.ci = v_noise.reshape(np.asarray(mc.ci).shape)
    assert_stored_roots(mc3, ref_ci=[v])
    # a ten-percent mixing is a different state: refuses
    v_mix = v + 0.1 * partner
    mc4 = readout_mc(mol, mf, mo_full, NCAS, NELECAS, 1)
    mc4.ci = v_mix.reshape(np.asarray(mc.ci).shape)
    with pytest.raises(ValueError):
        assert_stored_roots(mc4, ref_ci=[v])


def test_class_rows_carry_the_compression_meter():
    """e_per_norm = energy/norm; a zero norm carries no meter (nan)."""
    from atomterm.perturb import ClassRow

    row = ClassRow("Sir", "(0)'", norm=2.0, energy=-4.0)
    assert row.e_per_norm == pytest.approx(-2.0)
    assert ClassRow("Sr", "(-1)'", norm=0.0, energy=0.0).e_per_norm != \
        ClassRow("Sr", "(-1)'", norm=0.0, energy=0.0).e_per_norm  # nan
