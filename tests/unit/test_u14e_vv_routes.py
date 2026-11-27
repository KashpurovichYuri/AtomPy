"""The routes of the V&V library (scripts/u14e/vv/routes.py) on the
frozen carbon reference's 2p^2 window: every solution is a certified
object (its orbitals orthonormal, its CI stationary in ITS OWN
orbitals, in its solvers' core fold), the inconsistent pair that
mispriced six table rows cannot pass the certificate, the free and
the pinned routes agree on the frozen window and on the state average,
the two pricing adapters agree, the Hamiltonian of the certificate
matches the determinant-basis construction, and the row conventions
follow the route (no SCF; tests/carbon_reference.py)."""

import numpy as np
import pytest

from atomterm.operators import build_operators
from atomterm.operators import ci_in_detbasis_order
from atomterm.windows import lowdin

from vv import cases as ladder
from vv.pricing import price_classes, price_pyscf
from vv.routes import (State, assert_consistent, energy_expectation,
                       free_casci, free_casscf, free_sa, pinned_casci,
                       pinned_sa)
from vv.table import Route
from vv.windows import window_from_shells

pytestmark = pytest.mark.fast

TERM = (1, "g")        # the 3P term of the 2p^2 window
SECTOR = (2, 0)
SHELLS = (("2p", 3),)


@pytest.fixture(scope="module")
def space(vv_substrate):
    """The 2p window on the purified core."""
    return ladder.visible_window(vv_substrate, SHELLS, SECTOR).space


@pytest.fixture(scope="module")
def frozen(space, vv_substrate):
    """(the free frozen solution priced, the pinned one priced)."""
    return (free_casci(vv_substrate, space).price(),
            pinned_casci(vv_substrate, space, TERM).price())


# ------------------------------------------------------ the birth certificate

def test_every_state_is_stationary_in_its_own_orbitals(frozen, space):
    for priced in frozen:
        cert = priced.solution.certificate
        assert cert["ortho"] < 1e-13
        assert cert["residual"] < 1e-12
        # re-certified by hand, in each route's own core fold
        fold = "raw" if priced.solution.route is Route.FREE else "averaged"
        again = assert_consistent([p.state for p in priced.states], space, fold)
        assert again["residual"] < 1e-12


def test_energy_expectation_matches_the_determinant_basis_hamiltonian(
        frozen, space, vv_substrate):
    # an independent construction of the same Hamiltonian: the
    # determinant-basis H (Slater-Condon over explicit determinants)
    # through the CI-order bridge, against the FCI feed of the
    # certificate -- both folding the SAME core columns (left to its
    # default, build_operators folds the canonical core, a different
    # Hamiltonian by 6e-07 Eh here)
    st = frozen[0].states[0].state
    parts = build_operators(vv_substrate.mol, vv_substrate.mf, st.active,
                            st.ncore, *SECTOR, c_core=st.core)
    w = ci_in_detbasis_order(st.ci.ravel(), st.ncas, *SECTOR)
    e_det = float(np.real(w @ parts["H"] @ w))
    e_fci = energy_expectation(st.mol, st.mf, st.active, st.core, SECTOR,
                               st.ci, fold="averaged")
    assert abs(e_det - e_fci) < 1e-10
    assert abs(e_fci - st.e_cas) < 1e-12


def test_each_route_is_certified_in_its_solvers_core_fold(vv_substrate):
    # the CANONICAL carbon core is polarized by the open shell (s/d_z^2
    # mixing): PySCF's fold (raw) and the pinned solvers' fold (the
    # spherical average) then differ by the polarization energy
    # (observed 2.9e-05 Eh) -- each route passes in its own fold only
    core = np.asarray(vv_substrate.mf.mo_coeff)[:, :2]
    raw = window_from_shells(vv_substrate.mol, vv_substrate.mf, SHELLS,
                             core, SECTOR).space
    free = free_casci(vv_substrate, raw)
    pinned = pinned_casci(vv_substrate, raw, TERM)
    for sol, own, other in ((free, "raw", "averaged"),
                            (pinned, "averaged", "raw")):
        states = [s for s in sol.states]
        assert assert_consistent(states, raw, own)["residual"] < 1e-12
        with pytest.raises(ValueError, match="not stationary"):
            assert_consistent(states, raw, other)


def test_the_start_window_pairing_is_refused(vv_substrate, space):
    # the defect's shape: the CI vectors of the CONVERGED window placed
    # on the START window's columns.  Clause (1) catches it as written
    # (the set is not orthonormal); clause (2) catches it even when
    # the foreign columns are re-orthonormalized
    sol = pinned_sa(vv_substrate, space, TERM)
    converged = sol.states[0].mo_coeff
    start = space.active
    # the optimizer moved the window (otherwise the test is vacuous)
    cos = np.linalg.svd(start.T @ vv_substrate.mf.get_ovlp()
                        @ sol.states[0].active, compute_uv=False)
    assert cos.min() < 0.99

    def franken(mo):
        out = []
        for st in sol.states:
            mc = object.__new__(type(st.mc))
            mc.__dict__.update(st.mc.__dict__)
            mc.mo_coeff = mo
            out.append(State(mc, st.root, st.irrep))
        return out

    swapped = np.array(converged)
    swapped[:, 2:5] = start
    with pytest.raises(ValueError, match="not S-orthonormal"):
        assert_consistent(franken(swapped), space, "averaged")
    s = vv_substrate.mf.get_ovlp()
    orthonormal_set = lowdin(swapped, s)
    with pytest.raises(ValueError, match="not stationary in its own orbitals"):
        assert_consistent(franken(orthonormal_set), space, "averaged")
    assert assert_consistent(franken(converged), space, "averaged")["residual"] < 1e-12


def test_the_partition_must_be_the_active_spaces(space, vv_substrate):
    sol = free_casci(vv_substrate, space)
    # the same shells, another sector: a legitimate space, not the one
    # the solution was solved on
    core = space.core
    other = window_from_shells(vv_substrate.mol, vv_substrate.mf, SHELLS,
                               core, (1, 1)).space
    with pytest.raises(ValueError, match="partition"):
        assert_consistent(list(sol.states), other)


# ------------------------------------------------ the routes against each other

def test_the_free_and_the_pinned_routes_agree_on_the_frozen_window(frozen):
    free, pinned = frozen
    assert len(pinned.states) == 3                 # the three 3P components
    assert abs(free.states[0].state.e_ref - pinned.states[0].state.e_ref) < 1e-10
    assert np.ptp([p.state.e_ref for p in pinned.states]) < 1e-10
    # the same E(2) through the two pricing adapters, too
    assert abs(free.states[0].e2 - np.mean(pinned.e2)) < 1e-9
    assert all(abs(p.l2 - 2.0) < 1e-9 and abs(p.s2 - 2.0) < 1e-9
               for p in list(free.states) + list(pinned.states))


def test_the_state_averages_land_the_same_ensemble(vv_substrate, space):
    free = free_sa(vv_substrate, space, nroots=3)
    pinned = pinned_sa(vv_substrate, space, TERM)
    # PySCF's state average IS the mean of its root energies
    assert abs(free.ensemble - np.mean([s.e_ref for s in free.states])) < 1e-10
    assert abs(free.ensemble - pinned.ensemble) < 1e-8
    # the variational ordering: the single-root optimum lies at or
    # below the ensemble (it may break the term to get there)
    assert free_casscf(vv_substrate, space).states[0].e_ref <= free.ensemble + 1e-12


def test_the_single_root_block_carries_its_bare_vector(vv_substrate, space):
    # mrpt's load_ci returns mc.ci as it is when nroots == 1: a
    # single-root block must hold the (nsa, nsb) vector, not a list
    sol = pinned_casci(vv_substrate, space, TERM)
    for st in sol.states:
        assert np.shape(st.ci) == (3, 1)       # C(3,2) alpha strings, 1 beta
        assert price_classes(st).e_corr < 0


def test_the_pricing_adapters_agree(vv_substrate, space):
    sol = free_sa(vv_substrate, space, nroots=3)
    for st in sol.states:
        assert abs(price_pyscf(st) - price_classes(st).e_corr) < 1e-9
    by_kernel, by_classes = sol.price(), sol.price(how="classes")
    assert by_kernel.solution.route is Route.FREE and by_kernel.states[0].classes is None
    assert len(by_classes.states[0].classes) == 8      # the Koopmans classes
    assert by_classes.states[0].classes[0].label == "Sr"


def test_the_optimizer_leaks_its_core_and_it_is_measured(vv_substrate, space):
    # a space used as built does not leak; the optimizer's returned
    # core block sits off the frozen span at its convergence noise
    assert free_casci(vv_substrate, space).certificate["core_leak"] < 1e-13
    pinned = pinned_sa(vv_substrate, space, TERM)
    assert 0.0 <= pinned.certificate["core_leak"] < 1e-4
    assert 0.0 <= pinned.certificate["readout_shift"] < 1e-5
    # the readout describes the optimizer's solution: its component
    # mean is the ensemble (each state's core bookkeeping is the core
    # its readout folded -- a purified-core constant sits 1e-06 off)
    assert pinned.certificate["ensemble_dev"] < 1e-8


# ------------------------------------------------------------ the row conventions

def test_a_free_row_carries_its_first_root_and_a_pinned_row_the_component_mean(
        vv_substrate, space):
    free = free_sa(vv_substrate, space, nroots=3).price().row("(x)")
    pinned = pinned_sa(vv_substrate, space, TERM).price().row("(x)")
    assert (free.method, free.route) == ("sa-casscf", Route.FREE)
    assert (pinned.method, pinned.route) == ("sa-casscf", Route.PINNED)
    assert free.key == "sa-casscf|free@pur|(x)"
    assert pinned.key == "sa-casscf|pinned@pur|(x)"
    sol = free_sa(vv_substrate, space, nroots=3)
    priced = sol.price()
    assert free.e_method == pytest.approx(priced.states[0].state.e_ref)
    assert free.e2 == pytest.approx(priced.states[0].e2)
    assert free.e_total == pytest.approx(free.e_method + free.e2)
    p = pinned_sa(vv_substrate, space, TERM).price()
    assert pinned.e_method == pytest.approx(p.solution.ensemble)
    assert pinned.e2 == pytest.approx(np.mean(p.e2))
    assert pinned.e_total == pytest.approx(pinned.e_method + pinned.e2)


def test_occupations_and_the_empty_active_orbitals(vv_substrate, space):
    # the frozen 3P ground of the 2p^2 window (px^1 py^1, Ms = 1) has
    # its p_z empty: one active orbital below PySCF's threshold
    from vv.channels import Channel, channel_block, virtual_spectrum
    from vv.pricing import (FRAC_OCC_THRESHOLD, occupations, partition_as_solved,
                            resorted)
    from vv.windows import window_with_channels

    st = free_casci(vv_substrate, space).states[0]
    occ = occupations(st)
    assert occ.shape == (3,) and np.all(np.diff(occ) <= 1e-12)
    assert occ.sum() == pytest.approx(2.0)                     # the active electrons
    assert occ[:2] == pytest.approx([1.0, 1.0]) and occ[2] < FRAC_OCC_THRESHOLD
    assert resorted(st) == 1                                  # observed: p_z at 0.0


def test_an_empty_active_orbital_is_repartitioned_by_the_prelude(vv_substrate,
                                                                space):
    # THE PARTITION hazard, measured: a non-coupling s orbital beside
    # the window leaves the frozen energy and the term unchanged, yet
    # PySCF's default prelude moves E(2) by ~100 mEh (a second empty
    # slot, filled by the next-lowest virtual); with the window priced
    # as solved the two readings differ by < 1 mEh (observed 0.57)
    from vv.channels import Channel, channel_block, virtual_spectrum
    from vv.pricing import partition_as_solved, resorted
    from vv.windows import window_with_channels

    j = next(r["index"] for r in virtual_spectrum(vv_substrate) if r["l"] == "s")
    built = ladder.visible_window(vv_substrate, (("2p", 3),), SECTOR)
    extended = window_with_channels(
        vv_substrate, built, [channel_block(vv_substrate, Channel("4s", (j,), ("3s",)))])
    a = free_casci(vv_substrate, space).states[0]
    b = free_casci(vv_substrate, extended.space).states[0]
    assert abs(b.e_ref - a.e_ref) < 1e-10                 # a spectator
    assert resorted(a) == 1 and resorted(b) == 2
    default = price_pyscf(b) - price_pyscf(a)
    assert default > 0.05                                  # observed +0.1007 Eh
    with partition_as_solved():
        solved = price_pyscf(b) - price_pyscf(a)
        assert abs(price_pyscf(a) - price_classes(a).e_corr) < 1e-9
    assert abs(solved) < 1e-3                              # observed -0.00057 Eh
    # the context restores PySCF's convention
    assert resorted(b) == 2


def test_meters_without_pricing(vv_substrate, space):
    sol = free_sa(vv_substrate, space, nroots=3)
    meters = sol.meters()
    assert len(meters) == 3
    assert all(abs(l2 - 2.0) < 1e-9 and abs(s2 - 2.0) < 1e-9 for l2, s2 in meters)
    assert len(sol.meters(roots=[1])) == 1
