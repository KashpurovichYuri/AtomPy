"""The route A/route B anchor: PySCF's wfnsym-pinned CI against our
exact-matrix irrep blocks, on identical folded integrals (the frozen
D2h carbon twin), plus the symmetry-lane runners on the same
substrate."""

import numpy as np
import pytest
from pyscf import ao2mo, fci
from pyscf.lib.exceptions import WfnSymmetryError

from atomterm.casscf import run_irrep_casci, run_irrep_sa, run_sa_casscf
from atomterm.irreps import (
    IRREP_IDS,
    IRREP_NAMES,
    det_irreps,
    orb_irreps,
)
from atomterm.operators import (
    build_operators,
    ci_in_detbasis_order,
    window_closure,
)
from atomterm.windows import build_full_mo, occupied_core_idx

pytestmark = pytest.mark.fast

_SECTORS = ((2, 0), (1, 1))  # the triplet and the alpha-beta p^2 blocks


@pytest.fixture(scope="module")
def twin(carbon_d2h):
    """The twin problem over the 2p AO shell of the frozen D2h
    reference: our folded h1 (the m-scalar core fold), the window
    eri in PySCF's own feed, and our determinant-basis matrices."""
    _, mol, mf, act_idx = carbon_d2h
    c_act = np.eye(mol.nao_nr())[:, act_idx]
    orbsyms = orb_irreps(mol, c_act)
    eri = ao2mo.kernel(mol, c_act)
    h1 = None
    mats = {}
    for na, nb in _SECTORS:
        parts = build_operators(mol, mf, c_act, 2, na, nb)
        mats[(na, nb)] = {k: parts[k] for k in ("H", "L2", "S2")}
        h1 = parts["h1"]  # sector-independent: the fold sees only the core
    return mol, mf, c_act, h1, eri, orbsyms, mats


def _block_solver(twin, na, nb, irrep, nroots):
    """Route B: the symmetry-pinned CI on the shared integrals."""
    mol, _, _, h1, eri, orbsyms, _ = twin
    solver = fci.direct_spin1_symm.FCI(mol)
    solver.orbsym = orbsyms
    solver.wfnsym = irrep
    solver.nroots = nroots
    solver.conv_tol = 1e-12
    e, ci = solver.kernel(h1, eri, 3, (na, nb))
    return np.atleast_1d(np.asarray(e, dtype=float)), np.asarray(ci)


def test_block_energies_match_our_matrices(twin):
    # per irrep: their pinned roots == eigenvalues of our H on the
    # same block -- same integrals, two independent solutions
    # (observed max |dE| 4e-16)
    mol, _, _, _, _, orbsyms, mats = twin
    for (na, nb), m in mats.items():
        ids = det_irreps(orbsyms, na, nb)
        for ir in np.unique(ids):
            sel = ids == ir
            ours = np.linalg.eigvalsh(m["H"][np.ix_(sel, sel)])
            theirs, _ = _block_solver(twin, na, nb, IRREP_NAMES[ir],
                                      int(sel.sum()))
            assert np.allclose(ours, np.sort(theirs), atol=1e-8), (na, nb, ir)


def test_sector_spectrum_is_block_union(twin):
    # the partition is a BLOCK structure: our full-sector spectrum ==
    # the union of the per-block spectra, so every off-diagonal
    # inter-irrep element of H vanishes (observed 1e-15)
    _, _, _, _, _, orbsyms, mats = twin
    for (na, nb), m in mats.items():
        ids = det_irreps(orbsyms, na, nb)
        blocks = np.concatenate([
            np.linalg.eigvalsh(m["H"][np.ix_(ids == ir, ids == ir)])
            for ir in np.unique(ids)
        ])
        assert np.allclose(np.linalg.eigvalsh(m["H"]), np.sort(blocks),
                           atol=1e-10), (na, nb)


def test_ground_states_term_pure(twin):
    # the triplet block ground states read L=1, S=1 by OUR meters on
    # THEIR vectors (observed deviation 2e-15)
    _, _, _, _, _, orbsyms, mats = twin
    m = mats[(2, 0)]
    ids = det_irreps(orbsyms, 2, 0)
    for ir in np.unique(ids):
        _, ci = _block_solver(twin, 2, 0, IRREP_NAMES[ir], 1)
        w = ci_in_detbasis_order(ci, 3, 2, 0)
        assert abs(float(np.real(w @ m["L2"] @ w)) - 2.0) < 1e-9, ir
        assert abs(float(np.real(w @ m["S2"] @ w)) - 2.0) < 1e-9, ir


def test_component_degeneracy(twin):
    # on the m-scalar fold the three L=1 components are degenerate
    # (observed spread 0.0 -- exact to the last bit)
    _, _, _, _, _, orbsyms, _ = twin
    ids = det_irreps(orbsyms, 2, 0)
    e = [_block_solver(twin, 2, 0, IRREP_NAMES[ir], 1)[0][0]
         for ir in np.unique(ids)]
    assert np.ptp(e) < 1e-9


def test_reachability_contract(twin):
    # the blocks partition the sector; irreps with no determinant
    # raise at the solver (the parity gate)
    _, _, _, _, _, orbsyms, mats = twin
    for (na, nb), m in mats.items():
        ids = det_irreps(orbsyms, na, nb)
        sizes = {int(ir): int((ids == ir).sum()) for ir in np.unique(ids)}
        assert sum(sizes.values()) == m["H"].shape[0]
        for ir in set(range(8)) - set(sizes):
            with pytest.raises(WfnSymmetryError):
                _block_solver(twin, na, nb, IRREP_NAMES[ir], 1)


def test_runner_matches_anchor(twin):
    # run_irrep_casci(L=1, g) reproduces the anchor: per-block pinned
    # roots == our restricted-H eigenvalues, meters exact (observed
    # 4e-16 / 2e-15)
    mol, mf, c_act, _, _, orbsyms, mats = twin
    m = mats[(2, 0)]
    ids = det_irreps(orbsyms, 2, 0)
    roots = run_irrep_casci(mol, mf, c_act, 2, 2, 0, 1, "g")
    assert set(roots) == {IRREP_NAMES[ir] for ir in np.unique(ids)}
    for irrep, (e, ci) in roots.items():
        sel = ids == IRREP_IDS[irrep]
        ours = np.linalg.eigvalsh(m["H"][np.ix_(sel, sel)])
        assert np.allclose(ours, np.sort(e), atol=1e-8), irrep
        w = ci_in_detbasis_order(ci, 3, 2, 0)
        assert abs(float(np.real(w @ m["L2"] @ w)) - 2.0) < 1e-9, irrep


def test_sa_head_to_head(twin, carbon, tmp_path):
    # The two SA philosophies on the same physical problem, each in
    # its native habitat.  Route A (symmetry-blind) runs the
    # production driver on the C1 substrate: PySCF's SA on a
    # symmetry-adapted molecule is itself wfnsym-pinned (its wrapper
    # forces one irrep per solver), so the unpinned sector average
    # exists only on C1.  The variational problem is identical --
    # same basis, same raw-shell start window.
    mol, mf, c_act, *_ = twin
    _, mol1, mf1, act1 = carbon

    route_b = run_irrep_sa(mol, mf, c_act, 2, 2, 0, 1, "g")
    assert route_b["converged"]
    mo1 = build_full_mo(mf1, occupied_core_idx(mf1, 2),
                        np.eye(mol1.nao_nr())[:, act1])
    route_a = run_sa_casscf(mf1, 3, (2, 0), mo1, [1 / 3] * 3,
                            str(tmp_path / "sa_a.chk"))

    ens_a = float(np.mean(route_a.energies))  # uniform over the set
    assert abs(ens_a - route_b["e_tot"]) < 1e-8  # observed 2e-12

    # the finding: route B's irrep-blocked rotations close the window
    # to the floor; route A's free C1 optimizer closes only to its
    # convergence tolerance (observed A: lz 2.3e-6, lperp 2.8e-6;
    # B: lz 1.0e-15, lperp 3.1e-8)
    cert_a = window_closure(mol1, route_a.mo_coeff[:, 2:5])
    assert cert_a["lz"] < 1e-5 and cert_a["lperp"] < 1e-5
    cert_b = window_closure(mol1, route_b["mo"][:, 2:5])
    assert cert_b["lz"] < 1e-10 and cert_b["lperp"] < 1e-6

    ops = build_operators(mol, mf, route_b["mo"][:, 2:5], 2, 2, 0)
    for _, (_, ci) in route_b["roots"].items():
        w = ci_in_detbasis_order(ci, 3, 2, 0)
        assert abs(float(np.real(w @ ops["L2"] @ w)) - 2.0) < 1e-8
