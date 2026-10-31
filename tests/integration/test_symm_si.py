"""The symmetry lane on the second p^2 atom: the runner anchor and
the SA head-to-head on Si/cc-pVDZ (the C pins live in the unit tier
on the frozen twin; Si converges a fresh reference once per
session)."""

import numpy as np
import pytest
from pyscf import gto, scf

from atomterm.casscf import run_irrep_casci, run_irrep_sa, run_sa_casscf
from atomterm.irreps import IRREP_IDS, det_irreps, orb_irreps
from atomterm.operators import (
    build_operators,
    ci_in_detbasis_order,
    window_closure,
)
from atomterm.system import AtomConfig
from atomterm.windows import build_full_mo, occupied_core_idx, select_shell_ao

pytestmark = pytest.mark.slow

_NCORE = 6  # Si: the Ne-like core; 3p^2 active


@pytest.fixture(scope="module")
def si_pair():
    """Both substrates of one Si reference: the C1 mol (route A) and
    its D2h twin (route B), each converged natively."""
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    mol1 = gto.M(atom=f"{cfg.symbol} 0 0 0", basis=cfg.basis, spin=cfg.spin,
                 verbose=0)
    mf1 = scf.ROHF(mol1)
    mf1.verbose = 0
    mf1.conv_tol = 1e-12
    mf1.kernel()
    assert mf1.converged
    mol2 = gto.M(atom=f"{cfg.symbol} 0 0 0", basis=cfg.basis, spin=cfg.spin,
                 symmetry="D2h", verbose=0)
    mf2 = scf.ROHF(mol2)
    mf2.verbose = 0
    mf2.conv_tol = 1e-12
    mf2.kernel()
    assert mf2.converged
    return mol1, mf1, mol2, mf2


@pytest.fixture(scope="module")
def si_lane(si_pair):
    """The 3p-shell twin problem: window, orbsyms, our matrices."""
    mol1, mf1, mol2, mf2 = si_pair
    c_act = np.eye(mol2.nao_nr())[:, select_shell_ao(mol2, "3p")]
    orbsyms = orb_irreps(mol2, c_act)
    ops = build_operators(mol2, mf2, c_act, _NCORE, 2, 0)
    return mol1, mf1, mol2, mf2, c_act, orbsyms, ops


def test_si_runner_anchor(si_lane):
    # the pinned component roots == eigenvalues of our H on the same
    # blocks, and OUR meters read the term exactly on THEIR vectors
    # (C's observed agreement: 4e-16 energies, 2e-15 purity)
    _, _, mol2, mf2, c_act, orbsyms, ops = si_lane
    roots = run_irrep_casci(mol2, mf2, c_act, _NCORE, 2, 0, 1, "g")
    ids = det_irreps(orbsyms, 2, 0)
    for irrep, (e, ci) in roots.items():
        sel = ids == IRREP_IDS[irrep]
        ours = np.linalg.eigvalsh(ops["H"][np.ix_(sel, sel)])
        assert np.allclose(ours, np.sort(e), atol=1e-8), irrep
        w = ci_in_detbasis_order(ci, 3, 2, 0)
        assert abs(float(np.real(w @ ops["L2"] @ w)) - 2.0) < 1e-9, irrep
        assert abs(float(np.real(w @ ops["S2"] @ w)) - 2.0) < 1e-9, irrep


def test_si_sa_head_to_head(si_lane, tmp_path):
    # the C finding holds on the second p^2 atom: identical ensembles
    # (observed 2e-12 on C), route B's window closed to the floor
    # (observed C: lz 1.0e-15 / lperp 3.1e-8) against route A's
    # optimizer-tolerance closure (observed C: 2.3e-6 / 2.8e-6)
    mol1, mf1, mol2, mf2, c_act, _, _ = si_lane

    route_b = run_irrep_sa(mol2, mf2, c_act, _NCORE, 2, 0, 1, "g")
    assert route_b["converged"]
    mo1 = build_full_mo(mf1, occupied_core_idx(mf1, _NCORE),
                        np.eye(mol1.nao_nr())[:, select_shell_ao(mol1, "3p")])
    route_a = run_sa_casscf(mf1, 3, (2, 0), mo1, [1 / 3] * 3,
                            str(tmp_path / "si_sa.chk"))

    assert abs(float(np.mean(route_a.energies)) - route_b["e_tot"]) < 1e-7
    # (observed Si |dEns| ~3e-11 over three fresh draws)

    active = slice(_NCORE, _NCORE + 3)
    cert_a = window_closure(mol1, route_a.mo_coeff[:, active])
    # observed Si: lz 1.5e-5, lperp 2.5e-5 (C ran 2.3e-6 / 2.8e-6 --
    # the free optimizer's closure tracks its tolerance, per atom);
    # one draw grazed 5e-5, hence the headroom
    assert cert_a["lz"] < 1e-4 and cert_a["lperp"] < 1e-4
    cert_b = window_closure(mol1, route_b["mo"][:, active])
    # observed Si: lz and lperp both 5.1e-8 (C ran 1.0e-15 -- its
    # barely-rotating SA stays at the raw-shell floor; Si's converged
    # rotations close at the optimizer level, still far under route A)
    assert cert_b["lz"] < 1e-6 and cert_b["lperp"] < 1e-6
