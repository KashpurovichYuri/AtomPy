"""The symmetry lane on the second p^2 atom: the runner anchor and
the SA head-to-head on Si/cc-pVDZ (the C pins live in the unit tier
on the frozen twin; Si converges a fresh reference once per
session)."""

import numpy as np
import pytest
from pyscf import gto, scf

from atomterm.casscf import (constrained_casscf, run_irrep_casci,
                             run_irrep_sa)
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
    """Both substrates of one Si reference: the C1 mol (the free
    route) and its D2h twin (the pinned route), each converged
    natively."""
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    mol_c1 = gto.M(atom=f"{cfg.symbol} 0 0 0", basis=cfg.basis, spin=cfg.spin,
                 verbose=0)
    mf_c1 = scf.ROHF(mol_c1)
    mf_c1.verbose = 0
    mf_c1.conv_tol = 1e-12
    mf_c1.kernel()
    assert mf_c1.converged
    mol_d2h = gto.M(atom=f"{cfg.symbol} 0 0 0", basis=cfg.basis, spin=cfg.spin,
                 symmetry="D2h", verbose=0)
    mf_d2h = scf.ROHF(mol_d2h)
    mf_d2h.verbose = 0
    mf_d2h.conv_tol = 1e-12
    mf_d2h.kernel()
    assert mf_d2h.converged
    return mol_c1, mf_c1, mol_d2h, mf_d2h


@pytest.fixture(scope="module")
def si_lane(si_pair):
    """The 3p-shell twin problem: window, orbsyms, our matrices."""
    mol_c1, mf_c1, mol_d2h, mf_d2h = si_pair
    c_act = np.eye(mol_d2h.nao_nr())[:, select_shell_ao(mol_d2h, "3p")]
    orbsyms = orb_irreps(mol_d2h, c_act)
    ops = build_operators(mol_d2h, mf_d2h, c_act, _NCORE, 2, 0)
    return mol_c1, mf_c1, mol_d2h, mf_d2h, c_act, orbsyms, ops


def test_si_runner_anchor(si_lane):
    # the pinned component roots == eigenvalues of our H on the same
    # blocks, and OUR meters read the term exactly on THEIR vectors
    # (C's observed agreement: 4e-16 energies, 2e-15 purity)
    _, _, mol_d2h, mf_d2h, c_act, orbsyms, ops = si_lane
    roots = run_irrep_casci(mol_d2h, mf_d2h, c_act, _NCORE, 2, 0, 1, "g")
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
    # (observed 2e-12 on C), the pinned route's window closed to the floor
    # (observed C: lz 1.0e-15 / lperp 3.1e-8) against the free route's
    # optimizer-tolerance closure (observed C: 2.3e-6 / 2.8e-6)
    mol_c1, mf_c1, mol_d2h, mf_d2h, c_act, _, _ = si_lane

    pinned_sa = run_irrep_sa(mol_d2h, mf_d2h, c_act, _NCORE, 2, 0, 1, "g")
    assert pinned_sa["converged"]
    mo_start_c1 = build_full_mo(mf_c1, occupied_core_idx(mf_c1, _NCORE),
                        np.eye(mol_c1.nao_nr())[:, select_shell_ao(mol_c1, "3p")])
    free_sa = constrained_casscf(mf_c1, 3, (2, 0),
                                 str(tmp_path / "si_sa.chk"),
                                 mo_full=mo_start_c1,
                                 weights=[1 / 3] * 3, name="sa-casscf")

    assert abs(float(np.mean(free_sa.energies)) - pinned_sa["e_tot"]) < 1e-7
    # (observed Si |dEns| ~3e-11 over three fresh draws)

    active = slice(_NCORE, _NCORE + 3)
    cert_free = window_closure(mol_c1, free_sa.mo_coeff[:, active])
    # observed Si: lz 1.5e-5, lperp 2.5e-5 (C ran 2.3e-6 / 2.8e-6 --
    # the free optimizer's closure tracks its tolerance, per atom);
    # one draw grazed 5e-5, hence the headroom
    assert cert_free["lz"] < 1e-4 and cert_free["lperp"] < 1e-4
    cert_pinned = window_closure(mol_c1, pinned_sa["mo"][:, active])
    # observed Si: lz and lperp both 5.1e-8 (C ran 1.0e-15 -- its
    # barely-rotating SA stays at the raw-shell floor; Si's converged
    # rotations close at the optimizer level, still far under the free route)
    assert cert_pinned["lz"] < 1e-6 and cert_pinned["lperp"] < 1e-6
