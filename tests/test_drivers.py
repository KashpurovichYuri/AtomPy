"""Tests of atomterm.drivers: orbital selection, block-safe window
construction, S-metric weights, and the fractional-occupation reference.

The build_full_mo tests are the regression guard for the ticket-09 bug:
a naive complement projector silently returned a NON-orthonormal MO matrix
whenever the raw window AOs overlapped the canonical core (Si: ~0.3) --
and everything fed that matrix (CASCI tolerated it; the CASSCF optimizer
diverged to sub-variational energies).
"""

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from atomterm.system import AtomConfig
from atomterm.drivers import (
    build_full_mo,
    run_frac_occ_scf,
    run_rohf,
    select_shell_ao,
    select_valence_ao,
    shell_weights,
)


def _s_ortho_err(mo, s):
    """max |C^dag S C - I| of an orbital block (S-orthonormality error)."""
    mo = np.asarray(mo)
    return float(np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max())


def test_select_valence_ao_concatenates_in_pattern_order(carbon):
    _, mol, _, act_idx = carbon
    assert select_valence_ao(mol, ["2p"]) == act_idx
    assert select_valence_ao(mol, ["2s", "2p"]) == \
        select_shell_ao(mol, "2s") + select_shell_ao(mol, "2p")


def test_select_valence_ao_rejects_unknown_pattern(carbon):
    _, mol, _, _ = carbon
    with pytest.raises(ValueError, match="9f"):
        select_valence_ao(mol, ["9f"])


def test_shell_weights_partition_unity_and_group_by_l(carbon):
    _, mol, mf, act_idx = carbon
    s = mf.get_ovlp()
    c = np.eye(mol.nao_nr())[:, act_idx]  # raw 2p AOs: exactly l = 1
    names, w = shell_weights(mol, c, s)
    assert names == ["s", "p", "d"]
    assert np.allclose(w.sum(axis=1), 1.0, atol=1e-12)
    assert np.allclose(w[:, names.index("p")], 1.0, atol=1e-12)


@pytest.fixture(scope="session")
def build_full_mo_refs(tmp_path_factory):
    """ROHF per atom, once per session (chkfiles under pytest tmp, never
    leaked across sessions -- stale caches would mask physics changes)."""
    out = []
    for cfg, patterns in [
        (AtomConfig(), ["2p"]),
        (AtomConfig("Si", "ccpvdz", spin=2), ["3p"]),
    ]:
        chk = str(tmp_path_factory.mktemp(f"bfm_{cfg.symbol}") / "rohf.chk")
        mol, mf = run_rohf(cfg, chk)
        out.append((cfg, mol, mf, patterns))
    return out


def test_build_full_mo_blocks_orthonormal_and_span_preserving(build_full_mo_refs):
    """The ticket-09 regression guard, on both atoms.

    Si is the case that broke: raw 3p AOs overlap the canonical 2p core
    MOs by ~0.3, so the block-safe construction has real work to do.
    """
    for cfg, mol, mf, patterns in build_full_mo_refs:
        s = mf.get_ovlp()
        act_idx = select_valence_ao(mol, patterns)
        ncore = (mol.nelectron - 2) // 2

        mo = build_full_mo(mf, list(range(ncore)), act_idx)
        # (nao, nao) and fully S-orthonormal -- the old bug returned
        # (nao, nao+6) with err ~ 1 on Si
        assert mo.shape == (mol.nao_nr(), mol.nao_nr()), cfg.symbol
        assert _s_ortho_err(mo, s) < 1e-10, cfg.symbol

        # core SPAN preserved: the block projectors agree with the plain
        # canonical core projector (Löwdin within the span must not move it)
        c_can = mf.mo_coeff[:, :ncore]
        c_new = mo[:, :ncore]
        p_can = c_can @ c_can.conj().T @ s
        p_new = c_new @ c_new.conj().T @ s
        assert np.abs(p_new - p_can).max() < 1e-8, cfg.symbol

        # window l-purity survives the core-overlap correction
        c_win = mo[:, ncore : ncore + len(act_idx)]
        names, w = shell_weights(mol, c_win, s)
        assert w[:, names.index("p")].min() > 0.99, cfg.symbol


@pytest.fixture(scope="session")
def si_frac_chk(tmp_path_factory):
    return tmp_path_factory.mktemp("frac_si") / "frac.chk"


def test_frac_occ_scf_open_block_ordering(si_frac_chk):
    """Average-of-configuration reference for Si: the fractional block must
    land on the 3p shell -- closed [Ne]+3s below (the ordering assumption
    of run_frac_occ_scf, first nontrivial atom after carbon).

    Runs twice deliberately: first call = fresh SCF, second = chkfile cache
    hit. Both returned objects must carry the ensemble occupation hook --
    re-deriving occupations from energies must reproduce the stored mo_occ
    (the ticket-10 finding: the hook used to be missing on cache hits).
    """
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    Path(si_frac_chk).unlink(missing_ok=True)
    for k in range(2):  # 1st pass: fresh; 2nd pass: cache hit
        mol_f, mf_f = run_frac_occ_scf(
            dataclasses.replace(cfg, spin=0), str(si_frac_chk),
            nelec_open=2, n_open=3)
        assert mf_f.converged, k
        occ = np.asarray(mf_f.mo_occ)  # the occupations the SCF converged with
        n_closed = (mol_f.nelectron - 2) // 2
        assert np.allclose(occ[:n_closed], 2.0), k
        assert np.allclose(occ[n_closed : n_closed + 3], 2.0 / 3.0), k
        assert np.allclose(occ[n_closed + 3 :], 0.0), k
        # the hook is part of the object's identity on BOTH paths
        assert np.allclose(mf_f.get_occ(mf_f.mo_energy), occ), k
        # the fractional block sits on the 3p AOs (l-purity of the open shell)
        s = mf_f.get_ovlp()
        names, w = shell_weights(mol_f, mf_f.mo_coeff[:, n_closed : n_closed + 3], s)
        assert w[:, names.index("p")].min() > 0.99, k


def test_sc_rohf_reference_is_3d1_4s2_with_argon_core(sc_rohf):
    """Ticket 11's TM reference risks, pinned: the ROHF converges to the
    3d1 4s2 occupation with the canonical [:9] core = [Ar] exactly -- no
    occupation flipping, no 3d charge in the folded core."""
    cfg, mol, mf, act_idx = sc_rohf
    assert mf.converged
    occ = np.asarray(mf.mo_occ)
    ncore = 9
    assert np.allclose(occ[:10], 2.0) and occ[10] == 1.0
    # the 9 canonical core orbitals are the 9 lowest, and carry no d charge
    assert mf.mo_energy[:ncore].max() < mf.mo_energy[ncore]
    s = mf.get_ovlp()
    names, w = shell_weights(mol, mf.mo_coeff[:, :ncore], s)
    assert w[:, names.index("d")].max() < 1e-4  # measured ~5e-5, numerical


def test_sc_frac_occ_selects_open_shell_by_projection(sc_rohf, tmp_path):
    """The fact-15(c) invariant: on a TM the fractional block MUST be
    selected by character. Energy-ordered filling converges -- to a
    p-containing ensemble (4p intrudes into the 3d/4s block); the
    open_idx projection picks {3d, 4s} exactly."""
    cfg, mol, mf, act_idx = sc_rohf
    n_open = len(act_idx)

    chk = tmp_path / "frac_proj.chk"
    mol_p, mf_p = run_frac_occ_scf(
        cfg, str(chk), nelec_open=3, n_open=n_open, open_idx=act_idx)
    assert mf_p.converged
    frac_cols = np.where(np.isclose(np.asarray(mf_p.mo_occ), 0.5))[0]
    s = mf_p.get_ovlp()
    names, w = shell_weights(mol_p, mf_p.mo_coeff[:, frac_cols], s)
    wd = dict(zip(names, w.mean(axis=0)))
    assert abs(wd.get("d", 0.0) - 5 / 6) < 1e-3
    assert abs(wd.get("s", 0.0) - 1 / 6) < 1e-3
    assert wd.get("p", 0.0) < 1e-3

    chk_e = tmp_path / "frac_energy_ordered.chk"
    mol_e, mf_e = run_frac_occ_scf(
        cfg, str(chk_e), nelec_open=3, n_open=n_open)  # legacy energy order
    assert mf_e.converged  # it converges -- to the WRONG ensemble:
    frac_cols = np.where(np.isclose(np.asarray(mf_e.mo_occ), 0.5))[0]
    s = mf_e.get_ovlp()
    names, w = shell_weights(mol_e, mf_e.mo_coeff[:, frac_cols], s)
    assert w[:, names.index("p")].mean() > 0.3  # 4p-contaminated block
