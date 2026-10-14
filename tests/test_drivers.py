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


def _s_ortho_err(mol_or_mf_mo, s):
    mo, = [np.asarray(mol_or_mf_mo)]
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


def test_build_full_mo_blocks_orthonormal_and_span_preserving():
    """The ticket-09 regression guard, on both atoms.

    Si is the case that broke: raw 3p AOs overlap the canonical 2p core
    MOs by ~0.3, so the block-safe construction has real work to do.
    """
    for cfg, patterns in [
        (AtomConfig(), ["2p"]),
        (AtomConfig("Si", "ccpvdz", spin=2), ["3p"]),
    ]:
        chk = f"/tmp/opencode/test_drivers_{cfg.symbol}.chk"
        mol, mf = run_rohf(cfg, chk)
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


def test_frac_occ_scf_open_block_ordering():
    """Average-of-configuration reference for Si: the fractional block must
    land on the 3p shell -- closed [Ne]+3s below (the ordering assumption
    of run_frac_occ_scf, first nontrivial atom after carbon).

    Runs twice deliberately: first call = fresh SCF, second = chkfile cache
    hit. Both returned objects must carry the ensemble occupation hook --
    re-deriving occupations from energies must reproduce the stored mo_occ
    (the ticket-10 finding: the hook used to be missing on cache hits).
    """
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    path = "/tmp/opencode/test_drivers_Si_frac.chk"
    Path(path).unlink(missing_ok=True)
    for k in range(2):  # 1st pass: fresh; 2nd pass: cache hit
        mol_f, mf_f = run_frac_occ_scf(
            dataclasses.replace(cfg, spin=0), path, nelec_open=2, n_open=3)
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
