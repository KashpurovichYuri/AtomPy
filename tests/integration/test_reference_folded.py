"""The folded reference (:func:`atomterm.reference.run_folded`): the
confined core frozen, the valence ROHF self-consistent in its field --
the lane's reference class.  Two core partitions on the carbon fixture
(1s-only: a closed valence pair + the p opens -- the Ti shape;
1s2s: a pure-open valence); one real valence SCF per partition
(seconds)."""

import numpy as np
import pytest

from atomterm.core import closed_shell_core_energy, confined_core_orbitals
from atomterm.reference import _closed_field, run_folded
from atomterm.windows import lowdin, select_shell_ao, select_span, span_projector

pytestmark = pytest.mark.fast


def _span(mol, mf, patterns):
    rows: list[int] = []
    for p in patterns:
        rows.extend(select_shell_ao(mol, p))
    return lowdin(np.eye(mol.nao_nr())[:, sorted(set(rows))],
                  mf.get_ovlp())


def test_folded_reference_freezes_confined_core_bitwise(carbon):
    """The frozen block: the confined solve's columns, BIT-identical,
    first in the assembled MO matrix."""
    _, mol, mf, _ = carbon
    span = _span(mol, mf, ("1s", "2s"))
    _, fm = run_folded(mf, span, 2, ("l0",), ("l1",))
    c_conf, _ = confined_core_orbitals(mf, span, 2)
    assert fm.converged
    assert np.array_equal(fm.mo_coeff[:, :2], c_conf)


def test_folded_reference_sits_above_full_rohf(carbon):
    """The variational bound: the frozen core constrains the full
    same-filling problem (both partitions)."""
    _, mol, mf, _ = carbon
    for shells, npairs in ((("1s",), 1), (("1s", "2s"), 2)):
        _, fm = run_folded(mf, _span(mol, mf, shells), npairs,
                           ("l0",), ("l1",))
        assert fm.e_tot >= mf.e_tot - 1e-9


def test_folded_reference_energy_is_the_explicit_bookkeeping(carbon):
    """E = E_conf + Tr[D_v h_eff] + 1/2 Tr[D_v V[D_v]], recomputed
    from the returned object alone (no solve internals)."""
    _, mol, mf, _ = carbon
    span = _span(mol, mf, ("1s",))
    _, fm = run_folded(mf, span, 1, ("l0",), ("l1",))
    c_conf, _ = confined_core_orbitals(mf, span, 1)
    dm = fm.make_rdm1()
    if dm.ndim == 3:            # the ROHF convention: the per-spin pair
        dm = dm.sum(axis=0)     # (the spin-summed density)
    dm_conf = 2.0 * (c_conf @ c_conf.conj().T)
    d_v = dm - dm_conf
    h_eff = mf.get_hcore() + _closed_field(mf, dm_conf)
    e = (closed_shell_core_energy(mf, c_conf)
         + float(np.einsum("ij,ji->", h_eff, d_v).real)
         + 0.5 * float(np.einsum("ij,ji->", _closed_field(mf, d_v),
                                 d_v).real))
    assert fm.e_tot == pytest.approx(e, abs=1e-9)


def test_folded_reference_filling_and_orthonormality(carbon):
    """The assembled MOs: S-orthonormal; the frozen pairs first; the
    electron count exact; the opens sit in the open span."""
    _, mol, mf, _ = carbon
    span = _span(mol, mf, ("1s",))
    _, fm = run_folded(mf, span, 1, ("l0",), ("l1",))
    s = mf.get_ovlp()
    assert np.allclose(fm.mo_coeff.conj().T @ s @ fm.mo_coeff,
                       np.eye(mol.nao_nr()), atol=1e-10)
    assert (fm.mo_occ == 2).sum() == 2   # 1 frozen + 1 closed valence
    assert (fm.mo_occ == 1).sum() == 2   # the p opens
    assert fm.mo_occ.sum() == mol.nelectron
    opens = np.flatnonzero(fm.mo_occ == 1)
    w_open = np.einsum("ai,ab,bi->i", fm.mo_coeff.conj()[:, opens],
                       span_projector(s, select_span(mol, "l1")),
                       fm.mo_coeff[:, opens]).real
    assert np.all(w_open > 0.99)


def test_folded_reference_pick_rule_is_declared(carbon):
    """The pick rules are DECLARED per atom-gauge (both certified
    forms, neither universal -- and on the uncontracted twins the
    folded landscape is draw-chaotic, so the converged-filling
    certificate judges every draw).  An unknown rule is rejected; the
    three-level form converges on the carbon fixture to its own
    certified field (a DIFFERENT basin than the gates default's,
    measured -- both hold the declared filling and the variational
    bound)."""
    _, mol, mf, _ = carbon
    span = _span(mol, mf, ("1s",))
    with pytest.raises(ValueError, match="picks"):
        run_folded(mf, span, 1, ("l0",), ("l1",), picks="energy-only")
    _, fm_a = run_folded(mf, span, 1, ("l0",), ("l1",))
    _, fm_b = run_folded(mf, span, 1, ("l0",), ("l1",),
                         picks="three-level")
    for fm in (fm_a, fm_b):
        assert fm.converged
        assert fm.e_tot >= mf.e_tot - 1e-9   # the variational bound
        assert (fm.mo_occ == 2).sum() == 2   # 1 frozen + 1 closed valence
        assert (fm.mo_occ == 1).sum() == 2   # the p opens
    # the measured fields (the gates default is the deeper basin here)
    assert fm_a.e_tot == pytest.approx(-37.3925371924919, abs=1e-9)
    assert fm_b.e_tot == pytest.approx(-35.188412949352745, abs=1e-9)
