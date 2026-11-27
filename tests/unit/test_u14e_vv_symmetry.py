"""The D2h column algebra of the V&V library (scripts/u14e/vv/symmetry.py)
on the frozen carbon D2h twin: the impurity is the S-metric dominant
share, the cleanup is the identity on pure columns and refuses real
mixing, the per-irrep orthonormalization preserves purity and span, and
the core/irrep intersection reaches its tolerances in closed form or
refuses a vanished column (no SCF; tests/carbon_reference_d2h.py)."""

import numpy as np
import pytest

from atomterm.irreps import orb_irreps
from atomterm.windows import select_shell_ao

from vv import symmetry

pytestmark = pytest.mark.fast

CORE = 2  # the frozen carbon reference's closed shell: 1s, 2s


@pytest.fixture(scope="module")
def twin(carbon_d2h):
    """(mol, S, the core MOs, the 2p AO columns, the first d AO that
    shares the core's irrep) -- every column irrep-pure by symmetry."""
    _, mol, mf, act_idx = carbon_d2h
    mo = np.asarray(mf.mo_coeff)
    core = mo[:, :CORE]
    p = np.eye(mol.nao_nr())[:, act_idx]
    d_all = np.eye(mol.nao_nr())[:, select_shell_ao(mol, "3d")]
    d_same = next(d_all[:, j] for j in range(d_all.shape[1])
                  if orb_irreps(mol, d_all[:, j:j + 1])
                  == orb_irreps(mol, core[:, :1]))
    return mol, mf.get_ovlp(), core, p, d_same


def _unit(v, s):
    return v / np.sqrt(float(v @ s @ v))


def test_impurity_is_the_s_metric_dominant_share(twin):
    # two columns in different irreps are S-orthogonal; a mixture's
    # reported impurity is the smaller share of its S-metric norm --
    # a coordinate-only count would weight AO rows, not the metric
    mol, s, _, p, _ = twin
    u1, u2 = p[:, 0], p[:, 1]
    assert abs(float(u1 @ s @ u2)) < 1e-12
    n1, n2 = float(u1 @ s @ u1), float(u2 @ s @ u2)
    for a, b in ((1.0, 1e-3), (1.0, 0.5), (0.3, 0.7)):
        _, (imp,) = symmetry.symmetrize_columns(mol, s, (a * u1 + b * u2)[:, None],
                                                tol=1.0)
        assert abs(imp - min(a * a * n1, b * b * n2) / (a * a * n1 + b * b * n2)) < 1e-12


def test_cleanup_is_the_identity_on_pure_columns(twin):
    mol, s, _, p, _ = twin
    fixed, imp = symmetry.symmetrize_columns(mol, s, p)
    assert max(imp) < 1e-14
    assert np.abs(fixed - np.column_stack([_unit(p[:, j], s) for j in range(3)])).max() < 1e-12


def test_cleanup_removes_a_tail_and_refuses_real_mixing(twin):
    mol, s, _, p, _ = twin
    taint = p.copy()
    taint[:, 0] = p[:, 0] + 1e-3 * p[:, 1]
    fixed, imp = symmetry.symmetrize_columns(mol, s, taint)
    # the tail's S-metric share: AO columns are S-unit here, so a
    # 1e-3 tail carries exactly 1e-6 of the norm
    assert abs(imp[0] - 1e-6 / (1.0 + 1e-6)) < 1e-14
    assert np.abs(fixed[:, 0] - _unit(p[:, 0], s)).max() < 1e-12
    with pytest.raises(ValueError, match="irrep impurity"):
        symmetry.symmetrize_columns(mol, s, taint, tol=1e-12)


def test_per_irrep_orthonormalization_keeps_purity_and_span(twin):
    # a same-irrep pair that is not S-orthogonal is Loewdin'd within
    # its own group; the other irreps' columns are only normalized
    mol, s, core, p, _ = twin
    block = np.column_stack([p[:, 0], core[:, 0],
                             core[:, 0] + 0.3 * core[:, 1]])
    before = np.concatenate([orb_irreps(mol, p[:, :1]), orb_irreps(mol, core)])
    assert np.abs(block[:, 1:].T @ s @ block[:, 1:] - np.eye(2)).max() > 1e-2
    fixed, loss = symmetry.orthonormal_per_irrep(mol, s, block)
    assert loss < 1e-13
    assert np.abs(fixed.T @ s @ fixed - np.eye(3)).max() < 1e-12
    assert np.array_equal(orb_irreps(mol, fixed), before)
    assert np.abs(fixed[:, 0] - _unit(p[:, 0], s)).max() < 1e-12
    # the group still spans exactly the two core MOs' space (a
    # unitary overlap between the two orthonormal bases)
    assert np.linalg.svd(core.T @ s @ fixed[:, 1:], compute_uv=False).min() > 1 - 1e-10


def test_per_irrep_orthonormalization_reports_a_cross_irrep_tail(twin):
    # a tail in another irrep is NOT rotated away (the treatment never
    # leaves a group); it is measured as the purity loss
    mol, s, _, p, _ = twin
    block = p.copy()
    block[:, 0] = p[:, 0] + 1e-3 * p[:, 1]
    fixed, loss = symmetry.orthonormal_per_irrep(mol, s, block)
    assert abs(loss - 1e-6 / (1.0 + 1e-6)) < 1e-14        # the same share
    assert np.abs(fixed[:, 1] - _unit(p[:, 1], s)).max() < 1e-12


def test_purify_reaches_the_core_irrep_intersection(twin):
    # a column of the core's own irrep carrying a core tail: one closed
    # form delivers core-orthogonal, irrep-pure, orthonormal
    mol, s, core, p, d_same = twin
    block = np.column_stack([p[:, 0], p[:, 1], _unit(d_same, s) + 1e-5 * core[:, 1]])
    assert np.abs(core[:, 1] @ s @ block[:, 2]) > 1e-6     # the premise: a core tail
    fixed, rep = symmetry.purify_window(mol, s, core, block)
    assert rep["passes"] == 1
    assert rep["core_res"] < 1e-10 and rep["irrep_imp"] < 1e-10
    assert rep["ortho_dev"] < 1e-10
    assert np.abs(core.T @ s @ fixed).max() < 1e-10
    assert np.abs(fixed.T @ s @ fixed - np.eye(3)).max() < 1e-10
    assert list(orb_irreps(mol, fixed)) == list(orb_irreps(mol, p[:, :2])) + list(orb_irreps(mol, core[:, :1]))


def test_purify_matches_the_independent_core_projection(twin):
    # on a core-lean column of the core's irrep the closed form equals
    # the plain S-metric core projection of that column
    mol, s, core, _, d_same = twin
    raw = np.column_stack([_unit(d_same, s) + 0.5 * core[:, 0]])
    fixed, rep = symmetry.purify_window(mol, s, core, raw)
    w = d_same - core @ np.linalg.solve(core.T @ s @ core, core.T @ s @ d_same)
    assert np.abs(fixed[:, 0] - _unit(w, s)).max() < 1e-10
    assert rep["core_res"] < 1e-10


def test_purify_refuses_a_column_inside_the_core_span(twin):
    mol, s, core, p, _ = twin
    inside = core[:, 0] + 0.1 * core[:, 1]
    with pytest.raises(symmetry.VanishedColumnError, match="vanished into the core"):
        symmetry.purify_window(mol, s, core, np.column_stack([p[:, 0], inside]))


def test_purify_refuses_genuinely_mixed_input(twin):
    mol, s, core, p, _ = twin
    mixed = p[:, 0] + 0.1 * p[:, 1]
    with pytest.raises(ValueError, match="irrep impurity"):
        symmetry.purify_window(mol, s, core, np.column_stack([p[:, 2], mixed]))
