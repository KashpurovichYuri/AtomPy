"""Unit pins of atomterm.core: the confined closed core on a fold
span (the multi-column fold's construction) and the closed-shell core
energy's two forms -- the columnwise bookkeeping and the
projected-density form (independent summation orders) -- on the
frozen carbon reference (no SCF)."""

import numpy as np
import pytest
from pyscf import ao2mo

from atomterm.core import closed_shell_core_energy, confined_core_orbitals
from atomterm.operators import build_operator_parts
from atomterm.windows import closed_electrons, lowdin, select_shell_ao

pytestmark = pytest.mark.fast


def _span(mol, mf, bands):
    """The Loewdin-orthonormal raw-AO span of the given bands (the
    fold-span construction, ``cells._fold_span``'s)."""
    rows = sorted(r for b in bands for r in select_shell_ao(mol, b))
    return lowdin(np.eye(mol.nao_nr())[:, rows], mf.get_ovlp())


def _density_form_energy(mol, mf, c):
    """The closed-shell energy of doubly-occupied ``c`` columns in the
    AO density form -- the INDEPENDENT construction: with
    ``D = c c^dag`` (the pair density, AO metric),

        E = 2 tr(D h) + sum_mnls D_mn D_ls (2 (mn|ls) - (ml|ns))

    (a different summation order than the columnwise J/K loop of
    :func:`closed_shell_core_energy`)."""
    d = c @ c.conj().T
    e = 2.0 * float(np.trace(d @ mf.get_hcore()))
    eri = mol.intor("int2e", aosym="s1")            # chemist (mn|ls)
    e += float(np.einsum("mn,ls,mnls->", d, d, 2.0 * eri, optimize=True))
    e -= float(np.einsum("mn,ls,mlns->", d, d, eri, optimize=True))
    return e


def test_confined_core_is_the_span_on_a_one_to_one_realization(carbon):
    """n_pairs == span width returns the span ITSELF (the certified
    1:1 raw-shell fold reproduced bit-exactly -- the identity, not a
    re-computation)."""
    _, mol, mf, _ = carbon
    c_span = _span(mol, mf, ("1s", "2s"))
    n_pairs = closed_electrons(("1s", "2s")) // 2
    assert n_pairs == c_span.shape[1]
    c_core, weights = confined_core_orbitals(mf, c_span, n_pairs)
    assert c_core is c_span
    assert np.all(weights == 1.0)


def test_confined_core_returns_the_closed_mos_on_the_whole_space(carbon):
    """The degenerate anchor: the whole AO space as the span confines
    to the closed MOs THEMSELVES (the projection is the identity), so
    the confined core equals the reference's closed block and the
    folded operators are IDENTICAL to the canonical-prefix fold --
    exact by construction, any deviation is a bug."""
    _, mol, mf, _ = carbon
    closed = [i for i, o in enumerate(mf.mo_occ) if o == 2.0]
    c_closed = np.asarray(mf.mo_coeff[:, closed])
    c_span = lowdin(np.eye(mol.nao_nr()), mf.get_ovlp())
    c_core, weights = confined_core_orbitals(mf, c_span, len(closed))
    assert np.allclose(c_core, c_closed, atol=1e-12)
    assert np.allclose(weights, 1.0, atol=1e-12)

    act = np.eye(mol.nao_nr())[:, select_shell_ao(mol, "2p")]
    a = build_operator_parts(mol, mf, act, len(closed), 2, 0,
                             c_core=c_closed)
    b = build_operator_parts(mol, mf, act, len(closed), 2, 0,
                             c_core=c_core)
    assert np.abs(a["h1"] - b["h1"]).max() < 1e-12


def test_confined_core_energy_equals_the_density_form(carbon):
    """The projected-density certification: on a MULTI-column
    realization (the 1s2s chemistry realized by all three s bands --
    two pairs over a three-column span) the columnwise bookkeeping
    over the confined columns equals the AO density form, an
    independent summation order; and the 1:1 columnwise form equals it
    too (one formula, both regimes)."""
    _, mol, mf, _ = carbon
    c_span = _span(mol, mf, ("1s", "2s", "3s"))
    assert c_span.shape[1] == 3
    c_core, weights = confined_core_orbitals(mf, c_span, 2)
    assert c_core.shape[1] == 2
    # the s bands capture C's closed pair to numerical completeness
    assert weights.min() > 0.99, weights
    e_cols = closed_shell_core_energy(mf, c_core)
    assert abs(e_cols - _density_form_energy(mol, mf, c_core)) < 1e-9

    c_11 = _span(mol, mf, ("1s", "2s"))
    assert abs(closed_shell_core_energy(mf, c_11)
               - _density_form_energy(mol, mf, c_11)) < 1e-9


def test_confined_core_refuses_the_ambiguous_cut_and_the_impossible(carbon):
    """The refusals: more pairs than columns; more pairs than the
    closed supply; a cut without a gap (the whole-space span at
    n_pairs=1 sees two closed MOs at weight 1 -- the fold's shells are
    ambiguous there, and the construction says so)."""
    _, mol, mf, _ = carbon
    c_span = _span(mol, mf, ("1s", "2s"))
    with pytest.raises(ValueError, match="holds columns"):
        confined_core_orbitals(mf, c_span, 3)
    whole = lowdin(np.eye(mol.nao_nr()), mf.get_ovlp())
    with pytest.raises(ValueError, match="not clean"):
        confined_core_orbitals(mf, whole, 1)
    with pytest.raises(ValueError, match="supplies"):
        confined_core_orbitals(
            mf, lowdin(np.eye(mol.nao_nr())[:, :8], mf.get_ovlp()), 6)
