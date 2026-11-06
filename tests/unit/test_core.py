"""Unit pins of atomterm.core: the closed core ON a fold span (the
multi-column fold's mini-Roothaan construction) and the closed-shell
core energy's two forms -- the columnwise bookkeeping and the
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
    re-computation; the mini-Roothaan never runs)."""
    _, mol, mf, _ = carbon
    c_span = _span(mol, mf, ("1s", "2s"))
    n_pairs = closed_electrons(("1s", "2s")) // 2
    assert n_pairs == c_span.shape[1]
    c_core, report = confined_core_orbitals(mf, c_span, n_pairs)
    assert c_core is c_span
    assert report == {"cycles": 0, "comm": 0.0}


def test_confined_core_is_self_consistent_and_variational(carbon):
    """The mini-Roothaan's exact-by-construction identities on a
    MULTI-column realization (the 1s2s chemistry realized by all
    three s bands -- two pairs over a three-column span): the
    converged density COMMUTES with its own Fock (the
    self-consistency equation itself), the two energy forms agree
    (the RHF trace form vs the columnwise bookkeeping), and the
    construction is VARIATIONAL -- the bigger span's core sits at or
    below the 1:1 raw-span core's energy."""
    _, mol, mf, _ = carbon
    n_pairs = 2
    c_span = _span(mol, mf, ("1s", "2s", "3s"))
    assert c_span.shape[1] == 3
    c_core, report = confined_core_orbitals(mf, c_span, n_pairs)
    assert c_core.shape[1] == n_pairs
    assert report["comm"] < 1e-10 and report["cycles"] >= 1

    # the self-consistency equation, in the SPAN basis (the invariance
    # the Roothaan delivers: the occupied space is invariant under the
    # SPAN-RESTRICTED Fock -- the full-space Fock has off-span blocks
    # by construction).  The span basis is S-ORTHONORMAL: the
    # coefficient rep of the core columns rides the S metric
    # (``c^dag S c_core``, not the Euclidean ``c^dag c_core``).
    h = mf.get_hcore()
    eri = mol.intor("int2e", aosym="s1")
    sm = mf.get_ovlp()
    d_ao = 2.0 * c_core @ c_core.conj().T
    f_ao = (h
            + np.einsum("ls,mnls->mn", d_ao, eri, optimize=True)
            - 0.5 * np.einsum("ls,mlns->mn", d_ao, eri, optimize=True))
    f_sp = c_span.conj().T @ f_ao @ c_span
    occ_sp = c_span.conj().T @ sm @ c_core
    d_sp = 2.0 * occ_sp @ occ_sp.conj().T
    comm = f_sp @ d_sp - d_sp @ f_sp
    assert np.abs(comm).max() < 1e-9

    # the two energy forms (independent summation orders)
    e_cols = closed_shell_core_energy(mf, c_core)
    e_rhf = 0.5 * float(np.trace(d_ao @ (h + f_ao)))
    assert abs(e_cols - e_rhf) < 1e-9
    assert abs(e_cols - _density_form_energy(mol, mf, c_core)) < 1e-9

    # variational monotonicity: the superset span's core is not above
    # the 1:1 raw-span core ({1s, 2s} columns sit inside the span)
    c_11, _ = confined_core_orbitals(mf, _span(mol, mf, ("1s", "2s")),
                                     n_pairs)
    assert e_cols <= closed_shell_core_energy(mf, c_11) + 1e-12


def test_confined_core_energy_equals_the_density_form(carbon):
    """The projected-density certification on the 1:1 regime too: one
    formula, both regimes."""
    _, mol, mf, _ = carbon
    c_11 = _span(mol, mf, ("1s", "2s"))
    assert abs(closed_shell_core_energy(mf, c_11)
               - _density_form_energy(mol, mf, c_11)) < 1e-9


def test_confined_core_refuses_the_impossible(carbon):
    """The refusal: more pairs than the span holds columns."""
    _, mol, mf, _ = carbon
    c_span = _span(mol, mf, ("1s", "2s"))
    with pytest.raises(ValueError, match="holds columns"):
        confined_core_orbitals(mf, c_span, 3)
