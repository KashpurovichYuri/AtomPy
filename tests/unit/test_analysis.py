"""Tests of atomterm.analysis: term-symbol helpers and configuration labels."""

import numpy as np
import pytest

from atomterm.analysis import (
    configuration_summary,
    quantum_number,
    spin_quantum_number,
    term_symbol,
)
from atomterm.detbasis import determinant_basis

pytestmark = pytest.mark.fast


def test_quantum_number_round_trips_measured_momenta():
    """quantum_number inverts l(l+1) and term_symbol encodes (2S+1)L."""
    assert quantum_number(2.0) == 1
    assert quantum_number(6.0) == 2
    assert quantum_number(0.0) == 0
    assert quantum_number(2.0 + 1e-9) == 1
    assert term_symbol(2.0, 2.0) == "3P"
    assert term_symbol(0.0, 6.0) == "1D"
    assert term_symbol(0.0, 0.0) == "1S"
    assert term_symbol(2.0, 12.0) == "3F"
    assert term_symbol(0.0, 20.0) == "1G"


def test_letter_table_follows_spectroscopic_convention():
    """The sequence skips J, then the second P and second S: after
    O (L = 11) comes Q (L = 12) -- the review-wave catch: the table
    reused P at L = 12. f^3 d^1 reaches L = 8 (the uranium 5L ground),
    f^4 terms reach L = 12."""
    assert term_symbol(0.0, 72.0) == "1L"      # L = 8, the U ground
    assert term_symbol(0.0, 110.0) == "1N"     # L = 10
    assert term_symbol(0.0, 132.0) == "1O"     # L = 11
    assert term_symbol(0.0, 156.0) == "1Q"     # L = 12: not P
    assert term_symbol(0.0, 380.0) == "1Y"     # L = 19, table end - 1
    with pytest.raises(ValueError, match="beyond the letter table"):
        term_symbol(0.0, 99 * 100.0)           # L = 99: refuse, not wrap


def test_half_integer_spin_labels():
    """Odd-electron states: S = 1/2, 3/2 give multiplicities 2 and 4.

    Rounding S to the integer grid relabels doublets/quartets as singlets/
    triplets with 1e-8-level <S^2> noise deciding which -- first exposed
    by Sc (ticket 11), the project's first doublet.
    """
    assert spin_quantum_number(0.75) == 0.5
    assert spin_quantum_number(0.75 + 1e-9) == 0.5
    assert spin_quantum_number(0.75 - 1e-9) == 0.5
    assert spin_quantum_number(3.75) == 1.5
    assert spin_quantum_number(2.0) == 1.0
    assert term_symbol(0.75, 6.0) == "2D"   # Sc ground term
    assert term_symbol(0.75 + 1e-9, 6.0) == "2D"
    assert term_symbol(0.75 - 1e-9, 6.0) == "2D"
    assert term_symbol(3.75, 12.0) == "4F"


def test_configuration_summary_reads_shell_configurations(carbon):
    """On a raw-shell window the label IS the determinant's AO occupation."""
    _, mol, mf, act_idx = carbon
    s = mf.get_ovlp()
    c_win = np.eye(mol.nao_nr())[:, act_idx]
    dets = determinant_basis(3, 2, 0)

    # ground 3P determinant |px py>: pure p^2
    ci = np.zeros(len(dets))
    ci[dets.index((0, 1))] = 1.0
    labels = configuration_summary(mol, s, [ci], c_win, na=2, nb=0)
    assert labels[0].startswith("p^2.00")
    assert "(|c|=1.00)" in labels[0]

    # a 50/50 mix of two determinants: both components reported, p^2 each
    ci2 = np.zeros(len(dets))
    ci2[0] = ci2[1] = 2 ** -0.5
    labels = configuration_summary(mol, s, [ci2], c_win, na=2, nb=0)
    assert labels[0].count("p^2.00") == 2

    # an s-admixed window column shows up as s-mass in the label
    c_pol = c_win.copy()
    c_pol[:, 0] = 0.9 * c_win[:, 0] + 0.3 * mf.mo_coeff[:, 1]  # 2s admixture
    labels = configuration_summary(mol, s, [ci], c_pol, na=2, nb=0)
    assert "s^" in labels[0]


def test_configuration_summary_maps_pyscf_ci_order(carbon):
    """The ordering bridge inside configuration_summary, on the smallest
    sector where the orders differ (norb=4, na=2: PySCF's ascending bit
    strings give (0,1),(0,2),(1,2),(0,3),(1,3),(2,3); detbasis is
    lexicographic (0,1),(0,2),(0,3),...). A vector on PySCF address 3 =
    detbasis (0,3) = the 2s+2pz determinant must label s^1 p^1 -- without
    the bridge it would read detbasis index 3 = (1,2) = p^2."""
    _, mol, mf, act_idx = carbon
    from atomterm.windows import select_valence_ao

    sp = select_valence_ao(mol, ["2s", "2p"])
    c_sp = np.eye(mol.nao_nr())[:, sp]
    dets = determinant_basis(4, 2, 0)
    assert len(dets) == 6

    ci = np.zeros(6)
    ci[3] = 1.0                                   # PYSCF address 3
    labels = configuration_summary(mol, mf.get_ovlp(), [ci], c_sp, 2, 0)
    assert "s^1.00" in labels[0] and "p^1.00" in labels[0]


def test_analyze_rejects_a_complex_hamiltonian():
    """The reality guard: a Hamiltonian with non-noise imaginary part is
    a construction error upstream, not physics to diagonalize."""
    import pytest
    from atomterm.analysis import analyze

    n = 4
    h = np.eye(n) + 1e-3j * np.ones((n, n))
    with pytest.raises(ValueError, match="should be real"):
        analyze(h, np.eye(n), np.eye(n), np.eye(n), "test")


def test_report_carries_the_sz_spectrum_line():
    """The delivered Sz matrix enters the report as one compressed line:
    its spectrum IS the sector's M_S (Sz = M_S*1 exactly), so the line
    documents the quantum number and certifies nothing. No sz matrix ->
    no line (analyze stays backward compatible)."""
    from atomterm.analysis import analyze, format_report
    from atomterm.detbasis import determinant_basis, slater_condon_matrix
    from atomterm.operators import sz_operator

    def report(norb, na, nb, pass_sz):
        dets = determinant_basis(norb, na, nb)
        sz_mb = slater_condon_matrix(dets, one_body=sz_operator(norb))
        h = slater_condon_matrix(dets, one_body=np.eye(2 * norb))
        kw = {"sz": sz_mb} if pass_sz else {}
        res = analyze(h, sz_mb, 2.0 * np.eye(len(dets)),
                      0.75 * np.eye(len(dets)),
                      f"Ms={(na - nb) / 2:+g}", **kw)
        return res, format_report(res, "test")

    # Ms = +1 sector: three determinants, all at M_S = +1
    res, text = report(3, 2, 0, pass_sz=True)
    assert res.sz_spectrum == [1.0] * 3
    assert "Sz eigenvalues (M_S content): [+1] x3" in text
    assert "exactly M_S·1 in a fixed sector" in text

    # Ms = 0 sector: the zero prints unsigned
    _, text0 = report(3, 1, 1, pass_sz=True)
    assert "Sz eigenvalues (M_S content): [0] x9" in text0

    # without the matrix the line is omitted, not broken
    _, text_none = report(3, 2, 0, pass_sz=False)
    assert "Sz eigenvalues" not in text_none


def test_shell_occupations_electron_scale(carbon):
    """n_l partitions the active electrons EXACTLY: a p² determinant
    over the raw 2p window gives n_p = 2; a column rotated between 2s
    and a 2p member (exactly S-orthogonal, different l) splits the
    pair's electrons s:1 / p:1. The S-CONTRACTION is load-bearing: a
    same-l non-orthogonal mix (2s carrying a 1s tail -- cc-pVDZ s
    shells overlap) must still split by the S-metric weights, which
    the bare AO diagonal of C gamma C^dag gets wrong."""
    from atomterm.analysis import shell_occupations
    from atomterm.windows import select_valence_ao

    _, mol, mf, act_idx = carbon
    s = mf.get_ovlp()
    c_p = np.eye(mol.nao_nr())[:, act_idx]          # the three 2p AOs
    dets = determinant_basis(3, 1, 1)
    k = dets.index((0, 3))                          # p0-alpha p0-beta
    ci = np.zeros(len(dets))
    ci[k] = 1.0

    n = shell_occupations(mol, s, c_p, ci, 1, 1)
    assert n["p"] == pytest.approx(2.0, abs=1e-12)
    assert n["s"] == pytest.approx(0.0, abs=1e-12)
    assert sum(n.values()) == pytest.approx(2.0, abs=1e-12)

    # the mixed first column takes the pair: s:1, p:1
    mixed = c_p.copy()
    v = np.zeros(mol.nao_nr())
    v[select_valence_ao(mol, ["2s"])[0]] = np.sqrt(0.5)
    v[act_idx[0]] = np.sqrt(0.5)
    mixed[:, 0] = v
    n2 = shell_occupations(mol, s, mixed, ci, 1, 1)
    assert n2["s"] == pytest.approx(1.0, abs=1e-12)
    assert n2["p"] == pytest.approx(1.0, abs=1e-12)

    # same-l non-orthogonality: an s-column with a 1s tail, S-renormed,
    # crossed with 2p0 -- the split stays s:1 / p:1 (the S weights),
    # while the bare-diagonal reading would drift
    vs = np.zeros(mol.nao_nr())
    vs[select_valence_ao(mol, ["2s"])[0]] = 1.0
    vs[select_valence_ao(mol, ["1s"])[0]] = 0.3
    vs /= np.sqrt(vs @ s @ vs)                       # S-normalized
    mixed2 = c_p.copy()
    w = np.zeros(mol.nao_nr())
    w[:] = vs / np.sqrt(2)
    w[act_idx[0]] = np.sqrt(0.5)                     # s-block and p0
    mixed2[:, 0] = w
    n3 = shell_occupations(mol, s, mixed2, ci, 1, 1)
    assert n3["s"] == pytest.approx(1.0, abs=1e-12)
    assert n3["p"] == pytest.approx(1.0, abs=1e-12)
    assert sum(n3.values()) == pytest.approx(2.0, abs=1e-12)
