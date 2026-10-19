"""Tests of atomterm.analysis: term-symbol helpers and configuration labels."""

import numpy as np

from atomterm.analysis import (
    configuration_summary,
    quantum_number,
    spin_quantum_number,
    term_symbol,
)

import pytest
pytestmark = pytest.mark.fast
from atomterm.detbasis import determinant_basis


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
