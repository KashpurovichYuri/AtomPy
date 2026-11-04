"""The UHF reference class on a real SCF: run_uhf's convergence
discipline and the contamination meter on the campaign's smallest
atom, against the ROHF anchor at the same forced occupations."""

import numpy as np
import pytest

from atomterm.cases import CAMPAIGN
from atomterm.reference import ReferenceSpec, run_rohf, run_uhf, \
    uhf_contamination


@pytest.fixture(scope="module")
def carbon_references(tmp_path_factory):
    """(mf_rohf, mf_uhf) for the campaign's C/cc-pVDZ configuration:
    ROHF plain, UHF with the per-spin character forcing (1s/2s pairs,
    2p open) -- the same ReferenceSpec vocabulary both formalisms."""
    cfg = CAMPAIGN["C"].atom_config()
    ref = ReferenceSpec(occ_closed=("1s", "2s"), occ_open=("2p",))
    root = tmp_path_factory.mktemp("c_uhf")
    _, mf_r = run_rohf(cfg, str(root / "rohf.chk"))
    _, mf_u = run_uhf(cfg, str(root / "uhf.chk"), ref=ref)
    return mf_r, mf_u


def test_uhf_reference_converges_at_or_below_rohf(carbon_references):
    """UHF contains the ROHF solution space: converged, at or below
    the ROHF energy at the same forced occupations (equal when C does
    not break)."""
    mf_r, mf_u = carbon_references
    assert mf_u.converged
    assert mf_u.e_tot <= mf_r.e_tot + 1e-9


def test_uhf_reference_meter_reads(carbon_references):
    """The contamination meter on the converged reference: <S^2> at
    or above the pure-triplet floor (contamination only adds), one
    corresponding-orbital pair per beta electron inside (0, 1]."""
    _, mf_u = carbon_references
    s2, sigma = uhf_contamination(mf_u)
    assert s2 >= 2.0 - 1e-9
    assert sigma.size == 2
    assert np.all(sigma <= 1.0 + 1e-12)
    assert np.all(sigma > 0.0)


def test_uhf_reference_cache_trust(tmp_path_factory):
    """The chkfile round-trip: a second run_uhf on the cached chk
    loads it (no SCF) and returns the same energy -- the
    reproducibility unit, UHF-shaped (mo_coeff (2, nao, nao))."""
    cfg = CAMPAIGN["C"].atom_config()
    ref = ReferenceSpec(occ_closed=("1s", "2s"), occ_open=("2p",))
    chk = str(tmp_path_factory.mktemp("c_uhf_cache") / "uhf.chk")
    _, first = run_uhf(cfg, chk, ref=ref)
    _, second = run_uhf(cfg, chk, ref=ref)
    assert second.converged and second.e_tot == first.e_tot
    assert np.asarray(second.mo_coeff).shape[0] == 2
