"""cells.readout_lane: the lane store's readout twin -- the store's
own re-solve reproduces the stored lane record (energies strict, the
meters within the re-solve noise band) and the blocks price through
``perturb.nevpt2_rows_on``; a drifted store refuses.  All on the
FROZEN D2h carbon twin (NO SCF): the CI lane cell is written once per
module into a repo-local ``*.tmp`` store (removed after)."""

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from atomterm.cases import CAMPAIGN
from atomterm.cells import CellKey, load_lane_cell, readout_lane, \
    run_lane_cell

pytestmark = pytest.mark.fast

WORKDIR = Path("elements") / "C" / "lane-twin.tmp"
ADDR = "c_2p@1s2s_2-0_casci_ao_-_pur_avg"


@pytest.fixture(scope="module")
def lane_store(carbon_d2h):
    """The CI lane cell written once into the repo-local tmp store."""
    _, mol, mf, _ = carbon_d2h
    spec = CAMPAIGN["C"]
    key = CellKey.from_address(ADDR)
    WORKDIR.mkdir(parents=True, exist_ok=True)
    run_lane_cell(key, mol, mf, spec.term,
                  shell_bands=spec.shell_bands,
                  shell_radials=spec.shell_radials, workdir=WORKDIR)
    yield key
    shutil.rmtree(WORKDIR)


def test_lane_readout_twin_reproduces_the_store(lane_store, carbon_d2h):
    """Each block's re-solved roots equal the stored record's own
    slice (per-block anchoring; the meters assertion's silence is
    the other half of the invariant)."""
    _, mol, mf, _ = carbon_d2h
    payload, blocks = readout_lane(lane_store, mol, mf, workdir=WORKDIR)
    stored = load_lane_cell(lane_store, WORKDIR)
    i0 = 0
    for mc in blocks.values():
        mult = mc.fcisolver.nroots
        assert np.allclose(np.asarray(mc.e_cas).ravel(),
                           np.asarray(stored["energies"][i0:i0 + mult]),
                           atol=1e-11, rtol=0.0)
        i0 += mult
    assert i0 == len(stored["energies"])


def test_lane_blocks_price_through_the_seam(lane_store, carbon_d2h):
    """Each block prices through nevpt2_rows_on: the class sum is the
    row's e_corr and e_total = e_ref + e_corr (the payload
    identities, at unit scale)."""
    from atomterm.perturb import nevpt2_rows_on

    _, mol, mf, _ = carbon_d2h
    _, blocks = readout_lane(lane_store, mol, mf, workdir=WORKDIR)
    n = 0
    for mc in blocks.values():
        for row in nevpt2_rows_on(mc):
            assert row.e_corr == pytest.approx(
                sum(c.energy for c in row.classes), abs=1e-12)
            assert row.e_total == pytest.approx(
                row.e_ref + row.e_corr, abs=1e-14)
            n += 1
    assert n == 3


def test_lane_readout_refuses_drift(lane_store, carbon_d2h):
    """A tampered energy record refuses the whole call (fail-fast,
    before any PT2)."""
    from pyscf import lib

    from atomterm.cells import cell_chk

    _, mol, mf, _ = carbon_d2h
    stored = load_lane_cell(lane_store, WORKDIR)
    bad = list(stored["energies"])
    bad[0] += 1e-6
    lib.chkfile.save(str(cell_chk(lane_store, WORKDIR)),
                     f"atomterm/{ADDR}/energies", bad)
    with pytest.raises(ValueError, match="lane re-solve refusal"):
        readout_lane(lane_store, mol, mf, workdir=WORKDIR)
