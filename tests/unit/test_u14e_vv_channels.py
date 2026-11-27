"""The HF-virtual channels of the V&V library (scripts/u14e/vv/
channels.py) on the frozen carbon reference: the virtual spectrum
reads each virtual's angular character (purity to the l(l+1) grid,
energy order), a channel's block is an m-complete, l-pure,
S-orthonormal radial set, a mixed or empty signal is refused or
reported, and a channel set beside a base window leaves the window
orthonormal and rotation-closed while the frozen ground stays the
same term at or below the base's energy (no SCF;
tests/carbon_reference.py)."""

import numpy as np
import pytest

from atomterm.windows import shell_weights

from vv import cases as ladder
from vv.channels import Channel, channel_block, virtual_spectrum
from vv.reference import Substrate
from vv.routes import free_casci
from vv.windows import window_from_shells, window_with_channels

pytestmark = pytest.mark.fast

SECTOR = (2, 0)


@pytest.fixture(scope="module")
def p_channel(vv_substrate):
    """The cc-pVDZ second p shell: the three p-type virtuals above
    +0.5 Eh (found by character and energy, as the inventory is
    declared; the fourth p-type virtual, at +0.19 Eh, is the open
    shell's empty 2p partner, not a channel)."""
    spec = virtual_spectrum(vv_substrate)
    cols = tuple(r["index"] for r in spec
                 if r["l"] == "p" and r["de"] > 0.5)
    assert len(cols) == 3                       # one p shell of virtuals
    # the layer names BOTH p contractions: on a contracted basis the
    # second shell's orbitals spread over the two contractions of one
    # primitive set (a radial declared on the 3p band alone captures
    # only 0.15..0.32 of them)
    return Channel("3p", cols, ("2p", "3p"), tail=0.5)   # the members are split by the open shell


def test_the_virtual_spectrum_reads_character_and_purity(vv_substrate):
    spec = virtual_spectrum(vv_substrate, below=2.0)
    occ = np.asarray(vv_substrate.mf.mo_occ)
    assert all(occ[r["index"]] == 0 for r in spec)          # virtuals only
    assert [r["de"] for r in spec] == sorted(r["de"] for r in spec)
    assert all(r["de"] > 0 for r in spec)
    # every virtual of this basis is l-pure: <L^2> sits on the grid
    for r in spec:
        l = "spdf".index(r["l"])
        assert abs(r["l2"] - l * (l + 1)) < 1e-3
        assert r["weight"] > 0.999
    assert {r["l"] for r in spec} == {"s", "p", "d"}
    assert all(r["r2"] > 0 for r in spec)


def test_a_channel_block_is_a_pure_complete_orthonormal_radial_set(
        vv_substrate, p_channel):
    mol, mf = vv_substrate.mol, vv_substrate.mf
    s = mf.get_ovlp()
    block, lines = channel_block(vv_substrate, p_channel)
    assert block.shape == (mol.nao, 3)
    assert np.abs(block.T @ s @ block - np.eye(3)).max() < 1e-12
    _, w = shell_weights(mol, block, s)
    assert np.abs(w[:, 1] - 1.0).max() < 1e-10                 # l-pure
    assert any("m-averaged radial" in x for x in lines)
    # the radial lies on the declared bands' rows, and it carries
    # the signal: each virtual is captured
    rows = [i for i, lab in enumerate(mol.ao_labels())
            if lab.split()[-1].startswith(("2p", "3p"))]
    off = [i for i in range(mol.nao) if i not in rows]
    assert np.abs(block[off]).max() == 0.0
    virt = np.asarray(mf.mo_coeff)[:, list(p_channel.columns)]
    capture = np.sum((block.T @ s @ virt) ** 2, axis=0)
    assert capture.min() > 0.9


def test_a_fused_shell_declared_on_one_construction_is_not_captured(
        vv_substrate, p_channel):
    # the lesson behind the fixture's layer: the 3p band alone
    # carries only part of the second p shell (observed 0.15..0.32)
    mol, mf = vv_substrate.mol, vv_substrate.mf
    s = mf.get_ovlp()
    narrow = Channel("3p", p_channel.columns, ("3p",), tail=0.5)
    block, _lines = channel_block(vv_substrate, narrow)
    virt = np.asarray(mf.mo_coeff)[:, list(p_channel.columns)]
    capture = np.sum((block.T @ s @ virt) ** 2, axis=0)
    assert capture.max() < 0.5


def test_an_s_channel_is_the_virtual_s_s_sector(vv_substrate):
    # l-pure BY CONSTRUCTION: a C1 reference's s orbital carries a
    # small l>0 admixture and the closure meter is first order in its
    # amplitude -- the channel takes the s sector only
    mol, mf = vv_substrate.mol, vv_substrate.mf
    s = mf.get_ovlp()
    j = next(r["index"] for r in virtual_spectrum(vv_substrate)
             if r["l"] == "s")
    raw = np.asarray(mf.mo_coeff)[:, j]
    block, lines = channel_block(vv_substrate, Channel("4s", (j,), ("3s",)))
    s_rows = np.array([lab.split()[-1].endswith("s")
                       for lab in mol.ao_labels()])
    expect = raw * s_rows
    expect = expect / np.sqrt(float(expect @ s @ expect))
    assert np.abs(block[:, 0] - expect).max() == 0.0
    assert np.abs(block[~s_rows]).max() == 0.0
    assert float(block[:, 0] @ s @ block[:, 0]) == pytest.approx(1.0)
    _, w = shell_weights(mol, block, s)
    assert np.abs(w[0, 0] - 1.0) < 1e-12
    assert any("s sector" in x for x in lines)


def test_a_raw_c1_s_orbital_is_not_closed_but_its_s_sector_is(vv_substrate):
    # the measured defect and its cure, on the frozen reference: a
    # 5e-06 l>0 admixture moves the closure meter to ~1e-03; the s
    # sector closes at roundoff
    from atomterm.operators import window_closure

    mol, mf = vv_substrate.mol, vv_substrate.mf
    j = next(r["index"] for r in virtual_spectrum(vv_substrate)
             if r["l"] == "s")
    raw = np.asarray(mf.mo_coeff)[:, j:j + 1]
    block, _lines = channel_block(vv_substrate, Channel("4s", (j,), ("3s",)))
    cl_raw = max(window_closure(mol, raw).values())
    cl_pure = max(window_closure(mol, block).values())
    assert cl_raw > 100 * max(cl_pure, 1e-14)
    assert cl_pure < 1e-12


def test_a_signal_without_pure_columns_is_refused(vv_substrate, p_channel):
    impure = Channel("3p", p_channel.columns, ("3p",), pure=1.01)
    with pytest.raises(ValueError, match="none of its columns"):
        channel_block(vv_substrate, impure)
    mixed = Channel("3d", p_channel.columns, ("3d",))
    with pytest.raises(ValueError, match="carries >=0.97 d-weight"):
        channel_block(vv_substrate, mixed)


def test_a_channel_beside_a_base_window(vv_substrate, p_channel):
    base = ladder.visible_window(vv_substrate, (("2p", 3),), SECTOR)
    block, lines = channel_block(vv_substrate, p_channel)
    built = window_with_channels(vv_substrate, base, [(block, lines)])
    space = built.space
    assert (space.ncore, space.ncas, space.sector) == (2, 6, SECTOR)
    assert space.certify(vv_substrate.mf) < 1e-12
    cl = space.closure(vv_substrate.mol)
    assert max(cl["lz"], cl["lperp"]) < 1e-12
    assert built.report[:1] == base.report and lines[0] in built.report
    # the base's columns come first, untouched
    assert np.abs(space.active[:, :3] - base.space.active).max() < 1e-14
    # physics: the frozen ground stays the term, at or below the base
    a = free_casci(vv_substrate, base.space).price().states[0]
    b = free_casci(vv_substrate, space).price().states[0]
    assert abs(b.l2 - 2.0) < 1e-9 and abs(b.s2 - 2.0) < 1e-9
    assert b.state.e_ref <= a.state.e_ref + 1e-10
