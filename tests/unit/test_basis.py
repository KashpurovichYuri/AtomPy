"""Unit pins of atomterm.basis: the prior axis's constructions (the
density-block family) -- the metric eigenproblem against scipy's
independent generalized solve, the member-summed occupation-weighted
block, and series_window on synthetic ladders over the uncontracted
carbon substrate (planted radials, the rank-1 refusal, the closure
and capture certificates; no SCF anywhere)."""

import numpy as np
import pytest

from scipy.linalg import eigh as gen_eigh

from atomterm.basis import (SeriesSpec, metric_radials, run_series,
                            series_window, valence_block)
from atomterm.system import AtomConfig
from atomterm.windows import member_chains, select_shell_ao
from helpers import uncontracted_carbon_mole

pytestmark = pytest.mark.fast


# --------------------------------------------- the metric eigenproblem

def test_metric_radials_matches_the_independent_generalized_solve():
    """The transported construction against scipy's generalized eigh
    (B, s_r^-1) -- the capture-maximizing problem: the same top-k
    eigenvalues and projectors, the columns s_r-orthonormal.  The
    PLANTED-RADIAL anchor (series_window's tests) is the independent
    direction check -- two constructions sharing a flipped conjugation
    would cancel here and only the anchor would catch it."""
    rng = np.random.default_rng(3)
    n = 6
    a = rng.normal(size=(n, n))
    s_r = a @ a.T + np.eye(n)
    b = rng.normal(size=(n, n))
    block = 0.5 * (b + b.T) @ (0.5 * (b + b.T)).T
    v, lam = metric_radials(block, s_r, 3)
    assert np.abs(v.T @ s_r @ v - np.eye(3)).max() < 1e-9
    w = gen_eigh(block, np.linalg.inv(s_r), eigvals_only=True)[::-1]
    assert np.allclose(lam[:3], w[:3], atol=1e-10)
    q = gen_eigh(block, np.linalg.inv(s_r))[1][:, -3:][:, ::-1]
    v_ref = np.linalg.solve(s_r, q)        # the capture-maximizer form
    cross = v.T @ s_r @ v_ref              # subspace cosines
    assert np.abs(np.abs(cross) - np.eye(3)).max() < 1e-9


def test_metric_radials_recovers_the_planted_capture_maximizer():
    """The direction anchor on a nontrivial metric: a rank-1 block
    built from an s_r-normalized planted radial returns EXACTLY that
    radial (the capture-maximizer) -- the ``(block, s_r)`` flip would
    return the s_r^{-1}-distorted direction instead."""
    rng = np.random.default_rng(11)
    n = 4
    a = rng.normal(size=(n, n))
    s_r = a @ a.T + np.eye(n)
    r = rng.normal(size=n)
    r /= np.sqrt(r @ s_r @ r)
    v, lam = metric_radials(6 * np.outer(r, r), s_r, 1)
    assert np.allclose(np.abs(v[:, 0]), np.abs(r), atol=1e-10)
    assert lam[1] / lam[0] < 1e-12


def test_metric_radials_rank1_block_carries_one_radial():
    """The single-configuration premise: a rank-1 block's spectrum
    collapses after the leading eigenvalue -- the second radial of
    ONE configuration would be noise."""
    rng = np.random.default_rng(5)
    n = 5
    a = rng.normal(size=(n, n))
    s_r = a @ a.T + np.eye(n)
    u = rng.normal(size=n)
    _v, lam = metric_radials(np.outer(u, u), s_r, 1)
    assert lam[0] > 0
    assert lam[1] / lam[0] < 1e-12


# ------------------------------------------------------- the signal block

class _FakeMF:
    """A minimal SCF stand-in: columns + occupations (the block math
    needs nothing else)."""

    def __init__(self, c, occ):
        self.mo_coeff, self.mo_occ = c, occ


def test_valence_block_sums_members_weights_and_cuts_the_core():
    """The block's three contracts on a hand-built substrate: rows
    beyond ncol never enter, occupations weight the columns, and the
    member summation makes the block m-scalar (identical on every
    member's chain rows)."""
    mol = uncontracted_carbon_mole()
    s = mol.intor_symmetric("int1e_ovlp")
    chain, s_r = member_chains(mol, "2p", ("2p", "3p"), s)
    n = mol.nao_nr()
    # two valence columns: a planted radial replicated on the m-set
    # (occ 2) and a member-zero-only column (occ 1) -- the m-symmetry
    # detector
    radial = np.array([0.8, 0.2])
    c_sym = np.zeros((n, 1))
    for m in range(3):
        c_sym[chain[m], 0] = radial
    c_asym = np.zeros((n, 1))
    c_asym[chain[0], 0] = [0.3, 0.1]
    c = np.hstack([c_sym, c_asym])
    # one core column entirely inside the 1s/2s rows (past ncol=2 it
    # would poison the block if the cut failed)
    core = np.zeros((n, 1))
    core[select_shell_ao(mol, "1s")[0], 0] = 1.0
    mf = _FakeMF(np.hstack([core, c]), [2.0, 2.0, 1.0])
    b = valence_block(mf, chain, ncol=1)
    assert b.shape == (2, 2)
    assert np.allclose(b, b.T)
    # the symmetric column contributes w=2 x its member-summed outer
    # product; the asymmetric one w=1 x member-0 only; the core never
    want = (2.0 * 3 * np.outer(radial, radial)
            + 1.0 * np.outer([0.3, 0.1], [0.3, 0.1]))
    assert np.allclose(b, want, atol=1e-12)


# --------------------------------------------------------- the window

def _planted_config(mol, layer, shells, seed):
    """One synthetic series configuration: one CORE column (the
    valence cut's foil -- past ncol it must never enter the block)
    then per shell a planted radial on the declared bands, replicated
    on the m-set, doubly occupied (the block family reads the
    occupations, not the convergence)."""
    rng = np.random.default_rng(seed)
    n = mol.nao_nr()
    core = np.zeros((n, 1))
    core[select_shell_ao(mol, "1s")[0], 0] = 1.0
    cols, occs = [core], [2.0]
    for sh in shells:
        _chain, s_r = member_chains(mol, sh, layer[sh],
                                    mol.intor_symmetric("int1e_ovlp"))
        radial = rng.normal(size=len(_chain[0]))
        radial /= np.sqrt(radial @ s_r @ radial)
        c = np.zeros((n, 1))
        for m in range(len(_chain)):
            c[_chain[m], 0] = radial
        cols.append(c)
        occs.append(2.0)
    mf = _FakeMF(np.hstack(cols), occs)
    return mf


def test_series_window_recovers_the_head_radial_and_closes():
    """The @r form's exact anchor: one configuration over a multi-band
    family, k=1 -- the distilled radial IS the planted one on every
    member (capture 1, spectrum [1, ~0]), the window S-orthonormal
    and rotation-CLOSED (the one-particle certificate)."""
    from atomterm.operators import window_closure

    mol = uncontracted_carbon_mole()
    s = mol.intor_symmetric("int1e_ovlp")
    layer = {"2p": ("2p", "3p")}
    head = _planted_config(mol, layer, ["2p"], seed=7)
    cols, rep = series_window(mol, [head], s, ["2p"], layer,
                              {"2p": 1}, ncol=1)
    assert cols.shape == (mol.nao_nr(), 3)
    assert abs(rep["2p"]["spectrum"][0] - 1.0) < 1e-10
    assert rep["2p"]["spectrum"][1] < 1e-10
    assert abs(max(max(c) for c in rep["2p"]["capture"]) - 1.0) < 1e-8
    assert np.abs(cols.T @ s @ cols - np.eye(3)).max() < 1e-9
    lz, lperp = window_closure(mol, cols).values()
    assert lz < 1e-12 and lperp < 1e-12


def test_series_window_ladder_average_carries_the_second_radial():
    """The @s form's certificate: two configurations with DIFFERENT
    planted radials -- the average's spectrum carries a real second
    eigenvalue (k=2 certified, the discarded fraction small), and the
    k=2 window captures BOTH configurations' columns (the average's
    span holds both radials) while k=1 cannot."""
    mol = uncontracted_carbon_mole()
    s = mol.intor_symmetric("int1e_ovlp")
    layer = {"2p": ("2p", "3p")}
    series = [_planted_config(mol, layer, ["2p"], seed=1),
              _planted_config(mol, layer, ["2p"], seed=2)]
    cols2, rep2 = series_window(mol, series, s, ["2p"], layer,
                                {"2p": 2}, ncol=1)
    assert cols2.shape == (mol.nao_nr(), 6)
    assert rep2["2p"]["spectrum"][1] > 1e-2     # a real spread radial
    assert rep2["2p"]["discarded"] < 1e-10      # 2 bands, k=2: all kept
    for cfg in rep2["2p"]["capture"]:
        assert min(cfg) > 0.999
    _cols1, rep1 = series_window(mol, series, s, ["2p"], layer,
                                 {"2p": 1}, ncol=1)
    # the mean radial misses each configuration's radial direction:
    # some configuration's capture drops below the k=2 readout
    assert min(min(c) for c in rep1["2p"]["capture"]) < 0.999


def test_series_window_refuses_the_noise_second_radial():
    """The kept-signal certificate: a single configuration (rank-1
    block) asked for k=2 refuses -- the second radial would be
    numerical noise, the @r2 disease."""
    mol = uncontracted_carbon_mole()
    s = mol.intor_symmetric("int1e_ovlp")
    layer = {"2p": ("2p", "3p")}
    head = _planted_config(mol, layer, ["2p"], seed=7)
    with pytest.raises(ValueError, match="exceeds the radial signal"):
        series_window(mol, [head], s, ["2p"], layer, {"2p": 2}, ncol=1)


def test_series_window_refuses_k_beyond_the_chain_width():
    """The budget's floor: k cannot exceed the realization's radial
    width (the declared bands hold fewer directions)."""
    mol = uncontracted_carbon_mole()
    s = mol.intor_symmetric("int1e_ovlp")
    layer = {"2p": ("2p", "3p")}
    series = [_planted_config(mol, layer, ["2p"], seed=1),
              _planted_config(mol, layer, ["2p"], seed=2)]
    with pytest.raises(ValueError, match="exceeds its chain width"):
        series_window(mol, series, s, ["2p"], layer, {"2p": 3}, ncol=1)


# ---------------------------------------------------------- the series

def test_run_series_head_is_the_trusted_reference(tmp_path):
    """The head contract: the FIRST ladder entry is the case's own
    reference -- its mf returns untouched (no SCF, no chk written for
    it), the head's ENERGY is checked against the declaration (the
    head-must-be-the-case-reference convention as a checked
    invariant), and a spin-mismatched head is a refused declaration."""
    mol = uncontracted_carbon_mole()
    _m, head = None, _planted_config(mol, {"2p": ("2p", "3p")},
                                     ["2p"], seed=7)
    head.e_tot = -37.5                      # the declared reference's
    cfg = AtomConfig("C", "x", spin=2)
    series = (SeriesSpec("p2", ("1s", "2s"), ("2p",), 2, -37.5),)
    out = run_series(cfg, tmp_path, "C", series, head)
    assert out == {"p2": head}
    wrong_energy = (SeriesSpec("p2", ("1s", "2s"), ("2p",), 2, -38.0),)
    with pytest.raises(ValueError, match="not the series head"):
        run_series(cfg, tmp_path, "C", wrong_energy, head)
    bad = (SeriesSpec("p2", ("1s", "2s"), ("2p",), 0, -37.5),)
    with pytest.raises(ValueError, match="mismatched series"):
        run_series(cfg, tmp_path, "C", bad, head)


def test_run_series_threads_the_declaring_atom_wholesale(
        tmp_path, monkeypatch):
    """The members' configs inherit the declaring atom WHOLESALE
    (``dataclasses.replace``, only the spin varies): the ECP pair, the
    charge, and the artifact tag ride along.  The constructor call this
    replaced dropped them -- an ECP atom's members would build
    all-electron at 8x the electron count (the U14e wave's seam fix,
    2026-11-17); a stub pins the seam because the honest run is an
    SCF (campaign compute, not a unit pin)."""
    from atomterm import basis as _basis_mod
    from atomterm.system import AtomConfig

    head = type("Mf", (), {"e_tot": -50.5})()
    cfg = AtomConfig("U", "x", spin=4, ecp="ecp-block", tag="U_tag")
    series = (SeriesSpec("head", ("1s", "2s"), ("2p",), 4, -50.5),
              SeriesSpec("f4", ("1s", "2s"), ("2p",), 4, -50.0))
    seen = {}

    class _Stop(Exception):
        pass

    def _capture(cfg_x, _chk, ref=None):
        seen["cfg"], seen["ref"] = cfg_x, ref
        raise _Stop

    monkeypatch.setattr(_basis_mod, "run_rohf", _capture)
    with pytest.raises(_Stop):
        run_series(cfg, tmp_path, "U14e", series, head)
    mx = seen["cfg"]
    assert (mx.symbol, mx.basis, mx.spin) == ("U", "x", 4)
    assert mx.ecp == "ecp-block"            # the pair rides along
    assert mx.tag == "U_tag" and mx.charge == 0
    # the gate is the member's own declared basin energy
    assert seen["ref"].occ_e_ref == -50.0


def test_run_series_refuses_bad_self_seed_declarations(
        tmp_path, monkeypatch):
    """The self-seed's declaration contract refuses at the seam that
    consumes it: a seed on the head (it is never re-run -- inert is a
    wrong declaration), unknown seed forms, a half-declared pair
    (seed/open_counts declare each other), counts misaligned with the
    spans, and a count sum off the member's spin -- all BEFORE any
    run (the stub below must never fire)."""
    from atomterm import basis as _basis_mod

    head_mf = type("Mf", (), {"e_tot": -37.5})()
    cfg = AtomConfig("C", "x", spin=2)
    head = SeriesSpec("p2", ("1s", "2s"), ("2p",), 2, -37.5)

    class _Stop(Exception):
        pass

    def _stop(*_a, **_k):
        raise _Stop

    monkeypatch.setattr(_basis_mod, "run_rohf", _stop)
    with pytest.raises(ValueError, match="declares a seed"):
        run_series(cfg, tmp_path, "C", (head._replace(seed="self"),),
                   head_mf)
    cases = [
        ("declare each other",
         (head, SeriesSpec("s2", ("1s",), ("2s", "2p"), 2, -37.0,
                           seed="self"))),
        ("declare each other",
         (head, SeriesSpec("s2", ("1s",), ("2s", "2p"), 2, -37.0,
                           open_counts=(1, 1)))),
        ("unknown seed form",
         (head, SeriesSpec("s2", ("1s",), ("2s", "2p"), 2, -37.0,
                           seed="proj", open_counts=(1, 1)))),
        ("are the seed_basis thread's own",
         (head, SeriesSpec("s2", ("1s",), ("2s", "2p"), 2, -37.0,
                           seed="self", open_counts=(1, 1),
                           source_closed=("1s",)))),
        ("does not align",
         (head, SeriesSpec("s2", ("1s",), ("2s", "2p"), 2, -37.0,
                           seed="self", open_counts=(2,)))),
        ("spin is 2",
         (head, SeriesSpec("s2", ("1s",), ("2s", "2p"), 2, -37.0,
                           seed="self", open_counts=(2, 1)))),
    ]
    for match, series in cases:
        with pytest.raises(ValueError, match=match):
            run_series(cfg, tmp_path, "C", series, head_mf)


def test_run_series_threads_the_self_seed_from_the_head(
        tmp_path, monkeypatch):
    """seed="self" members MOM-run on a seed built from the declaring
    HEAD's own columns (mom_self_seed over the declared counts): the
    stub captures the (c_seed, occ_seed) pair and the ref -- NO
    seed_basis on it (the two-masters rule), the member's own gate,
    the chk at the series address."""
    from atomterm import basis as _basis_mod

    mol = uncontracted_carbon_mole()
    eye = np.eye(mol.nao_nr())
    occ = np.zeros(mol.nao_nr())
    occ[[0, 1]] = 2.0                       # the head's closed block
    occ[3] = 1.0                            # one open 2px
    head_mf = type("Mf", (), {"e_tot": -37.5, "mol": mol,
                              "mo_coeff": eye, "mo_occ": occ,
                              "mo_energy": np.linspace(-5, 5,
                                                       mol.nao_nr())})()
    cfg = AtomConfig("C", "x", spin=2)
    series = (SeriesSpec("p2", ("1s", "2s"), ("2p",), 2, -37.5),
              SeriesSpec("p2b", ("1s",), ("2p",), 2, -37.0,
                         seed="self", open_counts=(2,)))
    seen = {}

    class _Stop(Exception):
        pass

    def _capture(cfg_x, chk, ref=None, mom_seed=None):
        seen["cfg"], seen["chk"], seen["ref"], seen["seed"] = (
            cfg_x, chk, ref, mom_seed)
        raise _Stop

    monkeypatch.setattr(_basis_mod, "run_rohf", _capture)
    with pytest.raises(_Stop):
        run_series(cfg, tmp_path, "C", series, head_mf)
    c_seed, occ_seed = seen["seed"]
    assert list(occ_seed) == [2.0, 2.0, 1.0, 1.0]   # pairs first
    for j, want in enumerate((0, 1, 3, 4)):         # closed + the pick
        assert np.array_equal(c_seed[:, j], eye[:, want])
    assert seen["ref"].seed_basis is None   # one seed source per run
    assert seen["ref"].occ_e_ref == -37.0   # the member's own gate
    assert str(seen["chk"]) == (str(tmp_path / "C" / "rungs" / "series"
                                    / "p2b.chk"))
