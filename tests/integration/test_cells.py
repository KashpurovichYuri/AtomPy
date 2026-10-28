"""The cell layer: key addressing, the store contract, and the CI
ordering bridge.

The bridge is the load-bearing piece: PySCF FCI vectors and
``detbasis`` enumerations order determinants differently for norb >= 4
with >= 2 electrons in a block (they coincide on carbon's norb=3
window -- why the difference went unnoticed). The certificate test
verifies the mapping through the eigen-equation against an
independently built Slater-Condon H: generic integrals admit no
permutation automorphism, so a wrong mapping cannot satisfy it.
"""

import numpy as np
import pytest
import re

from atomterm.cells import (
    CellData,
    CellKey,
    cell_chk,
    ci_in_detbasis_order,
    diagnostics,
    grid_keys,
    load_cell,
    planned,
    pyscf_ci_permutation,
    spine_keys,
    window_grid,
)
from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.windows import span_overlap

# mixed tiers in one module (1:1 with src): the key/bridge/store tests
# are fast; the carbon run pins are slow -- marked per test, so the
# module-level mark would leak SCF fixtures into `pytest -m fast`

fast = pytest.mark.fast
slow = pytest.mark.slow


# --------------------------------------------------------------- the key

@fast
def test_address_round_trip():
    keys = [
        CellKey("C", "2p", (1, 1)),
        CellKey("C", "2p", (1, 1), fold="raw"),
        CellKey("c", "2p", (2, 0), method="sa", mask=None, substrate="pur"),
        CellKey("Ti", "3d4s", (3, 1), method="sa", start="ao",
                mask="radial", substrate="pur", fold="avg"),
        CellKey("Ti", "3d4s", (2, 2), method="casscf", start="ao",
                mask="virtmix"),
        CellKey("Sc", "3d4s", (1, 1), method="casscf", mask="frozen",
                fold="frac"),
        CellKey("U", "3f2d", (4, 2), method="sa", start="ao"),
    ]
    for k in keys:
        assert CellKey.from_address(k.address) == k
        assert k.address.count("_") == 7
    assert keys[0].address == "c_2p_1-1_casci_ao_-_pol_avg"
    assert keys[3].address == "ti_3d4s_3-1_sa_ao_radial_pur_avg"


@fast
def test_key_validation_and_canonicalization():
    # the pur-fold degeneracy rule: raw on pur canonicalizes to avg
    k = CellKey("C", "2p", (1, 1), substrate="pur", fold="raw")
    assert k.fold == "avg"
    # frac is a different REFERENCE -- meaningless on the purified span
    with pytest.raises(ValueError, match="frac"):
        CellKey("C", "2p", (1, 1), substrate="pur", fold="frac")
    # a mask needs an optimizer method
    with pytest.raises(ValueError, match="mask"):
        CellKey("C", "2p", (1, 1), method="casci", mask="radial")
    # the radial mask rotates member-labeled AO slots: ao-start only
    with pytest.raises(ValueError, match="radial mask"):
        CellKey("C", "2p", (1, 1), method="casscf", start="mo",
                mask="radial")
    with pytest.raises(ValueError, match="radial mask"):
        CellKey("C", "2pno1", (1, 1), method="casscf", start="no",
                mask="radial")
    # mo x pur: VALID since the window-projection builder (T1a) -- the
    # canonical columns are re-orthogonalized on the pur core span; the
    # pur-fold canonicalization applies to it like every pur key
    k = CellKey("Ti", "3d4s", (3, 1), method="casscf", start="mo",
                substrate="pur", fold="raw")
    assert k.fold == "avg"
    assert k.address == "ti_3d4s_3-1_casscf_mo_-_pur_avg"
    for bad in [dict(method="fci"), dict(start="aos"), dict(mask="relax"),
                dict(substrate="polarized"), dict(fold="averaged")]:
        with pytest.raises(ValueError, match="not in"):
            CellKey("C", "2p", (1, 1), **bad)
    with pytest.raises(ValueError, match="sector"):
        CellKey("C", "2p", (0, 1))
    with pytest.raises(ValueError, match="shell label"):
        CellKey("C", "3d_4s", (1, 1))
    # the window label parses to shell selectors
    assert CellKey("Ti", "3d4s", (3, 1)).shell_patterns == ["3d", "4s"]
    assert CellKey("C", "2p", (1, 1)).shell_patterns == ["2p"]


@fast
def test_no_window_grammar():
    """NO windows: '<source>no<k>' -- the group prefix is part of the
    window identity (differently-sized NO cells get different
    addresses); start and label must agree."""
    k = CellKey("Ti", "3d4sno2", (3, 1), start="no")
    assert k.no_source == ("3d4s", 2)
    assert k.shell_patterns == ["3d", "4s"]     # the SOURCE shells
    assert CellKey.from_address(k.address) == k
    assert k.address == "ti_3d4sno2_3-1_casci_no_-_pol_avg"
    with pytest.raises(ValueError, match="agree"):
        CellKey("Ti", "3d4s", (3, 1), start="no")
    with pytest.raises(ValueError, match="agree"):
        CellKey("Ti", "3d4sno2", (3, 1))        # ao start, no-tail label


@fast
def test_spine_keys_contract():
    """The spine = README · Methods verbatim: the three folds; frozen x
    purified under both objectives (the certified class); one free row
    per substrate (the basin detectors)."""
    keys = spine_keys("Ti", "3d4s", (3, 1))
    assert [k.address for k in keys] == [
        "ti_3d4s_3-1_casci_ao_-_pol_raw",
        "ti_3d4s_3-1_casci_ao_-_pol_avg",
        "ti_3d4s_3-1_casci_ao_-_pol_frac",
        "ti_3d4s_3-1_casscf_ao_frozen_pur_avg",
        "ti_3d4s_3-1_sa_ao_frozen_pur_avg",
        "ti_3d4s_3-1_casscf_ao_-_pol_avg",
        "ti_3d4s_3-1_sa_ao_-_pur_avg",
    ]


@fast
def test_window_grid_is_the_one_declaration():
    """The declared key-set derived from CaseSpec (T1b): the full grid in
    render order, the spine a subset of it, the named groups filters.
    mo x pur IS in the grid (the projection builder, T1a); radial stays
    ao-only; the NO start addresses the full-span group prefix with
    k = the window's shell count."""
    from atomterm.cases import CAMPAIGN

    spec = CAMPAIGN["Ti"]
    grid = grid_keys(spec)
    addrs = [k.address for k in grid]
    # the manifest count: 29 established + the 4 mo x pur cells
    assert len(grid) == 33
    assert len(set(addrs)) == 33
    assert addrs[:4] == [           # the CASCI fold block, pur third
        "ti_3d4s_3-1_casci_ao_-_pol_raw",
        "ti_3d4s_3-1_casci_ao_-_pol_avg",
        "ti_3d4s_3-1_casci_ao_-_pur_avg",
        "ti_3d4s_3-1_casci_ao_-_pol_frac",
    ]
    mopur = [a for a in addrs if a.split("_")[4] == "mo"
             and a.split("_")[6] == "pur"]
    assert mopur == [                   # method-major, like every block
        "ti_3d4s_3-1_casscf_mo_-_pur_avg",
        "ti_3d4s_3-1_casscf_mo_frozen_pur_avg",
        "ti_3d4s_3-1_sa_mo_-_pur_avg",
        "ti_3d4s_3-1_sa_mo_frozen_pur_avg",
    ]
    # no radial x mo/no keys ever (void at the key)
    assert not [a for a in addrs if "_radial_" in a and "_ao_" not in a]
    # the NO windows: full-span group prefixes, k = shell count
    assert "ti_3d4sno2_3-1_casscf_no_-_pol_avg" in addrs
    # the spine is a subset view; diagnostics = the virtmix filter
    spine = {k.address for k in spine_keys("Ti", "3d4s", (3, 1))}
    assert spine <= set(addrs)
    assert [k.mask for k in diagnostics(spec)] == ["virtmix"]
    # planned: the declared windows' full grids, same shape as the primary
    pl = planned(spec)
    assert len(pl) == 2 * len(grid)
    assert {k.window for k in pl} == {
        "3d4d", "3d4dno2", "3d4s4p4d", "3d4s4p4dno4"}
    assert all(k.address not in addrs for k in pl)   # disjoint windows


# ------------------------------------------------------------- the bridge

SECTORS = [(3, 1, 1), (4, 2, 1), (5, 3, 1), (6, 3, 2), (6, 4, 2)]


@fast
def test_ci_permutation_is_a_permutation():
    from math import comb

    for norb, na, nb in SECTORS:
        perm = pyscf_ci_permutation(norb, na, nb)
        assert len(perm) == comb(norb, na) * comb(norb, nb)
        assert np.array_equal(np.sort(perm), np.arange(len(perm)))


@fast
def test_orders_coincide_only_in_trivial_sectors():
    """The empirical fact that kept the bug invisible: carbon's norb=3
    window forces pyscf's integer-ascending string order to coincide
    with detbasis' lexicographic order; norb >= 4 breaks it."""
    assert np.array_equal(pyscf_ci_permutation(3, 1, 1), np.arange(9))
    assert not np.array_equal(pyscf_ci_permutation(4, 2, 1), np.arange(24))


def _random_eri(norb, rng):
    """A generic chemist-notation eri with the full 8-fold symmetry
    (average over the orbit under the pair transpositions) -- PySCF's
    contraction assumes it."""
    a = rng.standard_normal((norb,) * 4)
    seen, fr = [], [a]
    while fr:
        x = fr.pop()
        if any(np.array_equal(x, s) for s in seen):
            continue
        seen.append(x)
        for p in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
            fr.append(np.ascontiguousarray(x.transpose(p)))
    eri = sum(seen) / len(seen)
    for p in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
        assert np.allclose(eri, eri.transpose(p))
    return eri


@fast
@pytest.mark.parametrize("norb,na,nb", SECTORS)
def test_bridge_certified_by_eigen_equation(norb, na, nb):
    """THE bridge certificate: a PySCF-solved CI vector, mapped into
    detbasis order, satisfies the eigen-equation of the independently
    built Slater-Condon Hamiltonian. No ordering assumption anywhere --
    generic integrals admit no permutation automorphism."""
    from pyscf import fci
    from atomterm.operators import spin_block_one_body, spin_block_two_body

    rng = np.random.default_rng(20261027 + 10 * norb + na)
    h1 = rng.standard_normal((norb, norb))
    h1 = (h1 + h1.T) / 2
    eri = _random_eri(norb, rng)

    e_fci, ci = fci.direct_spin1.kernel(h1, eri, norb, (na, nb),
                                        conv_tol=1e-13)
    dets = determinant_basis(norb, na, nb)
    # the owners' spin-orbital expansion; detbasis eats the PHYSICIST
    # arrangement of the same eri (fact 3), pyscf the chemist one
    ham = slater_condon_matrix(
        dets, one_body=spin_block_one_body(h1),
        two_body=spin_block_two_body(eri.transpose(0, 2, 1, 3)))

    w = ci_in_detbasis_order(ci, norb, na, nb)
    e_mine = float(np.real(w.conj() @ ham @ w))
    assert abs(e_mine - e_fci) < 1e-9
    resid = ham @ w - e_fci * w
    assert np.linalg.norm(resid) / np.linalg.norm(ham) < 1e-10


@fast
def test_configuration_summary_uses_the_bridge(tmp_path):
    """The consumer fix: a pyscf-ordered unit amplitude must label the
    determinant at ITS string address, not dets[k]."""
    from pyscf import gto
    from pyscf.fci import cistring
    from atomterm.analysis import configuration_summary

    mol = gto.M(atom="C 0 0 0", basis="cc-pvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    # a mixed-l window (2s + 2p) so column identity matters
    cols = [1, 2, 3, 4]
    c_win = np.eye(mol.nao_nr())[:, cols]
    na, nb = 2, 1
    # pyscf address 5 of the (2,1)/(4) sector -> its determinant
    a_strs = cistring.make_strings(range(4), na)
    b_strs = cistring.make_strings(range(4), nb)
    ia, ib = divmod(5, len(b_strs))
    sa, sb = a_strs[ia], b_strs[ib]
    occ_a = [i for i in range(4) if sa >> i & 1]
    occ_b = [i for i in range(4) if sb >> i & 1]
    # expected label: the l-mass of exactly those columns
    from atomterm.drivers import shell_weights
    names, wl = shell_weights(mol, c_win, s)
    occ_cols = sorted(occ_a + occ_b)
    mass = {n: sum(wl[c, k] for c in occ_cols)
            for k, n in enumerate(names)}
    expected = "".join(f"{n}^{mass[n]:.2f}" for n in names if mass[n] > 0.01
                       ) + " (|c|=1.00)"
    v = np.zeros(len(a_strs) * len(b_strs))
    v[5] = 1.0
    [label] = configuration_summary(mol, s, [v], c_win, na, nb, top=1)
    assert label == expected


# -------------------------------------------------------- the store layer

@fast
def test_load_cell_synthetic_roundtrip(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1))  # any atom label; synthetic store
    chk = cell_chk(key, tmp_path)
    chk.parent.mkdir(parents=True)
    rng = np.random.default_rng(7)
    mo = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
    energies = np.array([-1.0, -0.5, -0.25])
    ci = [rng.standard_normal(9), rng.standard_normal(9)]
    for name, val in (("energies", energies),
                      ("ci_stack", np.asarray(ci)),
                      ("core_energy", 5.25), ("mo_coeff", mo),
                      ("ncas", 3)):
        lib.chkfile.save(str(chk), f"atomterm/{key.address}/{name}", val)

    data = load_cell(key, tmp_path)
    assert isinstance(data, CellData)
    assert data.key == key
    assert np.array_equal(data.mo_coeff, mo)
    assert np.array_equal(data.energies, energies)
    assert data.core_energy == 5.25
    assert data.ncas == 3
    assert len(data.ci) == 2
    assert np.allclose(data.ci[1], ci[1])


@fast
def test_load_cell_prefers_sa_roots_and_reports_missing(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1), method="sa")
    chk = cell_chk(key, tmp_path)
    chk.parent.mkdir(parents=True)
    for name in (key.address, f"{key.address}_roots"):
        lib.chkfile.save(str(chk), f"atomterm/{name}/energies",
                         np.array([-2.0 if "roots" in name else -1.0]))
        lib.chkfile.save(str(chk), f"atomterm/{name}/ci_stack",
                         np.zeros((1, 9)))
        lib.chkfile.save(str(chk), f"atomterm/{name}/core_energy", 1.0)
    lib.chkfile.save(str(chk), f"atomterm/{key.address}/mo_coeff",
                     np.eye(6))
    assert load_cell(key, tmp_path).energies[0] == -2.0

    with pytest.raises(FileNotFoundError, match="run it first"):
        load_cell(CellKey("Xe", "5p", (2, 0)), tmp_path)


# ------------------------------------------------- the run pins (carbon)

def _same_facts(a, b, tol=1e-12):
    """Two facts objects describe the same run: energies, meters, labels."""
    assert a.tag == b.tag
    assert np.allclose(a.energies, b.energies, atol=tol), (a.energies, b.energies)
    assert np.allclose(a.energies_pyscf, b.energies_pyscf, atol=tol)
    for f in ("comm_norm_h_lz", "max_de_vs_pyscf",
              "ortho_err", "n_dev", "win_lz", "win_lperp",
              "ground_l2", "ground_lz2"):
        assert abs(getattr(a, f) - getattr(b, f)) <= tol, (f, getattr(a, f),
                                                           getattr(b, f))
    assert a.ground_term == b.ground_term
    assert a.terms == b.terms
    # configuration labels modulo the |c| amplitudes: inside a degenerate
    # cluster the Davidson gauge is chaotic even on pinned inputs (the
    # run_casci straggle note) -- the l-mass content is the invariant
    strip = lambda cs: [re.sub(r" \(\|c\|=[0-9.]+\)", "", c) for c in cs]
    assert strip(a.configurations) == strip(b.configurations)
    if a.core_anisotropy is not None or b.core_anisotropy is not None:
        assert a.core_anisotropy == pytest.approx(b.core_anisotropy, abs=1e-12)
    assert a.analysis.comm_norm_h_l2 == pytest.approx(
        b.analysis.comm_norm_h_l2, abs=tol)


@slow
def test_run_cell_matches_direct_path(carbon, tmp_path):
    """run_cell(casci cell) == the direct run_casci + stage_facts path --
    the same construction the rungs registry uses, byte-comparable on a
    pinned reference (fact 10: CASCI on pinned MOs is deterministic)."""
    from atomterm.casscf import run_casci
    from atomterm.stages import stage_facts
    from atomterm.system import ncore_of
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import run_cell, cell_payload

    _, mol, mf, act_idx = carbon
    na = nb = 1
    ncore = ncore_of(mol.nelectron, 2)
    mo = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
    direct = run_casci(mf, len(act_idx), (na, nb), mo,
                       str(tmp_path / "direct.chk"), nroots=9, name="casci_avg")
    key = CellKey("C", "2p", (na, nb))       # casci/ao/-/pol/avg
    facts_direct = stage_facts(
        mol, mf, direct, direct.mo_coeff[:, ncore:ncore + len(act_idx)],
        ncore, na, nb, tag=key.address, spherical_core=True,
        anisotropy_ref=(mf, None))

    facts = run_cell(key, mf, workdir=tmp_path)
    assert facts.tag == key.address == "c_2p_1-1_casci_ao_-_pol_avg"
    _same_facts(facts, facts_direct)
    payload = cell_payload(key, facts)
    assert payload["cell"] == key.address

    # the cache hit re-analyzes the stored run -- same facts again
    _same_facts(run_cell(key, mf, workdir=tmp_path), facts)


@slow
def test_analyze_cell_matches_insession(carbon, tmp_path):
    """THE store certificate: analyze_cell (load_cell + facts) reproduces
    the in-session run exactly -- energies, meters, AND the stored-CI
    configurations (the bridge applied to real stored data)."""
    from atomterm.cells import analyze_cell, run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), fold="raw")   # the detector cell
    insession = run_cell(key, mf, workdir=tmp_path, fresh=True)
    reloaded = analyze_cell(key, mol, mf, workdir=tmp_path)
    _same_facts(reloaded, insession)
    assert reloaded.core_anisotropy == insession.core_anisotropy


@slow
def test_posthoc_properties_from_stored_run(carbon, tmp_path):
    """Post-hoc physics from the STORED run alone (load_cell + the
    bridge): the ground root's <L^2>/<S^2> from the stored PySCF CI and
    the energy from the stored window reproduce the payload. Gauge-safe
    observables only -- <Lz^2> is member-gauge-dependent inside a
    degenerate cluster (fact 16) and is NOT compared."""
    from atomterm.cells import cell_payload, load_cell, run_cell
    from atomterm.operators import build_operators
    from atomterm.system import ncore_of

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1))
    facts = run_cell(key, mf, workdir=tmp_path)
    data = load_cell(key, tmp_path)
    ncore = ncore_of(mol.nelectron, 2)
    c_act = data.mo_coeff[:, ncore:ncore + 3]
    ops = build_operators(mol, mf, c_act, ncore, 1, 1)
    w = ci_in_detbasis_order(data.ci[0], 3, 1, 1)
    l2 = float(np.real(w @ ops["L2"].real @ w))
    s2 = float(np.real(w @ ops["S2"].real @ w))
    e = float(np.real(w @ ops["H"] @ w)) + data.core_energy
    p = cell_payload(key, facts)
    assert abs(l2 - p["ground_l2"]) < 1e-9      # 2.0: gauge-safe in the 3P cluster
    assert abs(s2 - 2.0) < 1e-9                 # S=1 across the whole sector
    assert abs(e - p["e_ground"]) < 1e-9


@slow
def test_run_cell_no_start_same_subspace(carbon, tmp_path):
    """The NO start, end to end: on carbon's 2p source the term-averaged
    density's occupation groups form ONE 3-column m-set, so the no1
    window IS the source span in the NO gauge -- the CASCI must be
    basis-invariant (exact diagonalization inside one subspace; fact:
    CASCI is start-irrelevant, window-dependent)."""
    from atomterm.cells import load_cell, run_cell

    _, mol, mf, _ = carbon
    ao = run_cell(CellKey("C", "2p", (1, 1)), mf, workdir=tmp_path)
    no = run_cell(CellKey("C", "2pno1", (1, 1), start="no"), mf,
                  workdir=tmp_path)
    assert load_cell(CellKey("C", "2pno1", (1, 1), start="no"),
                     tmp_path).ncas == 3
    assert np.allclose(ao.energies, no.energies, atol=1e-9)
    assert no.win_lz == pytest.approx(0.0, abs=1e-12)   # m-set: closed span
    assert no.win_lperp == pytest.approx(0.0, abs=1e-12)


@slow
def test_sa_run_completes_the_roots_store(carbon, tmp_path):
    """store_name= pins the WRITER side of the SA store (the aa21071
    bug): the SA run returns the per-ROOTS pass, and the cell-layer
    completion must land under ``<address>_roots`` -- completing
    ``<address>`` instead pairs the ensemble's energies with the
    roots' CI, and the loader's preferred store then fails loud
    (missing ci_stack)."""
    from pyscf import lib

    from atomterm.cells import run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="sa")
    facts = run_cell(key, mf, workdir=tmp_path)
    data = load_cell(key, tmp_path)          # prefers <address>_roots
    assert np.allclose(data.energies, facts.energies_pyscf, atol=1e-10)
    chk = cell_chk(key, tmp_path)
    ci = lib.chkfile.load(str(chk), f"atomterm/{key.address}_roots/ci_stack")
    assert ci is not None and len(ci) == len(facts.energies_pyscf)


@slow
def test_frozen_mask_applies_to_the_mo_start(carbon, tmp_path):
    """Masks are START-INDEPENDENT (the 141c1b8 wiring): a frozen-window
    MO-start cell must converge WITH its start columns -- if the mask
    resolution slipped back into the ao branch, this cell would run
    free and the window would rotate away (carbon's radial relaxation,
    the theta* gain)."""
    from atomterm.cells import run_cell
    from atomterm.system import ncore_of

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="casscf", start="mo",
                  mask="frozen")
    run_cell(key, mf, workdir=tmp_path)
    ncore = ncore_of(mol.nelectron, 2)
    start = mf.mo_coeff[:, ncore:ncore + 3]
    conv = load_cell(key, tmp_path).mo_coeff[:, ncore:ncore + 3]
    s = mf.get_ovlp()
    for i in range(3):
        # |<start_i|S|conv_i>| == 1: the same column up to phase
        assert abs(start[:, i] @ (s @ conv[:, i])) == pytest.approx(
            1.0, abs=1e-8), f"window column {i} moved under a frozen mask"


@slow
def test_run_cell_mo_start_on_pur_substrate(carbon, tmp_path):
    """The mo x pur cell end to end (T1a): canonical window+virtual
    columns re-orthogonalized on the pur core through build_full_mo's
    one treatment. On carbon the s-core x p-window are exactly
    S-orthogonal (different l) -- the projection is the no-op case; the
    defect-scale case is Ti's, measured when the gated batch runs.
    Machinery pins: the S-orthonormal substrate, the pur core block,
    CASCI internal consistency."""
    from atomterm.cells import load_cell, run_cell
    from atomterm.core import purified_core_orbitals
    from atomterm.system import ncore_of

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), start="mo", substrate="pur")
    facts = run_cell(key, mf, workdir=tmp_path)
    assert facts.ortho_err < 1e-10
    assert facts.max_de_vs_pyscf < 1e-9
    ncore = ncore_of(mol.nelectron, 2)
    c_pur, _ = purified_core_orbitals(mol, mf, ncore)
    data = load_cell(key, tmp_path)
    s = mf.get_ovlp()
    err = np.abs(data.mo_coeff.conj().T @ s @ data.mo_coeff
                 - np.eye(mol.nao_nr())).max()
    assert err < 1e-10
    # substrate = the core SPAN (PySCF canonicalizes the inactive block
    # WITHIN its span -- the fact-6 mutation, span-preserving; it hits
    # the ao-pur cells identically): the stored core spans the pur core
    assert span_overlap(data.mo_coeff[:, :ncore], c_pur, s) == \
        pytest.approx(1.0, abs=1e-10)


@fast
def test_select_no_window_selects_by_group_number():
    """k counts GROUPS (window identity), but no_group_windows returns
    a COMPACT list -- prefixes too small to open the sector are
    skipped, so labels can be non-contiguous. Selection is by label,
    never by index (the Ti {3d,4s}/(3,1) lattice failure class: the
    1-column s prefix skipped, no2 = the full 6-column span)."""
    from atomterm.cells import _select_no_window
    from atomterm.windows import Window

    wins = [Window("no2", None, np.eye(6, 6), True),
            Window("no4", None, np.eye(6, 6), True)]   # no1/no3 skipped
    assert _select_no_window(wins, 2).label == "no2"
    assert _select_no_window(wins, 4).c_win.shape[1] == 6
    with pytest.raises(ValueError, match="usable group prefixes"):
        _select_no_window(wins, 1)
    with pytest.raises(ValueError, match="no3 is not among them"):
        _select_no_window(wins, 3)
