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
    diagnostics,
    grid_keys,
    lane_keys,
    lattice_keys,
    load_cell,
    load_lane_cell,
    planned,
    run_lane_cell,
    spine_keys,
    void_reason,
    window_grid,
)
from atomterm.cases import CAMPAIGN, scenarios_for
from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.operators import ci_in_detbasis_order, pyscf_ci_permutation
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
        CellKey("c", "2p", (2, 0), method="sa-casscf", mask=None, substrate="pur"),
        CellKey("Ti", "3d4s", (3, 1), method="sa-casscf", start="ao",
                mask="radial", substrate="pur", fold="avg"),
        CellKey("Ti", "3d4s", (2, 2), method="casscf", start="ao",
                mask="virtmix"),
        CellKey("Sc", "3d4s", (1, 1), method="casscf", mask="frozen",
                fold="frac"),
        CellKey("U", "3f2d", (4, 2), method="sa-casscf", start="ao"),
    ]
    for k in keys:
        assert CellKey.from_address(k.address) == k
        assert k.address.count("_") == 7
    assert keys[0].address == "c_2p_1-1_casci_ao_-_pol_avg"
    assert keys[3].address == "ti_3d4s_3-1_sa-casscf_ao_radial_pur_avg"


@fast
def test_docstring_example_address_parses():
    """The module docstring's decoded example is a REAL address: it
    round-trips through ``from_address`` (the illustration it replaced
    was invalid -- a plain window under a no start, and a 'free' mask
    token the grammar does not define)."""
    import atomterm.cells as cells_mod

    m = re.search(r"^    (ti_\S+)$", cells_mod.__doc__, re.M)
    assert m, "the decoded example address is missing from the docstring"
    assert CellKey.from_address(m.group(1)).address == m.group(1)


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
def test_lane_window_grammar():
    """A lane window carries its fold as a tag: the tag is part of the
    address (never collides with a free-route label), stays out of the
    span selectors, and resolves through the fold vocabulary; NO
    windows take no tag (their fold is the source's)."""
    k = CellKey("ti", "3d3p4d@semi", (5, 3), substrate="pur")
    assert CellKey.from_address(k.address) == k
    assert k.address == "ti_3d3p4d@semi_5-3_casci_ao_-_pur_avg"
    assert k.shell_patterns == ["3d", "3p", "4d"]
    assert k.fold_tag == "semi"
    assert k.fold_shells() == ("1s", "2s", "2p", "3s", "4s")
    with pytest.raises(ValueError, match="fold tag"):
        CellKey("Ti", "3d4sno2@semi", (3, 1), start="no")
    with pytest.raises(ValueError, match="fold tag"):
        CellKey("Ti", "3d4s", (3, 1)).fold_shells()
    with pytest.raises(ValueError, match="neither"):
        CellKey("Ti", "3d3p@nope", (5, 3)).fold_shells()


@fast
def test_lane_keys_derive_from_the_declaration():
    """The lane legs and their SECTORS derive from the declaration:
    activating a shell moves its electrons into the window, so the
    sector follows the fold (Ti semi grows (3,1) -> (5,3), the
    M-shell folds to (6,4); the Ar legs keep the case sector); Sc's
    five-leg family (d1: one electron fewer per M-shell leg)."""
    ti = {k.window: k.sector for k in lane_keys(CAMPAIGN["Ti"])}
    assert ti["3d4s@ar"] == (3, 1)
    assert ti["3d4d@ar"] == (3, 1)
    assert ti["3d3p@semi"] == (5, 3)
    assert ti["3d3p4d@semi"] == (5, 3)
    assert ti["3s3p3d@1s2s2p4s"] == (6, 4)
    assert ti["3d3p4s@1s2s2p3s"] == (6, 4)
    assert [(k.window, k.sector) for k in lane_keys(CAMPAIGN["Sc"])] == [
        ("3d4s@ar", (2, 1)), ("3d3p@semi", (4, 3)),
        ("3d3p4d@semi", (4, 3)), ("3d3p4s@1s2s2p3s", (5, 4)),
        ("3s3p3d@1s2s2p4s", (5, 4))]
    assert lane_keys(CAMPAIGN["U"]) == []


@fast
def test_void_reason_gates_scenariod_cells():
    """The run gate: the anchor-first scenario voids EXACTLY the
    planned windows' optimizer cells (both label forms); the CASCI
    anchors and the primary grid run.  A void cell stays addressable
    (the key constructs; only the run is gated)."""
    spec = CAMPAIGN["Ti"]
    sc = scenarios_for(spec)
    for k in planned(spec):
        assert (void_reason(k, sc) is not None) == (k.method != "casci"), \
            k.address
    assert all(void_reason(k, sc) is None for k in grid_keys(spec))
    reason = void_reason("ti_3d4d_3-1_casscf_ao_-_pol_avg", sc)
    assert reason.startswith("anchor-first") and "2026-11-02" in reason


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
        "ti_3d4s_3-1_sa-casscf_ao_frozen_pur_avg",
        "ti_3d4s_3-1_casscf_ao_-_pol_avg",
        "ti_3d4s_3-1_sa-casscf_ao_-_pur_avg",
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
        "ti_3d4s_3-1_sa-casscf_mo_-_pur_avg",
        "ti_3d4s_3-1_sa-casscf_mo_frozen_pur_avg",
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


@fast
def test_lattice_keys_derive_from_the_declaration():
    """The lattice derives from CaseSpec.lattice_windows at the
    declared lattice_sector: two constructions per span (AO -- exactly
    closed; MO -- the m-delocalized control); undeclared atoms derive
    nothing; the sector is a declaration, not a constant (Sc's valence
    series differs from Ti's free-ion series)."""
    from atomterm.cases import CAMPAIGN

    assert lattice_keys(CAMPAIGN["C"]) == []
    assert lattice_keys(CAMPAIGN["Si"]) == []
    assert [k.address for k in lattice_keys(CAMPAIGN["Ti"])] == [
        "ti_3d_2-0_casci_ao_-_pol_avg",
        "ti_3d_2-0_casci_mo_-_pol_avg",
        "ti_3d4p_2-0_casci_ao_-_pol_avg",
        "ti_3d4p_2-0_casci_mo_-_pol_avg",
        "ti_3d4p4d_2-0_casci_ao_-_pol_avg",
        "ti_3d4p4d_2-0_casci_mo_-_pol_avg",
    ]
    # the valence series at its own sector: CAS(3,·), Ms = 1/2
    first = lattice_keys(CAMPAIGN["Sc"])[0]
    assert first.address == "sc_3d_2-1_casci_ao_-_pol_avg"
    assert [k.sector for k in lattice_keys(CAMPAIGN["Sc"])] == [(2, 1)] * 16


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
    # arrangement of the same eri; pyscf takes the chemist one
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
    from atomterm.windows import shell_weights
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
def test_load_cell_prefers_sa_casscf_roots_and_reports_missing(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1), method="sa-casscf")
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
    pinned reference (CASCI on pinned MOs is deterministic)."""
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
    degenerate cluster and is NOT compared."""
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
def test_sa_casscf_run_completes_the_roots_store(carbon, tmp_path):
    """store_name= pins the WRITER side of the SA store (the aa21071
    bug): the SA run returns the per-ROOTS pass, and the cell-layer
    completion must land under ``<address>_roots`` -- completing
    ``<address>`` instead pairs the ensemble's energies with the
    roots' CI, and the loader's preferred store then fails loud
    (missing ci_stack)."""
    from pyscf import lib

    from atomterm.cells import run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="sa-casscf")
    facts = run_cell(key, mf, workdir=tmp_path)
    data = load_cell(key, tmp_path)          # prefers <address>_roots
    assert np.allclose(data.energies, facts.energies_pyscf, atol=1e-10)
    chk = cell_chk(key, tmp_path)
    ci = lib.chkfile.load(str(chk), f"atomterm/{key.address}_roots/ci_stack")
    assert ci is not None and len(ci) == len(facts.energies_pyscf)


@slow
def test_frozen_mask_applies_to_the_mo_start(carbon, tmp_path):
    """Masks are START-INDEPENDENT: a frozen-window
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
    # WITHIN its span after the kernel (span-preserving; it hits
    # the ao-pur cells identically): the stored core spans the pur core
    assert span_overlap(data.mo_coeff[:, :ncore], c_pur, s) == \
        pytest.approx(1.0, abs=1e-10)


@fast
def test_pick_no_window_selects_by_group_number():
    """k counts GROUPS (window identity), but no_group_windows returns
    a COMPACT list -- prefixes too small to open the sector are
    skipped, so labels can be non-contiguous. Selection is by label,
    never by index (the Ti {3d,4s}/(3,1) lattice failure class: the
    1-column s prefix skipped, no2 = the full 6-column span)."""
    from atomterm.cells import _pick_no_window
    from atomterm.windows import Window

    wins = [Window("no2", None, np.eye(6, 6), True),
            Window("no4", None, np.eye(6, 6), True)]   # no1/no3 skipped
    assert _pick_no_window(wins, 2).label == "no2"
    assert _pick_no_window(wins, 4).c_win.shape[1] == 6
    with pytest.raises(ValueError, match="usable group prefixes"):
        _pick_no_window(wins, 1)
    with pytest.raises(ValueError, match="no3 is not among them"):
        _pick_no_window(wins, 3)


@slow
def test_run_lane_cell_carbon_roundtrip(carbon_d2h, tmp_path):
    """The lane runner on the frozen D2h carbon twin: the declared leg
    hosts its term exactly (all 2L+1 components read L(L+1), S(S+1),
    degenerate), the payload round-trips the store, and the runner
    refuses non-lane keys, wrong sectors, and absent fold shells."""
    cfg, mol, mf, _ = carbon_d2h
    key = lane_keys(CAMPAIGN["C"])[0]
    assert (key.window, key.sector) == ("2p@1s2s", (2, 0))
    p = run_lane_cell(key, mol, mf, CAMPAIGN["C"].term, workdir=tmp_path)
    assert p["cell"] == key.address and p["term"] == [1, "g"]
    assert len(p["energies"]) == 3 and p["ptp"] < 1e-9
    for v in p["l2"]:
        assert abs(v - 2.0) < 1e-8          # 3P: L(L+1) = 2
    for v in p["s2"]:
        assert abs(v - 2.0) < 1e-8          # S(S+1) = 2
    # the store round-trips (and a fresh run is a cache hit)
    q = load_lane_cell(key, workdir=tmp_path)
    assert q["total"] == p["total"] and q["energies"] == p["energies"]
    assert run_lane_cell(key, mol, mf, CAMPAIGN["C"].term,
                         workdir=tmp_path) == p
    # refusals: a free-route key, a sector disagreeing with the fold,
    # a fold shell the basis lacks
    with pytest.raises(ValueError, match="lane cell"):
        run_lane_cell(CellKey("C", "2p", (2, 0)), mol, mf, (1, "g"),
                      workdir=tmp_path)
    with pytest.raises(ValueError, match="active electrons"):
        run_lane_cell(CellKey("C", "2p@1s2s", (3, 1), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)
    with pytest.raises(ValueError, match="no 5s shell"):
        run_lane_cell(CellKey("C", "2p@1s5s", (2, 0), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)


@fast
def test_lane_emitter_deterministic():
    """The ladder emitter's generated block: byte-identical across
    calls from the same payloads (the emitter REPORTS, the ladder pins
    certify), the treatment columns carry the route/method identity,
    ΔE anchors on the primary leg with eV twins, a scrambling leg
    reads the ⟨L²⟩ RANGE (no verdict column), a NaN commutator
    (past the dense-meter cap) reads an em-dash, a leg without a
    payload renders pending, and a missing anchor leg refuses."""
    import importlib.util
    from pathlib import Path as _Path

    loader = importlib.util.spec_from_file_location(
        "build_lane_tables",
        _Path(__file__).resolve().parents[2] / "scripts"
        / "build_lane_tables.py")
    mod = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(mod)

    cs = CAMPAIGN["Sc"]
    keys = lane_keys(cs)
    payloads = {k.address: {"total": -759.7274, "n_det": 90,
                            "l2": [6.0] * 5, "s2": [0.75] * 5,
                            "ptp": 1e-13,
                            "comm_norm_h_lz": 1e-14,
                            "comm_norm_h_l2": 1e-13} for k in keys}
    block = mod.section(cs, payloads)
    assert block == mod.section(cs, payloads)
    lines = [l for l in block.splitlines() if l.startswith("|")]
    assert len(lines) == 2 + len(keys), lines   # head + sep + the legs
    assert lines[0].startswith("| route | method | window | CAS | fold"
                               " | dets | E, Eh | ΔE, mEh | E, eV |"
                               " ΔE, eV |")
    assert lines[2].startswith("| pinned | CI | 3d4s | (3,6) | Ar |")
    assert "| 0 |" in lines[2]                  # the anchor's ΔE
    assert f"| −{759.7274 * 27.2114:.4f} |" in lines[2]   # the eV twin
    assert lines[3].startswith("| pinned | CI | 3d3p | (7,8) | semi |")
    # the extension legs render in declaration order after the pair
    assert [l.split("|")[3].strip() for l in lines[3:]] == [
        "3d3p", "3d3p4d", "3d3p4s", "3s3p3d"]

    no_payload = dict(payloads)
    del no_payload[keys[1].address]
    assert "pending" in mod.section(cs, no_payload)
    no_anchor = dict(payloads)
    del no_anchor[keys[0].address]
    with pytest.raises(SystemExit, match="anchor"):
        mod.section(cs, no_anchor)

    # a scrambling leg (the planned-window class): the meters speak --
    # the ⟨L²⟩ RANGE carries the scramble, no verdict column exists
    scramble = dict(payloads)
    scramble[keys[1].address] = {**payloads[keys[1].address],
                                 "l2": [6.0, 12.0, 6.0, 6.0, 30.0]}
    sec = mod.section(cs, scramble)
    assert "6.0–30.0" in sec and "not a" not in sec

    # past the dense-meter cap: the commutators read em-dash
    cap = dict(payloads)
    cap[keys[1].address] = {**payloads[keys[1].address],
                            "comm_norm_h_lz": float("nan"),
                            "comm_norm_h_l2": float("nan")}
    row = [l for l in mod.section(cs, cap).splitlines()
           if l.startswith("| pinned | CI | 3d3p ")][0]
    assert row.endswith("| — | — |")


@fast
def test_lattice_emitter_deterministic():
    """The lattice emitter's generated section: byte-identical across
    calls from the same payloads, the CAS label carries the DECLARED
    sector (the valence series reads CAS(3,·), not the free-ion
    CAS(2,·)), an unrun cell renders pending in its own column, and
    the row set equals the declared spans in declaration order."""
    import importlib.util
    from pathlib import Path as _Path

    loader = importlib.util.spec_from_file_location(
        "build_lattice_tables",
        _Path(__file__).resolve().parents[2] / "scripts"
        / "build_lattice_tables.py")
    mod = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(mod)

    cs = CAMPAIGN["Sc"]
    keys = lattice_keys(cs)
    payloads = {k.address: {"n_det": 90,
                            "energies_pyscf": [-759.7296],
                            "ground_l2": 6.0} for k in keys}
    sec = mod.section(cs, payloads)
    assert sec == mod.section(cs, payloads)
    lines = [l for l in sec.splitlines() if l.startswith("|")]
    assert lines[0] == mod.HEAD and lines[1] == mod.SEP   # well-formed
    assert lines[1].count("|---") == 7 and "||" not in lines[1]
    assert len(lines) == 2 + len(cs.lattice_windows)
    assert lines[2].startswith("| 3d | CAS(3,5) | 90 |")
    assert lines[3].startswith("| 3d4s | CAS(3,6) | 90 |")
    assert "sector (2, 1)" in sec

    # an unrun cell is pending IN ITS COLUMN (the row survives)
    partial = dict(payloads)
    del partial[keys[1].address]          # 3d's MO cell
    row = [l for l in mod.section(cs, partial).splitlines()
           if l.startswith("| 3d |")][0]
    cells_ = row.split("|")
    assert cells_[6].strip() == "pending" and cells_[7].strip() == ""
    assert cells_[4].strip() == "−759.7296"


@fast
def test_method_ladder_row_policy():
    """The methods renderer's ROW POLICY (the compression): run cells
    only -- pending and void cells leave the table; within a group the
    AO head absorbs MO/NO runs agreeing within EPS (named in its start
    cell), disagreement keeps its own row; a window with no run cells
    renders the POINTER line, not a pending grid; byte-deterministic."""
    import importlib.util
    from pathlib import Path as _Path

    loader = importlib.util.spec_from_file_location(
        "build_methods_tables",
        _Path(__file__).resolve().parents[2] / "scripts"
        / "build_methods_tables.py")
    mod = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(mod)

    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)          # the primary window's factor grid

    def payload(e):
        return {"energies_pyscf": [e], "ground_l2": 12.0,
                "ground_lz2": 4.0, "comm_norm_h_lz": 1e-9,
                "comm_norm_h_l2": 1e-8, "window": {"d": 1.0},
                "electrons": {"d": 2.0}}

    # every cell ran and the optimizer start-triplets agree: full
    # collapse (every group renders one AO row naming +mo,no)
    agree = {k.address: payload(-848.4) for k in keys}
    sec = mod.window_section("H", "3d4s", keys, agree, -848.4065)
    assert sec == mod.window_section("H", "3d4s", keys, agree, -848.4065)
    rows = [r for r in sec.splitlines() if r.startswith("|")]
    ao = [r for r in rows if r.startswith("| CASSCF, free (")][0]
    assert "(+mo,no)" in ao
    assert not any(", MO" in r or ", NO" in r for r in rows)  # collapsed

    # a DISAGREEING mo run keeps its own row (the basin class)
    basin = dict(agree)
    mo = next(k for k in keys if k.method == "casscf" and k.mask is None
              and k.substrate == "pol" and k.start == "mo")
    basin[mo.address] = payload(-848.35)
    rows = [r for r in mod.window_section(
        "H", "3d4s", keys, basin, -848.4065).splitlines()
        if r.startswith("|")]
    assert sum(", MO" in r for r in rows) == 1

    # pending cells leave the table; a void cell is not a run either
    sparse = {k.address: payload(-848.4) for k in keys
              if k.start == "ao"}
    sec = mod.window_section("H", "3d4s", keys, sparse, -848.4065)
    rows = [r for r in sec.splitlines() if r.startswith("|")]
    assert all("(+" not in r for r in rows)      # nothing to collapse
    void = dict(agree)
    void[keys[0].address] = {"void_reason": "scenario"}
    rows = [r for r in mod.window_section(
        "H", "3d4s", keys, void, -848.4065).splitlines()
        if r.startswith("|")]
    assert not any("raw" in r for r in rows)     # the void CASCI row is gone

    # no run cells at all: the POINTER line, no table
    sec = mod.window_section("H", "3d4d", window_grid(cs, "3d4d"),
                             {}, -848.4065)
    assert "regen_cells.py Ti --window 3d4d" in sec
    assert not [r for r in sec.splitlines() if r.startswith("|")]
