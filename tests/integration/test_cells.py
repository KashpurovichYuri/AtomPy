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
    route_keys,
    spine_keys,
    void_reason,
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
        CellKey("U14e", "9f10f11f12d13dno2", (4, 0), method="sa-casscf",
                start="no"),
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
    assert lane_keys(CAMPAIGN["U14e"]) == []


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
def test_u14e_chemistry_declaration():
    """The U14e render chemistry lives in the declaration: the window
    display names (chemical shells -- the uncontracted basis counts
    radial primitives), the fold shell sets (the (5,1) promoted
    window folds 6s6p -- four pairs -- the (4,0) ones 6s6p7s), the
    chemical filling note, and the primary-family suppression."""
    spec = CAMPAIGN["U14e"]
    assert spec.window_names["10f13d"] == "5f6d"
    assert spec.window_names["10f13d7s"] == "5f6d7s"
    assert spec.window_names["10f13d14p"] == "5f6d6p"
    assert spec.fold_chemistry == {5: "6s6p7s", 4: "6s6p"}
    assert spec.render_primary is False
    assert spec.anchor_note.startswith("forced filling: 6s²6p⁶7s²")


@fast
def test_void_reason_voids_the_u14e_optimizer_class():
    """U14e's widened gate: the multi-solution class is void at the key
    on EVERY window (the primary included -- ``run_casscf=False``'s
    spirit made explicit) and the frac folds are void with it (the
    frac reference has no converged basin on the ECP atom); exactly
    the 20 remaining CASCI anchors run (the primary's 5 + the three
    promoted windows' 15 -- the tier families retired from the
    declaration 2026-11-06)."""
    spec = CAMPAIGN["U14e"]
    sc = scenarios_for(spec)
    keys = grid_keys(spec) + planned(spec)
    assert len(keys) == 172
    for k in keys:
        void = void_reason(k, sc) is not None
        expect = k.method != "casci" or k.fold == "frac"
        assert void == expect, k.address
    assert sum(not void_reason(k, sc) for k in keys) == 20


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
    # + the moc block (2 CASCI + 8 optimizer keys)
    assert len(grid) == 43
    assert len(set(addrs)) == 43
    assert addrs[:4] == [           # the CASCI fold block, pur third
        "ti_3d4s_3-1_casci_ao_-_pol_raw",
        "ti_3d4s_3-1_casci_ao_-_pol_avg",
        "ti_3d4s_3-1_casci_ao_-_pur_avg",
        "ti_3d4s_3-1_casci_ao_-_pol_frac",
    ]
    assert addrs[4:6] == [          # the character-selected CASCI pair
        "ti_3d4s_3-1_casci_moc_-_pol_avg",
        "ti_3d4s_3-1_casci_moc_-_pur_avg",
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
def test_planned_sector_override_shrinks_the_fold():
    """A planned window's declared sector overrides the case default in
    its whole grid -- the promoted-pair windows' fold count follows the
    sector (ncore_of(14, 6) = 4: the full core minus the 7s pair);
    windows without the override stay at the default, and the NO-form
    keys inherit their window's sector."""
    keys = planned(CAMPAIGN["U14e"])
    by_window = {}
    for k in keys:
        by_window.setdefault(k.window, []).append(k)
    assert all(k.sector == (4, 0) for k in by_window["10f13d"])
    assert all(k.sector == (4, 0) for k in by_window["10f13d14p"])
    assert all(k.sector == (5, 1) for k in by_window["10f13d7s"])
    # the tier windows (no override) and the NO-forms keep their
    # window's sector
    assert all(k.sector == (4, 0)
               for w, ks in by_window.items() if w.startswith("9f")
               for k in ks)
    no7 = [k for k in by_window if k.startswith("10f13d7sno")]
    assert no7 and all(k.sector == (5, 1)
                       for k in by_window[no7[0]])


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


@fast
def test_sa_budget_longer_than_sector_falls_back(carbon, tmp_path):
    """A declared sa_weights block LONGER than the sector cannot host
    its roots (the occupation-group windows: CAS(3,2) cannot host a
    term's 2L+1) -- the run falls back to the MEASURED ground block
    (the default rule, window-local) instead of crashing the solver."""
    from atomterm.cells import run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="sa-casscf")
    facts = run_cell(key, mf, workdir=tmp_path,
                     budgets={"sa_weights": [1 / 11] * 11})
    n = len(facts.energies_pyscf)
    assert 1 < n < 11        # the measured block anchors, not the declared


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
    # the optimizer's macro-cycle count persists (the callback is the
    # instrument -- PySCF stores no count); CI runs read 0
    assert data.cycles > 0
    ci_key = CellKey("C", "2p", (1, 1), method="casci")
    run_cell(ci_key, mf, workdir=tmp_path)
    assert load_cell(ci_key, tmp_path).cycles == 0


@slow
def test_kernel_cell_payload_meters_match_the_dense_path(carbon, tmp_path):
    """The apply-based meters of a kernel-only payload read the STORED
    window through the run's own fold resolution: the ground ⟨L²⟩
    agrees with the dense analysis payload's, and the state residual
    is present at the certificate floor (the free route's big-cell
    columns; here pinned at a small cell where both paths live)."""
    from atomterm.casscf import run_casci
    from atomterm.cells import (cell_payload, kernel_cell_payload,
                                run_cell)
    from atomterm.windows import select_valence_ao

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="casci",
                  substrate="pur", fold="avg")
    facts = run_cell(key, mf, workdir=tmp_path)
    dense = cell_payload(key, facts)
    kern = kernel_cell_payload(key, mol, mf, workdir=tmp_path)
    assert kern["energies_pyscf"] == dense["energies_pyscf"]
    assert abs(kern["ground_l2"] - dense["ground_l2"]) < 1e-9
    # observed on the frozen carbon: 1.2e-16 / 5.0e-16 (the apply floor)
    assert 0.0 <= kern["res_h_lz"] < 1e-8
    assert 0.0 <= kern["res_h_l2"] < 1e-8
    # no substrate given: the kernel-only shape (NaN meters, no
    # residual) -- the free-standing energies class
    bare = kernel_cell_payload(key, workdir=tmp_path)
    assert "res_h_lz" not in bare and "res_h_l2" not in bare


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
        assert abs(v - 2.0) < 1e-8
    # the state residual: the closed window's term components are
    # symmetry eigenstates -- the residual sits at the apply floor
    # (observed on the frozen D2h twin: 1.1e-16 / 2.2e-16)
    assert 0.0 <= p["res_h_lz"] < 1e-8
    assert 0.0 <= p["res_h_l2"] < 1e-8
    # the store round-trips (and a fresh run is a cache hit)
    q = load_lane_cell(key, workdir=tmp_path)
    assert q["total"] == p["total"] and q["energies"] == p["energies"]
    assert "mo_coeff" not in q and "cycles" not in q  # CI: no optimizer
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
    with pytest.raises(ValueError, match="no 5s band"):
        run_lane_cell(CellKey("C", "2p@1s5s", (2, 0), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)


@slow
def test_run_lane_cell_multicolumn_fold(carbon_d2h, carbon, tmp_path):
    """The multi-column fold end to end: the SAME lane leg with the
    shell layer declaring the 1s2s chemistry over three s bands (two
    pairs over a three-column span) still hosts its term EXACTLY (the
    span's own closed-core Roothaan is m-scalar like the 1:1 raw
    span), while the core axes read the CHEMISTRY: the payload's pair
    count is 2 (never the span width 3) with the Roothaan report.  A
    declared band shared by window and fold refuses at the ROW level
    (the layer's own ambiguity, not the labels')."""
    cfg, mol, mf, _ = carbon_d2h
    key = lane_keys(CAMPAIGN["C"])[0]
    layer = {"2s": ("2s", "3s")}
    p = run_lane_cell(key, mol, mf, CAMPAIGN["C"].term, shell_bands=layer,
                      workdir=tmp_path)
    assert p["core_pairs"] == 2
    assert p["core_cycles"] >= 1 and p["core_comm"] < 1e-10
    assert p["ptp"] < 1e-9
    for v in p["l2"]:
        assert abs(v - 2.0) < 1e-8              # 3P: L(L+1) = 2 exact
    for v in p["s2"]:
        assert abs(v - 2.0) < 1e-8
    assert 0.0 <= p["res_h_lz"] < 1e-8
    # the 1:1 leg (no layer) stays the certified record: same term,
    # its own model core (the energies differ -- DIFFERENT folds)
    q = run_lane_cell(key, mol, mf, CAMPAIGN["C"].term, workdir=tmp_path,
                      fresh=True)
    assert q["core_pairs"] == 2 and q["core_cycles"] == 0
    assert q["total"] != p["total"]
    # the partition is CHEMICAL: a shell in both sides refuses...
    with pytest.raises(ValueError, match="sit in both"):
        run_lane_cell(CellKey("C", "2p2s@1s2s", (2, 0), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)
    # ...while shared BANDS are physical (the fold's 2s realized over
    # the window's 3s band): the projection owns the shared rows
    # -- the leg RUNS, its own model core
    p_shared = run_lane_cell(
        CellKey("C", "2p3s@1s2s", (2, 0), substrate="pur"), mol, mf,
        (1, "g"), shell_bands={"2s": ("2s", "3s")}, workdir=tmp_path)
    assert p_shared["core_pairs"] == 2
    assert p_shared["ptp"] < 1e-9
    for v in p_shared["l2"]:
        assert abs(v - 2.0) < 1e-8
    # a band shared by TWO FOLD SHELLS duplicates its rows -- the span
    # is the UNION (a duplicated row would make the Gram exactly
    # singular; pinned after the Scunc lane hit it live)
    p_union = run_lane_cell(
        CellKey("C", "2p@1s2s", (2, 0), substrate="pur"), mol, mf,
        (1, "g"), shell_bands={"1s": ("1s", "2s"), "2s": ("2s",)},
        workdir=tmp_path)
    assert p_union["core_pairs"] == 2
    assert p_union["ptp"] < 1e-9


@slow
def test_unc_lanes_load_exact_from_the_store():
    """The validation lanes re-verified from the COMMITTED stores (the
    heavy-pin pattern: fresh legs run in the campaign; the verdicts
    load and re-check their meters): Scunc hosts 2D and Tiunc 3F
    EXACTLY through the multi-column folds -- every component at its
    L(L+1) and S(S+1), degenerate to the ptp floor, the core pair
    count chemical (9), the Roothaan converged, and the dense
    adjudication's corrections recorded (the pinned solver's
    near-degenerate skips -- the branch the live runs exercised)."""
    from atomterm.cells import load_lane_cell

    cases = (
        ("Scunc", "scunc_3d4s@ar_2-1_casci_ao_-_pur_avg", 6.0, 0.75),
        ("Tiunc", "tiunc_3d4s@ar_3-1_casci_ao_-_pur_avg", 12.0, 2.0),
    )
    for key_token, addr, l2, s2 in cases:
        p = load_lane_cell(CellKey.from_address(addr))
        assert p["core_pairs"] == 9                    # the chemistry
        assert p["core_comm"] < 1e-9                   # the Roothaan
        assert p["ptp"] < 1e-9                         # degenerate
        for v in p["l2"]:
            assert abs(v - l2) < 1e-8
        for v in p["s2"]:
            assert abs(v - s2) < 1e-8
        assert 0.0 <= p["res_h_lz"] < 1e-8
        assert p["anchor_corrected"]                    # the live skip
        assert p["anchor_dev"] > 1e-3                   # honestly recorded


@slow
def test_confined_core_reproduces_the_ion_rhf(carbon):
    """The EXTERNAL anchor of the confined-core construction (the one
    certificate the internal identities cannot give): on the WHOLE AO
    space the mini-Roothaan over n_pairs pairs is the closed-shell
    ion of the same charge in the same basis -- its RHF energy, to
    the Roothaan convergence floor.  The frozen reference supplies
    the field pieces; the ion SCF is the independent construction."""
    from pyscf import gto, scf

    from atomterm.core import closed_shell_core_energy, confined_core_orbitals
    from atomterm.windows import lowdin

    _, mol, mf, _ = carbon
    ion_mol = gto.M(atom="C 0 0 0", basis="ccpvdz", charge=2, verbose=0)
    ion = scf.RHF(ion_mol)
    ion.conv_tol = 1e-12
    ion.kernel()
    assert ion.converged
    c_span = lowdin(np.eye(mol.nao_nr()), mf.get_ovlp())
    c_core, report = confined_core_orbitals(mf, c_span, 2)
    assert report["comm"] < 1e-10
    assert abs(closed_shell_core_energy(mf, c_core)
               - ion.e_tot) < 1e-8          # C2+ closed shell, 2 pairs


@slow
def test_route_axis_twins(carbon, carbon_d2h, tmp_path):
    """The route axis: route_keys spawns every leg's four twins (the
    substrate axis naming the route); the SPAN-IDENTITY certificate --
    the free CI twin (run_cell, the C1 driver) reproduces the pinned
    leg's ground component to solver tolerance (same window, same
    fold, different construction); run_cell's fold branch refuses the
    pur class (run_lane_cell's), raw/frac, non-ao starts, and wrong
    sectors; the pinned SA runner delivers the two-pass payload with
    its term exact."""
    from dataclasses import replace

    from atomterm.cells import cell_payload, run_cell
    _, mol_c1, mf_c1, _ = carbon
    _, mol_d2h, mf_d2h, _ = carbon_d2h
    cs = CAMPAIGN["C"]
    leg = lane_keys(cs)[0]

    # the derivation: 1 pinned CI + 4 twins per leg, routes by substrate
    rk = route_keys(cs)
    assert len(rk) == 5 and rk[0].address == leg.address
    tails = sorted(k.address[len(leg.address) - len("casci_ao_-_pur_avg"):]
                   for k in rk[1:])
    assert tails == ["casci_ao_-_pol_avg", "casscf_ao_-_pol_avg",
                     "sa-casscf_ao_-_pol_avg", "sa-casscf_ao_-_pur_avg"]

    # the span-identity certificate at CI level: the free driver's
    # ground TOTAL equals the pinned leg's TERM total (the components
    # are degenerate -- same span, same fold, different construction;
    # observed on the frozen twins: 7.1e-15)
    pinned = run_lane_cell(leg, mol_d2h, mf_d2h, cs.term, workdir=tmp_path)
    free = cell_payload(
        replace(leg, substrate="pol"),
        run_cell(replace(leg, substrate="pol"), mf_c1, mol_c1,
                 workdir=str(tmp_path)))
    assert abs(free["energies_pyscf"][0] - pinned["total"]) < 1e-8

    # the fold branch's refusals
    for kw, msg in (
            (dict(substrate="pur"), "run_lane_cell"),
            (dict(substrate="pol", fold="raw"), "raw/frac"),
            (dict(substrate="pol", start="mo"), "ao/pol/avg"),
            (dict(substrate="pol"), "active electrons")):
        k = (replace(leg, window="2p@1s2s", sector=(3, 1), **kw)
             if msg == "active electrons" else replace(leg, **kw))
        with pytest.raises(ValueError, match=msg):
            run_cell(k, mf_c1, mol_c1, workdir=str(tmp_path))

    # the pinned SA runner: the two-pass payload, term exact
    sa = run_lane_cell(replace(leg, method="sa-casscf"), mol_d2h, mf_d2h,
                       cs.term, workdir=tmp_path)
    assert "ensemble" in sa and sa["ptp"] < 1e-8  # observed: 1.3e-14
    # the SA TOTAL is the kernel ensemble (the optimizer rotates
    # same-irrep core/window columns -- the readout's absolutes are of
    # the run's own partition; only the kernel bookkeeping compares)
    assert sa["total"] == sa["ensemble"]
    assert all(abs(v - 2.0) < 1e-7 for v in sa["l2"])   # 3P exact
    assert sa["total"] < min(pinned["energies"]) + 1.0  # relaxation helps
    # the converged window PERSISTS: the stored payload round-trips
    # the full MO matrix (restart + orbital-level route comparisons)
    # and the optimizer's cycle count
    sa_q = load_lane_cell(replace(leg, method="sa-casscf"),
                          workdir=tmp_path)
    assert sa_q["mo_coeff"].shape == (mol_d2h.nao_nr(),) * 2
    assert np.array_equal(sa_q["mo_coeff"], sa["mo_coeff"])
    assert sa_q["cycles"] == sa["cycles"] > 0
    # the convergence precision persists beside the count (the conv
    # column's data; 0.0 is the legitimately dead-still final cycle)
    assert sa_q["conv_dE"] == sa["conv_dE"] >= 0.0


def ladder_module():
    """The methods-table builder, loaded once per call site (a script,
    not a package import)."""
    import importlib.util
    from pathlib import Path as _Path

    loader = importlib.util.spec_from_file_location(
        "build_methods_tables",
        _Path(__file__).resolve().parents[2] / "scripts"
        / "build_methods_tables.py")
    mod = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(mod)
    return mod


def _ladder_payload(e, roots=1):
    """A minimal run payload for ladder rendering."""
    return {"energies_pyscf": [e] * roots, "n_det": 90,
            "ground_l2": 12.0, "comm_norm_h_lz": 1e-9,
            "comm_norm_h_l2": 1e-8, "win_lz": 1e-15,
            "win_lperp": 2e-15}


@fast
def test_ladder_table_is_byte_deterministic():
    """THE LADDER (one table per atom, between the markers): the
    section builds byte-identically from the same store, once."""
    mod = ladder_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _ladder_payload(-848.4) for k in keys}
    text = "\n".join(mod.ladder_section(cs, agree, -848.4065))
    assert text == "\n".join(mod.ladder_section(cs, agree, -848.4065))
    assert text.count(mod.BEGIN) == 1 and text.count(mod.END) == 1


@fast
def test_ladder_anchors_lead_and_families_follow():
    """The ROHF anchor leads; the UHF anchor row follows when the
    store carries the record (DeltaE against the same zero, the meter
    in the note); the first family follows directly."""
    mod = ladder_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _ladder_payload(-848.4) for k in keys}

    def body(cells, uhf=None):
        return "\n".join(mod.ladder_section(cs, cells, -848.4065, uhf=uhf))

    rows = [r for r in body(agree).splitlines() if r.startswith("|")]
    assert rows[2].startswith("| ROHF (anchor) |")
    with_uhf = body(agree, uhf={"e_tot": -848.41, "s2": 4.02,
                                "sigma_min": 0.98})
    rows_uhf = [r for r in with_uhf.splitlines() if r.startswith("|")]
    assert rows_uhf[2].startswith("| ROHF (anchor) |")
    cells_uhf = [c.strip() for c in rows_uhf[3].split("|")[1:-1]]
    assert cells_uhf[0] == "UHF (anchor)"
    assert cells_uhf[7] == "−3.5"
    assert "⟨S²⟩=4.0200" in cells_uhf[17]
    assert rows_uhf[4].startswith("| CASCI |")  # the first family follows


@fast
def test_ladder_every_run_cell_a_row_no_starts_excluded():
    """EVERY run cell a row -- the four fold variants' CI rows beside
    their optimizer variants, the method bare + restriction only, the
    window cell in slash form; NO starts are not ladder rows; free
    optimizer rows carry no lane marker."""
    mod = ladder_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _ladder_payload(-848.4) for k in keys}
    text = "\n".join(mod.ladder_section(cs, agree, -848.4065))
    rows = [r for r in text.splitlines() if r.startswith("|")]
    ci = [r for r in rows if r.startswith("| CASCI | free | 3d4s / AO |")]
    assert len(ci) == 4            # the four fold variants' CI rows
    assert not any("/ NO |" in r for r in rows)   # NO starts excluded
    opt = [r for r in rows
           if r.startswith("| CASSCF | free | 3d4s / AO |")]
    assert len(opt) == 2 and not any("·lane" in r for r in opt)  # Ar + Ar·pur substrates


@fast
def test_ladder_meter_and_ptp_columns_render():
    """The closure columns carry the span certificate; ptp computes
    where a root set exists and stays dash on single roots; the conv
    column carries precision/count."""
    mod = ladder_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _ladder_payload(-848.4) for k in keys}

    def body(cells):
        return "\n".join(mod.ladder_section(cs, cells, -848.4065))

    text = body(agree)
    rows = [r for r in text.splitlines() if r.startswith("|")]
    opt = [r for r in rows
           if r.startswith("| CASSCF | free | 3d4s / AO |")]
    assert opt[0].split("|")[14].strip() == "1e-15"
    assert opt[0].split("|")[15].strip() == "2e-15"
    assert opt[0].split("|")[13].strip() == "—"  # the single root: no set
    sa = next(k for k in keys if k.method == "sa-casscf"
              and k.mask is None and k.start == "ao"
              and k.substrate == "pol")
    traced = dict(agree)
    traced[sa.address] = {**_ladder_payload(-848.4, roots=5),
                          "conv_dE": 2.3e-9, "cycles": 5}
    row = [r for r in body(traced).splitlines()
           if r.startswith("| SA-CASSCF | free | 3d4s / AO |")][0]
    assert "| 1e-9/5 |" in row and "| 0e+00 |" in row


@fast
def test_ladder_cycle_cap_note():
    """The cap note fires by IDENTITY on payloads carrying the cap
    (the store's ``max_cycle`` field -- the run's actual mc cap) and
    keeps the legacy ``== 200`` key for payloads predating the field;
    a mid-cap count never flags."""
    mod = ladder_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _ladder_payload(-848.4) for k in keys}
    sa = next(k for k in keys if k.method == "sa-casscf"
              and k.mask is None and k.start == "ao"
              and k.substrate == "pol")

    def row_of(extra):
        cells = dict(agree)
        cells[sa.address] = {**_ladder_payload(-848.4, roots=5),
                             "conv_dE": 2.3e-9, **extra}
        row = [r for r in "\n".join(
                   mod.ladder_section(cs, cells, -848.4065)).splitlines()
               if r.startswith("| SA-CASSCF | free | 3d4s / AO |")][0]
        return "cycle cap" in row

    assert row_of({"cycles": 5, "max_cycle": 5})       # identity fires
    assert not row_of({"cycles": 5, "max_cycle": 200})  # below the cap
    assert row_of({"cycles": 200})                     # legacy key fires
    assert not row_of({"cycles": 5})                   # legacy, mid-cap


@fast
def test_ladder_ci_only_store_renders_planned_pointers():
    """No optimizer run cells anywhere: the CI rows still render, the
    gated planned windows follow as pointers."""
    mod = ladder_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    text = "\n".join(mod.ladder_section(
        cs, {k.address: _ladder_payload(-848.4) for k in keys
             if k.method == "casci"}, -848.4065))
    assert "| CASCI | free | 3d4s / AO |" in text
    assert "regen_cells.py Ti --window 3d4d" in text
    assert "| CASSCF" not in text and "| SA-CASSCF" not in text


@fast
def test_method_ladder_lane_families():
    """The one table carries the lane families: pinned CI, free CI,
    the free CASSCF/SA, and the pinned SA in one family block; the
    tagged folds carry the lane marker; un-run lane families are
    absent; families are separated by the rule row."""
    mod = ladder_module()
    cs = CAMPAIGN["Sc"]

    def p(e, pinned=False, ci=False):
        d = {"n_det": 204490, "comm_norm_h_lz": float("nan"),
             "comm_norm_h_l2": float("nan"), "res_h_lz": 1e-14,
             "res_h_l2": 2e-14, "win_lz": 1e-15, "win_lperp": 1e-15}
        if pinned or ci:
            d.update({"total": e, "l2": [6.0] * 5, "ptp": 5e-8,
                      "energies": [e] * 5})
        else:
            d["energies_pyscf"] = [e]
            d["ground_l2"] = 6.0
        return d

    routes = route_keys(cs)
    fam = [k for k in routes if k.window == "3d3p4d@semi"]
    cells = {k.address: p(-759.85,
                          pinned=(k.substrate == "pur"
                                  and k.method == "sa-casscf"),
                          ci=(k.method == "casci"))
             for k in fam}
    text = "\n".join(mod.ladder_section(cs, cells, -759.7363))
    lines = text.splitlines()
    assert sum(l.startswith("| method |") for l in lines) == 1
    assert lines[3].startswith("| ROHF (anchor) |")
    rules = [l for l in lines
             if set(l.replace("|", "").replace(" ", "")) == {"-"}]
    assert len(rules) == 1    # the SEP only: the lone family leads
    assert sum(l.startswith("| CASCI | pinned | 3d3p4d / AO |")
               for l in lines) == 1
    assert sum(l.startswith("| CASCI | free | 3d3p4d / AO |")
               for l in lines) == 1
    assert sum(l.startswith("| CASSCF | free | 3d3p4d / AO |")
               for l in lines) == 1
    assert sum(l.startswith("| SA-CASSCF | pinned | 3d3p4d")
               for l in lines) == 1
    assert "·lane" in [l for l in lines
                       if l.startswith("| SA-CASSCF | pinned | 3d3p4d")
                       ][0]
    assert "3d3p / AO" not in text       # un-run families absent

