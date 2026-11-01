"""The symmetry lane on the open d shells -- the stress case, and
where the two routes genuinely diverge.

The multiplicity-2 branching goes live through the runner (Sc's 2D:
the Ag block hosts TWO components; Ti's 3F: each B_g hosts two), and
the anchor holds at machine precision as on the p shells.

The head-to-head, though, splits by atom.  On Sc (d1) the routes
agree as on C/Si.  On Ti (d2) the free SA escapes the declared
{3d,4s} window by a CORE<->ACTIVE SWAP, measured on the inactive
block (Ar s:3 p:6 -> s:4 p:5): the 4s pair rotates INTO the fold,
one semi-core 3p orbital -- with its electron pair, 0.998 pure by
AO content -- rotates INTO the active, buying the 3p-3d semi-core
correlation at the price of a torn m-set (certificates O(1),
multiplet split at mEh) while the roots still READ as the term --
membership alone is not the gate, the multiplet degeneracy is.
The pinned route's D2h-blocked rotations forbid the swap (the swap
is cross-irrep: 4s is Ag, 3p is u): it stays in the declared
identity, components degenerate to 4e-7, certificates at 1e-5 --
the guaranteed component set at the price of ~5 mEh of ensemble
energy.  The symmetry-vs-variation trade, measured.

Windows, sectors, and weights derive from the CAMPAIGN declaration;
the ground-term table (L, 2S, parity per atom) is test data, the
same kind of row test_atoms carries."""

import numpy as np
import pytest
from pyscf import ao2mo, fci, gto, scf

from atomterm.casscf import run_irrep_casci, run_irrep_sa, run_sa_casscf
from atomterm.cases import CAMPAIGN
from atomterm.core import closed_shell_core_energy
from atomterm.detbasis import slater_condon_matrix
from atomterm.irreps import IRREP_IDS, branch_d2h, det_irreps, orb_irreps
from atomterm.operators import (
    active_space_integrals,
    build_operator_parts,
    build_operators,
    ci_in_detbasis_order,
    window_closure,
)
from atomterm.system import sector_na_nb
from atomterm.windows import (
    build_full_mo,
    lowdin,
    occupied_core_idx,
    select_shell_ao,
    select_valence_ao,
    shell_weights,
)

pytestmark = pytest.mark.slow

# the ground term per atom: L, S(S+1), and the configuration parity
# (d^k -> g); 2S arrives from the CAMPAIGN row's spin
_TERMS = {
    "Sc": dict(L=2, s2=0.75, parity="g"),
    "Ti": dict(L=3, s2=2.0, parity="g"),
}

# the free route root-set membership threshold: within-term mixes read
# <L^2> = L(L+1) exactly, neighbouring terms differ by >= 2L+2, and
# the free route's optimizer-tolerance window defect shifts <L^2> at ~1e-3
_MEMBER_TOL = 0.5
# the multiplet-splitting gate between the two branches: a formed
# term splits at <= 1e-6 (observed Sc 3.1e-7, the pinned route 3.8e-7); the
# escaped/tilted draw splits at mEh (observed Ti 3.1e-3)
_SPLIT_TOL = 1e-5

# observed bands (fresh draws, 2026-10-31): the free/pinned
# ensembles agree to 1.1e-8 on the symmetric draw; the pinned route
# closes at 6.9e-6..1.4e-5
# (both atoms); the symmetric free-route window carries an f-direction
# drift the energy cannot feel (Sc: certificates 5.2e-3 with
# ensembles agreeing to 1e-8) -- hence the generous certificate band
_ENS_BAND = 1e-6
_CERT_FREE = 2e-2
_CERT_PINNED = 1e-4


def _root_vectors(ci):
    """A solver's return as a list of root vectors (nroots = 1 returns
    a bare array, nroots > 1 a stacked one)."""
    ci = np.asarray(ci)
    return [ci] if ci.ndim == 2 else list(ci)


@pytest.fixture(scope="module", params=sorted(_TERMS))
def tm_lane(request):
    """Both substrates of one reference (C1 = the free route, D2h
    twin = the pinned route) with the window, orbsyms, and our
    matrices -- everything but
    the term table derived from the CAMPAIGN declaration.  Any
    converged reference of the right charge works on either side:
    the whole active set sits in the window, so the fold sees only
    the closed noble-gas core."""
    name = request.param
    spec = CAMPAIGN[name]
    term = _TERMS[name]
    na, nb = sector_na_nb(spec.nelec_act, spec.spin)
    mols, mfs = [], []
    for symmetry in (None, "D2h"):
        mol = gto.M(atom=f"{spec.symbol} 0 0 0", basis=spec.basis,
                    spin=spec.spin, symmetry=symmetry, verbose=0)
        mf = scf.ROHF(mol)
        mf.verbose = 0
        mf.conv_tol = 1e-12
        mf.kernel()
        assert mf.converged, name
        mols.append(mol)
        mfs.append(mf)
    mol_c1, mol_d2h = mols
    mf_c1, mf_d2h = mfs
    c_act = np.eye(mol_d2h.nao_nr())[
        :, select_valence_ao(mol_d2h, list(spec.shell_patterns))]
    orbsyms = orb_irreps(mol_d2h, c_act)
    ops = build_operators(mol_d2h, mf_d2h, c_act, spec.ncore, na, nb)
    return name, spec, term, mol_c1, mf_c1, mol_d2h, mf_d2h, c_act, orbsyms, ops, na, nb


def test_tm_runner_anchor(tm_lane):
    # the pinned component roots == the LOWEST eigenvalues of our H on
    # the same blocks (the multiplicity-2 blocks live: the Ag block
    # holds two 2D components on Sc, each B_g two 3F components on
    # Ti); every root reads the term exactly by OUR meters -- the
    # runner's documented "lowest per block = the term" assumption,
    # verified, not assumed; the components are degenerate on the
    # m-scalar fold
    name, spec, term, _, _, mol_d2h, mf_d2h, c_act, orbsyms, ops, na, nb = tm_lane
    l2 = term["L"] * (term["L"] + 1)
    roots = run_irrep_casci(mol_d2h, mf_d2h, c_act, spec.ncore, na, nb,
                            term["L"], term["parity"])
    assert set(roots) == set(branch_d2h(term["L"], term["parity"])), name
    ids = det_irreps(orbsyms, na, nb)
    e_all = []
    for irrep, (e, ci) in roots.items():
        sel = ids == IRREP_IDS[irrep]
        ours = np.linalg.eigvalsh(ops["H"][np.ix_(sel, sel)])[:len(e)]
        assert np.allclose(ours, np.sort(e), atol=1e-8), (name, irrep)
        for k, vec in enumerate(_root_vectors(ci)):
            w = ci_in_detbasis_order(vec, c_act.shape[1], na, nb)
            # meters exact on their vectors (the C-anchor class: 2e-15
            # purity at the raw window)
            assert abs(float(np.real(w @ ops["L2"] @ w)) - l2) < 1e-9, \
                (name, irrep, k)
            assert abs(float(np.real(w @ ops["S2"] @ w)) - term["s2"]) < 1e-9, \
                (name, irrep, k)
        e_all.append(e)
    e_all = np.concatenate(e_all)
    # the components degenerate on the m-scalar fold (observed Sc
    # 1.2e-14, Ti 9.8e-15)
    assert np.ptp(e_all) < 1e-9, name


def test_tm_sa_head_to_head(tm_lane, tmp_path):
    # THE STRESS.  The free route (C1 SA) HOPES its lowest 2L+1 roots
    # are the term's M_L set; the pinned route GUARANTEES them.  The
    # honest gate on the open d shells is the multiplet DEGENERACY,
    # not membership alone: on the recorded Ti draw all 7 free-route
    # roots still read <L^2> = 12 while the free optimizer had
    # performed a CORE<->ACTIVE SWAP (measured below: the inactive
    # block's Ar s:3 p:6 -> s:4 p:5) -- the 4s pair into the fold,
    # one semi-core 3p orbital with its pair into the active,
    # certificates O(1), multiplet split 3.1e-3; the production
    # tables record the same basin.  The pinned route's
    # D2h-blocked rotations cannot follow: the swap is cross-irrep
    # (4s is Ag, 3p is u) -- it stays in the declared s+d identity,
    # degenerate and closed, ~5 mEh above.
    name, spec, term, mol_c1, mf_c1, mol_d2h, mf_d2h, c_act, _, _, na, nb = tm_lane
    ncomp = 2 * term["L"] + 1
    l2 = term["L"] * (term["L"] + 1)
    ncas = c_act.shape[1]
    s_c1 = mf_c1.get_ovlp()

    pinned_sa = run_irrep_sa(mol_d2h, mf_d2h, c_act, spec.ncore, na, nb,
                             term["L"], term["parity"])
    assert pinned_sa["converged"], name
    mo_start_c1 = build_full_mo(
        mf_c1, occupied_core_idx(mf_c1, spec.ncore),
        np.eye(mol_c1.nao_nr())[
            :, select_valence_ao(mol_c1, list(spec.shell_patterns))])
    free_sa = run_sa_casscf(mf_c1, ncas, (na, nb), mo_start_c1,
                            list(spec.sa_weights),
                            str(tmp_path / f"{name.lower()}_sa.chk"))
    ens_free = float(np.mean(free_sa.energies))
    active = slice(spec.ncore, spec.ncore + ncas)

    # the pinned route's component set and l-identity are by
    # construction (one pinned solver per branch irrep; g-irreps only
    # -- a p admixture is cross-irrep, unreachable); its degeneracy,
    # certificates, and term purity are observed bands
    e_pinned = np.concatenate([e for e, _ in pinned_sa["roots"].values()])
    assert np.ptp(e_pinned) < 1e-5, name  # observed 3.8e-7
    labs, w = shell_weights(mol_c1, pinned_sa["mo"][:, active], s_c1)
    assert dict(zip(labs, w.sum(axis=0))).get("p", 0.0) < 1e-6, name
    cert_pinned = window_closure(mol_c1, pinned_sa["mo"][:, active])
    assert cert_pinned["lz"] < _CERT_PINNED, name
    assert cert_pinned["lperp"] < _CERT_PINNED, name
    ops_pinned = build_operators(mol_d2h, mf_d2h, pinned_sa["mo"][:, active],
                                 spec.ncore, na, nb)
    for irrep, (_, ci) in pinned_sa["roots"].items():
        for k, vec in enumerate(_root_vectors(ci)):
            wl = ci_in_detbasis_order(vec, ncas, na, nb)
            assert abs(float(np.real(wl @ ops_pinned["L2"] @ wl)) - l2) < 1e-6, \
                (name, irrep, k)  # observed ~1.1e-10

    # the free route's draw, measured: membership + a formed
    # (degenerate) multiplet + the core<->active swap detector (the
    # inactive block's l-content against the Ar start s:3 p:6)
    ops_free = build_operators(mol_c1, mf_c1, free_sa.mo_coeff[:, active],
                               spec.ncore, na, nb)
    l2_free = []
    for ci in free_sa.ci:
        wl = ci_in_detbasis_order(ci, ncas, na, nb)
        l2_free.append(float(np.real(wl @ ops_free["L2"] @ wl)))
    member = sum(abs(v - l2) < _MEMBER_TOL for v in l2_free) == ncomp
    formed = np.ptp(free_sa.energies) < _SPLIT_TOL
    labs_i, w_i = shell_weights(mol_c1, free_sa.mo_coeff[:, :spec.ncore],
                                s_c1)
    inact = dict(zip(labs_i, w_i.sum(axis=0)))
    labs_a, w_a = shell_weights(mol_c1, free_sa.mo_coeff[:, active], s_c1)
    swapped = float(w_a[:, labs_a.index("p")].sum()) > 0.5
    cert_free = window_closure(mol_c1, free_sa.mo_coeff[:, active])

    def _assert_inact(s, p, tag):
        # the swap detector: the Ar start reads s:3 p:6; the swap
        # draws s:4 p:5 (4s in, one 3p out with its pair)
        assert abs(inact["s"] - s) < 0.05 and abs(inact["p"] - p) < 0.05, \
            (name, tag)

    if member and formed:
        # the symmetric draw (recorded on Sc): the core untouched
        assert not swapped, name
        _assert_inact(3, 6, "symmetric")
        assert abs(ens_free - pinned_sa["e_tot"]) < _ENS_BAND, name
        assert cert_free["lz"] < _CERT_FREE, name
        assert cert_free["lperp"] < _CERT_FREE, name
    else:
        # whichever divergence shape the draw took, the free route's
        # ensemble sits BELOW the pinned route's: by free variation in
        # the swap shape, by displacement in the interleave (the
        # intruder enters below the term root it pushes out).  The
        # free route's certificate is measured, not pinned: the
        # divergent window class is basin chaos, and optimizer output
        # is never pinned where chaos is the object
        assert ens_free < pinned_sa["e_tot"], name
        if swapped:
            # THE MECHANISM, pinned: the core<->active swap -- the
            # 4s pair into the fold, one semi-core 3p orbital (with
            # its pair) into the active (observed 0.998 pure 3p)
            _assert_inact(4, 5, "swap")
            rows_3p = select_shell_ao(mol_c1, "3p")
            col = free_sa.mo_coeff[
                :, spec.ncore + int(np.argmax(w_a[:, labs_a.index("p")]))]
            scol = s_c1 @ col
            assert float(np.dot(col[rows_3p], scol[rows_3p])) > 0.9, name
        else:
            # the interleave shape: the core intact
            _assert_inact(3, 6, "interleave")


# ---------------------------------------------------------------------------
# The CAS ladder: the pinned route's TERM-energy survey across the
# window families (declared / valence / semi-core / the M shell), CI
# level at frozen shells, ONE uniform construction: fold = Loewdin of
# the fold's raw AO shell columns, window placed against it through
# build_full_mo(c_core=) (raw columns carry the documented ~85 mEh
# window-definition error on core-overlapping windows), pinned solver
# per branch component, totals = component e_cas mean +
# closed_shell_core_energy(fold) -- the fold-independent bookkeeping
# ruler (fact 9's fixed-core-block identity), so CROSS-FOLD legs
# compare.  The energies are TERM energies, not sector grounds: each
# window's 2L+1 branch components are solved and averaged.  Everything
# is draw-free (raw AO spans on both sides of the fold), so the legs
# are deterministic.

# the fold spans: the noble-gas prefix, the semi-core fold (3p ACTIVE,
# 4s folded), and the two M-shell-family folds
_AR = ("1s", "2s", "2p", "3s", "3p")
_SEMI = ("1s", "2s", "2p", "3s", "4s")
_M3 = ("1s", "2s", "2p", "4s")
_S3 = ("1s", "2s", "2p", "3s")

# (window shells, fold shells, active electrons, tag) per atom; the
# Ti {3d,4d} leg is read AS DECLARED by the campaign (fold Ar, the 4s
# pair in the dead zone) -- it is not expected to host the term
_LADDER = {
    "Ti": [
        (("3d", "4s"), _AR, 4, "declared"),
        (("3d", "4d"), _AR, 4, "planned"),
        (("3d", "4s", "4p", "4d"), _AR, 4, "valence"),
        (("3d", "3p"), _SEMI, 8, "semicore"),
        (("3d", "3p", "4d"), _SEMI, 8, "semicore_d"),
        (("3s", "3p", "3d"), _M3, 10, "mshell"),
        (("3d", "3p", "4s"), _S3, 10, "semicore_4s"),
    ],
    "Sc": [
        (("3d", "4s"), _AR, 3, "declared"),
        (("3d", "3p"), _SEMI, 7, "semicore"),
    ],
}


def _ladder_leg(mol, mf, spin, win, fold_shells, nelec, L, parity,
                anchor=False):
    """One ladder leg under the uniform construction (see the section
    comment); meters are built per IREP BLOCK -- L2/S2/H are totally
    symmetric, so they never connect determinant irreps and the block
    restriction is exact.  With ``anchor``, the block H is
    eigen-solved and compared against the pinned solver (the direct
    block pattern's convention certificate)."""
    s = mf.get_ovlp()
    n = mol.nao_nr()
    fold_rows = np.concatenate(
        [select_shell_ao(mol, sh) for sh in fold_shells])
    c_fold = lowdin(np.eye(n)[:, fold_rows], s)
    nf = c_fold.shape[1]
    rows = select_valence_ao(mol, list(win))
    mo = build_full_mo(mf, list(range(nf)), rows, c_core=c_fold)
    c = mo[:, nf:nf + len(rows)]
    ncas = c.shape[1]
    na, nb = sector_na_nb(nelec, spin)
    h1, _ = active_space_integrals(mol, mf, c, nf, c_core=c_fold)
    eri = ao2mo.kernel(mol, c)
    orbsyms = orb_irreps(mol, c)
    parts = build_operator_parts(mol, mf, c, nf, na, nb, c_core=c_fold)
    dets = parts["dets"]
    ids = det_irreps(orbsyms, na, nb)
    es, l2s, s2s = [], [], []
    anchor_dev = 0.0
    for irrep, mult in branch_d2h(L, parity).items():
        sol = fci.direct_spin1_symm.FCI(mol)
        sol.orbsym, sol.wfnsym = orbsyms, IRREP_IDS[irrep]
        sol.nroots, sol.conv_tol = mult, 1e-12
        e, ci = sol.kernel(h1, eri, ncas, (na, nb))
        e = np.atleast_1d(e)[:mult]
        es.extend(float(x) for x in e)
        sel = ids == IRREP_IDS[irrep]
        block = [d for d, m in zip(dets, sel) if m]
        w_blk = []
        for vec in _root_vectors(ci):
            w = ci_in_detbasis_order(vec, ncas, na, nb)[sel]
            # the pinned vector lives entirely in its block
            assert abs(np.linalg.norm(w) - 1.0) < 1e-8, irrep
            w_blk.append(w)
        if anchor:
            hb = slater_condon_matrix(block, one_body=parts["H"][0],
                                      two_body=parts["H"][1])
            ours = np.linalg.eigvalsh(hb)[:mult]
            anchor_dev = max(anchor_dev,
                             float(np.max(np.abs(np.sort(e) - ours))))
            del hb
        for meters, op in ((l2s, "L2"), (s2s, "S2")):
            mb = slater_condon_matrix(block, one_body=parts[op][0],
                                      two_body=parts[op][1])
            for w in w_blk:
                meters.append(float(np.real(w @ mb @ w)))
            del mb
    es = np.array(es)
    return dict(total=float(es.mean() + closed_shell_core_energy(mf, c_fold)),
                ptp=float(np.ptp(es)), l2=l2s, s2=s2s, anchor=anchor_dev)


@pytest.fixture(scope="module")
def cas_ladder():
    """Every ladder leg on one D2h substrate per atom (the head-to-head
    fixture's D2h twin, rebuilt here so the ladder stands alone).
    ~15 min: the semicore_d leg solves 368k determinants."""
    out = {}
    for name in sorted(_LADDER):
        spec = CAMPAIGN[name]
        term = _TERMS[name]
        mol = gto.M(atom=f"{spec.symbol} 0 0 0", basis=spec.basis,
                    spin=spec.spin, symmetry="D2h", verbose=0)
        mf = scf.ROHF(mol)
        mf.verbose = 0
        mf.conv_tol = 1e-12
        mf.kernel()
        assert mf.converged, name
        for win, fold, nelec, tag in _LADDER[name]:
            out[(name, tag)] = _ladder_leg(
                mol, mf, spec.spin, win, fold, nelec,
                term["L"], term["parity"],
                anchor=(name == "Ti" and tag == "semicore"))
    return out


def test_cas_ladder_terms(cas_ladder):
    # every m-complete window hosts the term EXACTLY: all 2L+1
    # branch components read L(L+1) and S(S+1) to machine precision and
    # sit degenerate (the uniform placement + m-scalar fold; observed
    # ptp <= 3.3e-13); the planned {3d,4d} leg is exempt -- its
    # reading is pinned separately
    for (name, tag), leg in cas_ladder.items():
        if tag == "planned":
            continue
        term = _TERMS[name]
        assert len(leg["l2"]) == 2 * term["L"] + 1, (name, tag)
        for v in leg["l2"]:
            assert abs(v - term["L"] * (term["L"] + 1)) < 1e-8, (name, tag)
        for v in leg["s2"]:
            assert abs(v - term["s2"]) < 1e-8, (name, tag)
        assert leg["ptp"] < 1e-9, (name, tag)


def test_cas_ladder_ordering(cas_ladder):
    # the Ti ordering, margins observed (fresh 2026-10-31):
    # semicore_d -106.0 < mshell -61.8 < semicore_4s -59.8 <
    # semicore -58.9 < valence -44.4 < declared 0 < planned +177.5 mEh
    e = {tag: cas_ladder[("Ti", tag)]["total"]
         for _, _, _, tag in _LADDER["Ti"]}
    order = ["semicore_d", "mshell", "semicore_4s", "semicore",
             "valence", "declared", "planned"]
    for lo, hi in zip(order, order[1:]):
        # every adjacent gap exceeds 0.5 mEh (smallest observed 0.94)
        assert e[lo] < e[hi] - 5e-4, (lo, hi)
    # the 4s pair stays in the fold: giving it back the window buys
    # <= 2 mEh on top of the semi-core window (observed 0.94)
    assert abs(e["semicore_4s"] - e["semicore"]) < 2e-3


def test_cas_ladder_planned_not_a_term_window(cas_ladder):
    # the planned {3d,4d} window as declared (fold Ar, the 4s pair in
    # the dead zone) does NOT host the 3F term: its lowest g-block
    # roots scramble D (L^2=6) and H (L^2=30) components, with quintet
    # admixture -- and the reading sits far above the declared window
    leg = cas_ladder[("Ti", "planned")]
    l2 = leg["l2"]
    assert any(abs(v - 12.0) > 0.5 for v in l2), l2
    assert leg["total"] > cas_ladder[("Ti", "declared")]["total"] + 0.15


def test_cas_ladder_semicorescaling(cas_ladder):
    # the d-count scaling, MEASURED: the m-complete 3p semi-core window
    # pays the SAME ~60 mEh at d1 (Sc, observed -59.7) as at d2 (Ti,
    # observed -58.9) -- the semi-core correlation is d-count-
    # independent; what was d2-specific (the swap mechanism pinned
    # above) was only the TORN route to a piece of it
    gains = {}
    for name in ("Sc", "Ti"):
        gains[name] = (cas_ladder[(name, "semicore")]["total"]
                       - cas_ladder[(name, "declared")]["total"])
        assert -0.065 < gains[name] < -0.055, (name, gains[name])
    assert abs(gains["Ti"] - gains["Sc"]) < 5e-3, gains


def test_cas_ladder_anchor(cas_ladder):
    # the direct block pattern's conventions: our block-restricted H
    # reproduces the pinned solver's component eigenvalues (observed
    # 4.3e-14) -- the c_core fold path certified against PySCF
    assert cas_ladder[("Ti", "semicore")]["anchor"] < 1e-10
