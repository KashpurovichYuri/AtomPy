"""The symmetry lane on the open d shells -- the stress case, and
where the two routes genuinely diverge.

The multiplicity-2 branching goes live through the runner (Sc's 2D:
the Ag block hosts TWO components; Ti's 3F: each B_g hosts two), and
the anchor holds at machine precision as on the p shells.

The head-to-head, though, splits by atom.  On Sc (d1) the routes
agree as on C/Si.  On Ti (d2) the free route-A SA escapes the
declared {3d,4s} window into the recorded swap basin (4s out, a full
3p shell in: lower ensemble, certificates O(1), multiplet split at
mEh) while its roots still READ as the term -- membership alone is
not the gate, the multiplet degeneracy is.  Route B's D2h-blocked
rotations forbid the escape (the swap is cross-irrep: 4s is Ag, 3p
is u): it stays in the declared identity, components degenerate to
4e-7, certificates at 1e-5 -- the guaranteed component set at the
price of ~5 mEh of ensemble energy.  The symmetry-vs-variation
trade, measured.

Windows, sectors, and weights derive from the CAMPAIGN declaration;
the ground-term table (L, 2S, parity per atom) is test data, the
same kind of row test_atoms carries."""

import numpy as np
import pytest
from pyscf import gto, scf

from atomterm.casscf import run_irrep_casci, run_irrep_sa, run_sa_casscf
from atomterm.cases import CAMPAIGN
from atomterm.irreps import IRREP_IDS, branch_d2h, det_irreps, orb_irreps
from atomterm.operators import (
    build_operators,
    ci_in_detbasis_order,
    window_closure,
)
from atomterm.system import sector_na_nb
from atomterm.windows import (
    build_full_mo,
    occupied_core_idx,
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

# route A root-set membership threshold: within-term mixes read
# <L^2> = L(L+1) exactly, neighbouring terms differ by >= 2L+2, and
# route A's optimizer-tolerance window defect shifts <L^2> at ~1e-3
_MEMBER_TOL = 0.5
# the multiplet-splitting gate between the two branches: a formed
# term splits at <= 1e-6 (observed Sc 3.1e-7, route B 3.8e-7); the
# escaped/tilted draw splits at mEh (observed Ti 3.1e-3)
_SPLIT_TOL = 1e-5

# observed bands (fresh draws, 2026-10-31): route A/B ensembles agree
# to 1.1e-8 on the symmetric draw; route B closes at 6.9e-6..1.4e-5
# (both atoms); the symmetric route A window carries an f-direction
# drift the energy cannot feel (Sc: certificates 5.2e-3 with
# ensembles agreeing to 1e-8) -- hence the generous certificate band
_ENS_BAND = 1e-6
_CERT_A = 2e-2
_CERT_B = 1e-4


def _root_vectors(ci):
    """A solver's return as a list of root vectors (nroots = 1 returns
    a bare array, nroots > 1 a stacked one)."""
    ci = np.asarray(ci)
    return [ci] if ci.ndim == 2 else list(ci)


@pytest.fixture(scope="module", params=sorted(_TERMS))
def tm_lane(request):
    """Both substrates of one reference (C1 route A, D2h twin route
    B) with the window, orbsyms, and our matrices -- everything but
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
    mol1, mol2 = mols
    mf1, mf2 = mfs
    c_act = np.eye(mol2.nao_nr())[
        :, select_valence_ao(mol2, list(spec.shell_patterns))]
    orbsyms = orb_irreps(mol2, c_act)
    ops = build_operators(mol2, mf2, c_act, spec.ncore, na, nb)
    return name, spec, term, mol1, mf1, mol2, mf2, c_act, orbsyms, ops, na, nb


def test_tm_runner_anchor(tm_lane):
    # the pinned component roots == the LOWEST eigenvalues of our H on
    # the same blocks (the multiplicity-2 blocks live: the Ag block
    # holds two 2D components on Sc, each B_g two 3F components on
    # Ti); every root reads the term exactly by OUR meters -- the
    # runner's documented "lowest per block = the term" assumption,
    # verified, not assumed; the components are degenerate on the
    # m-scalar fold
    name, spec, term, _, _, mol2, mf2, c_act, orbsyms, ops, na, nb = tm_lane
    l2 = term["L"] * (term["L"] + 1)
    roots = run_irrep_casci(mol2, mf2, c_act, spec.ncore, na, nb,
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
    # THE STRESS.  Route A (free C1 SA) HOPES its lowest 2L+1 roots
    # are the term's M_L set; route B GUARANTEES them.  Measured on
    # the open d shells, the honest gate is the multiplet DEGENERACY,
    # not membership alone: on the recorded Ti draw all 7 route-A
    # roots still read <L^2> = 12 while the free optimizer had
    # escaped the declared window (the swap basin: 4s out, a full 3p
    # shell in -- l-weights p1.0 d5.0, certificates O(1), multiplet
    # split 3.1e-3; the production tables record the same basin).
    # Route B's D2h-blocked rotations cannot follow: the swap is
    # cross-irrep (4s is Ag, 3p is u) -- it stays in the declared
    # s+d identity, degenerate and closed, ~5 mEh above.
    name, spec, term, mol1, mf1, mol2, mf2, c_act, _, _, na, nb = tm_lane
    ncomp = 2 * term["L"] + 1
    l2 = term["L"] * (term["L"] + 1)
    ncas = c_act.shape[1]

    route_b = run_irrep_sa(mol2, mf2, c_act, spec.ncore, na, nb,
                           term["L"], term["parity"])
    assert route_b["converged"], name
    mo1 = build_full_mo(
        mf1, occupied_core_idx(mf1, spec.ncore),
        np.eye(mol1.nao_nr())[
            :, select_valence_ao(mol1, list(spec.shell_patterns))])
    route_a = run_sa_casscf(mf1, ncas, (na, nb), mo1,
                            list(spec.sa_weights),
                            str(tmp_path / f"{name.lower()}_sa.chk"))
    ens_a = float(np.mean(route_a.energies))
    active = slice(spec.ncore, spec.ncore + ncas)

    # route B's component set and l-identity are by construction (one
    # pinned solver per branch irrep; g-irreps only -- a p admixture
    # is cross-irrep, unreachable); its degeneracy, certificates, and
    # term purity are observed bands
    e_b = np.concatenate([e for e, _ in route_b["roots"].values()])
    assert np.ptp(e_b) < 1e-5, name  # observed 3.8e-7
    labs, w = shell_weights(mol1, route_b["mo"][:, active],
                            mf1.get_ovlp())
    assert dict(zip(labs, w.sum(axis=0))).get("p", 0.0) < 1e-6, name
    cert_b = window_closure(mol1, route_b["mo"][:, active])
    assert cert_b["lz"] < _CERT_B and cert_b["lperp"] < _CERT_B, name
    ops_b = build_operators(mol2, mf2, route_b["mo"][:, active],
                            spec.ncore, na, nb)
    for irrep, (_, ci) in route_b["roots"].items():
        for k, vec in enumerate(_root_vectors(ci)):
            wl = ci_in_detbasis_order(vec, ncas, na, nb)
            assert abs(float(np.real(wl @ ops_b["L2"] @ wl)) - l2) < 1e-6, \
                (name, irrep, k)  # observed ~1.1e-10

    # route A's draw, measured: membership + a formed (degenerate)
    # multiplet.  Both hold -> the routes solved the same problem and
    # the ensembles must agree; otherwise the divergence IS the
    # finding, documented -- never forced into one shape
    ops_a = build_operators(mol1, mf1, route_a.mo_coeff[:, active],
                            spec.ncore, na, nb)
    l2_a = []
    for ci in route_a.ci:
        wl = ci_in_detbasis_order(ci, ncas, na, nb)
        l2_a.append(float(np.real(wl @ ops_a["L2"] @ wl)))
    member = sum(abs(v - l2) < _MEMBER_TOL for v in l2_a) == ncomp
    formed = np.ptp(route_a.energies) < _SPLIT_TOL
    cert_a = window_closure(mol1, route_a.mo_coeff[:, active])
    if member and formed:
        assert abs(ens_a - route_b["e_tot"]) < _ENS_BAND, name
        assert cert_a["lz"] < _CERT_A and cert_a["lperp"] < _CERT_A, name
    else:
        # whichever divergence shape the draw took -- the escape
        # (members all present, multiplet split: the recorded Ti swap
        # basin) or the interleave (a lower term's root displacing a
        # component) -- route A's ensemble sits BELOW route B's: by
        # free variation in the first shape, by displacement in the
        # second (the intruder enters below the term root it pushes
        # out).  Route A's certificate is measured, not pinned: the
        # divergent window class is basin chaos, and optimizer output
        # is never pinned where chaos is the object
        assert ens_a < route_b["e_tot"], name
