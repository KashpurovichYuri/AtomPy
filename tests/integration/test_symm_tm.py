"""The symmetry lane on the open d shells -- the stress case, and
where the two routes genuinely diverge.

The multiplicity-2 branching goes live through the runner (Sc's 2D:
the Ag block hosts TWO components; Ti's 3F: each B_g hosts two), and
the anchor holds at machine precision as on the p shells.

The head-to-head, though, splits by atom.  On Sc (d1) the routes
agree as on C/Si.  On Ti (d2) the free SA escapes the
declared {3d,4s} window into the recorded swap basin (the 4s column
replaced by a single TORN p function -- one m, a pair relocated
into it: lower ensemble, certificates O(1), multiplet split at
mEh) while its roots still READ as the term -- membership alone is
not the gate, the multiplet degeneracy is.  The pinned route's D2h-blocked
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
# (both atoms); the symmetric the free route window carries an f-direction
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
    """Both substrates of one reference (C1 the free route, D2h twin route
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
    # are the term's M_L set; the pinned route GUARANTEES them.  Measured on
    # the open d shells, the honest gate is the multiplet DEGENERACY,
    # not membership alone: on the recorded Ti draw all 7 free-route
    # roots still read <L^2> = 12 while the free optimizer had
    # escaped the declared window (the swap basin: 4s out, a full 3p
    # shell in -- l-weights p1.0 d5.0, certificates O(1), multiplet
    # split 3.1e-3; the production tables record the same basin).
    # The pinned route's D2h-blocked rotations cannot follow: the swap is
    # cross-irrep (4s is Ag, 3p is u) -- it stays in the declared
    # s+d identity, degenerate and closed, ~5 mEh above.
    name, spec, term, mol_c1, mf_c1, mol_d2h, mf_d2h, c_act, _, _, na, nb = tm_lane
    ncomp = 2 * term["L"] + 1
    l2 = term["L"] * (term["L"] + 1)
    ncas = c_act.shape[1]

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

    # the pinned route's component set and l-identity are by construction (one
    # pinned solver per branch irrep; g-irreps only -- a p admixture
    # is cross-irrep, unreachable); its degeneracy, certificates, and
    # term purity are observed bands
    e_pinned = np.concatenate([e for e, _ in pinned_sa["roots"].values()])
    assert np.ptp(e_pinned) < 1e-5, name  # observed 3.8e-7
    labs, w = shell_weights(mol_c1, pinned_sa["mo"][:, active],
                            mf_c1.get_ovlp())
    assert dict(zip(labs, w.sum(axis=0))).get("p", 0.0) < 1e-6, name
    cert_pinned = window_closure(mol_c1, pinned_sa["mo"][:, active])
    assert cert_pinned["lz"] < _CERT_PINNED and cert_pinned["lperp"] < _CERT_PINNED, name
    ops_pinned = build_operators(mol_d2h, mf_d2h, pinned_sa["mo"][:, active],
                            spec.ncore, na, nb)
    for irrep, (_, ci) in pinned_sa["roots"].items():
        for k, vec in enumerate(_root_vectors(ci)):
            wl = ci_in_detbasis_order(vec, ncas, na, nb)
            assert abs(float(np.real(wl @ ops_pinned["L2"] @ wl)) - l2) < 1e-6, \
                (name, irrep, k)  # observed ~1.1e-10

    # the free route's draw, measured: membership + a formed (degenerate)
    # multiplet.  Both hold -> the routes solved the same problem and
    # the ensembles must agree; otherwise the divergence IS the
    # finding, documented -- never forced into one shape
    ops_free = build_operators(mol_c1, mf_c1, free_sa.mo_coeff[:, active],
                            spec.ncore, na, nb)
    l2_free = []
    for ci in free_sa.ci:
        wl = ci_in_detbasis_order(ci, ncas, na, nb)
        l2_free.append(float(np.real(wl @ ops_free["L2"] @ wl)))
    member = sum(abs(v - l2) < _MEMBER_TOL for v in l2_free) == ncomp
    formed = np.ptp(free_sa.energies) < _SPLIT_TOL
    cert_free = window_closure(mol_c1, free_sa.mo_coeff[:, active])
    if member and formed:
        assert abs(ens_free - pinned_sa["e_tot"]) < _ENS_BAND, name
        assert cert_free["lz"] < _CERT_FREE and cert_free["lperp"] < _CERT_FREE, name
    else:
        # whichever divergence shape the draw took -- the escape
        # (members all present, multiplet split: the recorded Ti swap
        # basin) or the interleave (a lower term's root displacing a
        # component) -- the free route's ensemble sits BELOW the pinned route's: by
        # free variation in the first shape, by displacement in the
        # second (the intruder enters below the term root it pushes
        # out).  The free route's certificate is measured, not pinned: the
        # divergent window class is basin chaos, and optimizer output
        # is never pinned where chaos is the object
        assert ens_free < pinned_sa["e_tot"], name
