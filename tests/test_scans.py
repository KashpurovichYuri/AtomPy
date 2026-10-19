"""Tests of atomterm.scans: frozen-orbital scans (ticket 16, Phase 1).

Everything here is optimizer-free, so the numbers are DETERMINISTIC given
the ROHF reference (no basin chaos): the pinned structure is the physics
of the radial direction, not an optimizer's opinion.
"""

import numpy as np
import pytest

from atomterm.analysis import analyze
from atomterm.drivers import run_rohf, select_shell_ao
from atomterm.operators import build_operators
from atomterm.scans import member_tilt_map, radial_scan, scan_state, scan_window
from atomterm.system import AtomConfig


@pytest.fixture(scope="session")
def si_rohf(tmp_path_factory):
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    return run_rohf(cfg, str(tmp_path_factory.mktemp("scan_si") / "rohf.chk"))


def test_scan_windows_are_s_orthonormal(carbon):
    """The explicit cos/sin construction yields valid orthonormal columns
    for any angles (different members are exactly S-orthogonal)."""
    _, mol, mf, _ = carbon
    fam_a, fam_b = select_shell_ao(mol, "2p"), select_shell_ao(mol, "3p")
    s = mf.get_ovlp()
    for thetas in ([0.0] * 3, [0.6] * 3, [0.7, 0.5, 0.6]):
        c_act = scan_window(mol, mf, thetas, fam_a, fam_b)
        err = np.abs(c_act.conj().T @ s @ c_act - np.eye(3)).max()
        assert err < 1e-13, thetas


def test_radial_scan_carbon_full_gain_at_exact_symmetry(carbon):
    """Fact 18 (a): the ONE-PARAMETER m-symmetric radial rotation recovers
    the entire CASSCF correlation gain (theta* ~ 0.61 rad -> ~ -0.2508 Eh,
    vs the free/radial CASSCF gain of fact 16) with [H, Lz] = 0 at EVERY
    point of the curve -- the gain needs NO symmetry breaking. At theta=0
    the scan reproduces the casci_avg fold to machine precision."""
    _, mol, mf, _ = carbon
    fam_a, fam_b = select_shell_ao(mol, "2p"), select_shell_ao(mol, "3p")
    ncore, na, nb = 2, 2, 0

    p0 = scan_state(mol, mf, scan_window(mol, mf, [0.0] * 3, fam_a, fam_b),
                    ncore, na, nb)
    ops = build_operators(mol, mf, np.eye(mol.nao_nr())[:, fam_a],
                          ncore, na, nb)
    assert p0.e_active == pytest.approx(
        float(np.linalg.eigvalsh(ops["H"])[0]), abs=1e-12)
    assert analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], "x").comm_norm_h_lz < 1e-13

    pts = radial_scan(mol, mf, ncore, na, nb, fam_a, fam_b,
                      np.linspace(0.4, 0.8, 9))
    best = min(pts, key=lambda p: p.e_active)
    gain = best.e_active - p0.e_active
    assert -0.26 < gain < -0.24, gain          # the full CASSCF gain
    assert best.comm_h_lz < 1e-10               # exact symmetry throughout
    assert all(p.e_total_minus_one_two < 1e-12 for p in pts)


def test_member_tilts_exactly_flat_via_ml_dodge(carbon):
    """Fact 18 (b): tilting ONE member by any amount costs EXACTLY zero
    energy -- the 3P multiplet's M_L components each avoid one member, so
    the CI finds an untouched component (rotational degeneracy = the flat
    direction of fact 16, here measured without an optimizer). Tilting
    two members costs real energy. [H, Lz] is blind to the m=0 member
    (pz: Lz pz = 0) and to equal tilts of a full +-m pair."""
    _, mol, mf, _ = carbon
    fam_a, fam_b = select_shell_ao(mol, "2p"), select_shell_ao(mol, "3p")
    ncore, na, nb = 2, 2, 0
    theta = 0.61

    tilts = member_tilt_map(mol, mf, ncore, na, nb, fam_a, fam_b, theta, 0.10)
    e_ref = scan_state(mol, mf, scan_window(mol, mf, [theta] * 3, fam_a, fam_b),
                       ncore, na, nb).e_active
    for k in ("tilt 0", "tilt 1", "tilt 2"):
        assert tilts[k].e_active == pytest.approx(e_ref, abs=1e-12), k
    for k in ("tilt 0,1", "tilt 0,2", "tilt 1,2"):
        assert tilts[k].e_active > e_ref + 1e-4, k   # real cost (µEh..mEh)
    # the Lz blind spots: m=0 member alone, and the +-m pair tilted equally
    assert tilts["tilt 2"].comm_h_lz < 1e-10          # member 2 = pz (m=0)
    assert tilts["tilt 0,1"].comm_h_lz < 1e-10        # equal +-pair tilts
    assert tilts["tilt 0"].comm_h_lz > 1e-5           # unequal tilts: seen


def test_radial_scan_silicon_matches_casscf_gain(si_rohf):
    """Fact 18 (a) on Si: same structure as carbon (theta* ~ 0.85,
    full CASSCF gain ~ -0.127 Eh at exact symmetry; single-member tilts
    exactly flat, pairs cost real energy)."""
    mol, mf = si_rohf
    fam_a, fam_b = select_shell_ao(mol, "3p"), select_shell_ao(mol, "4p")
    ncore, na, nb = 6, 2, 0

    p0 = scan_state(mol, mf, scan_window(mol, mf, [0.0] * 3, fam_a, fam_b),
                    ncore, na, nb)
    pts = radial_scan(mol, mf, ncore, na, nb, fam_a, fam_b,
                      np.linspace(0.7, 1.0, 7))
    best = min(pts, key=lambda p: p.e_active)
    gain = best.e_active - p0.e_active
    assert -0.13 < gain < -0.12, gain
    assert best.comm_h_lz < 1e-10

    tilts = member_tilt_map(mol, mf, ncore, na, nb, fam_a, fam_b,
                            best.theta, 0.10)
    assert tilts["tilt 0"].e_active == pytest.approx(best.e_active, abs=1e-12)
    assert tilts["tilt 0,1"].e_active > best.e_active + 1e-5


def test_sc_radial_structure_inverts(sc_rohf, tmp_path):
    """Fact 18 (c) -- Sc is the opposite of C/Si: the m-SYMMETRIC d-shell
    radial rotation gains NOTHING (the curve rises: radial diffusion of a
    one-electron d shell only loses nuclear attraction), while MEMBER-
    ASYMMETRIC tilts DESCEND (~ -0.3..-0.9 mEh at delta = 0.1): a genuine
    instability along the degenerate 2D manifold (the Jahn-Teller-like
    drive behind Sc's broken single-root solutions). The 4s column is
    held fixed (4s -> 5s is harmful and would mask the d structure)."""
    cfg, mol, mf, _ = sc_rohf
    d3, d4 = select_shell_ao(mol, "3d"), select_shell_ao(mol, "4d")
    c4s = np.eye(mol.nao_nr())[:, select_shell_ao(mol, "4s")].reshape(-1, 1)
    ncore, na, nb = 9, 2, 1

    def state(thetas):
        c_act = np.hstack([scan_window(mol, mf, thetas, d3, d4), c4s])
        return scan_state(mol, mf, c_act, ncore, na, nb)

    p0 = state([0.0] * 5)
    p_up = state([0.3] * 5)
    assert (p_up.e_active - p0.e_active) > 0.02          # symmetric: harmful
    assert p_up.comm_h_lz < 1e-10                        # but exactly closed

    tilts = member_tilt_map(mol, mf, ncore, na, nb, d3, d4, 0.0, 0.10,
                            extra_cols=c4s)
    e_t1 = tilts["tilt 0"].e_active
    assert -8e-4 < e_t1 - p0.e_active < -5e-5            # asymmetric descends
    assert tilts["tilt 0"].comm_h_lz > 0.1               # and breaks L
    e_t3 = tilts["tilt 0,1,2"].e_active
    assert e_t3 < e_t1 < p0.e_active                     # more tilt, lower E


@pytest.fixture(scope="session")
def ti_scan_ref(tmp_path_factory):
    """Ti reference for the frozen scans: the character-forced, basin-gated
    ROHF (conftest TI_REF) -- deterministic closure measurements."""
    from conftest import TI_REF

    cfg = AtomConfig("Ti", "ccpvdz", spin=2)
    mol, mf = run_rohf(cfg, str(tmp_path_factory.mktemp("scan_ti") / "rohf.chk"),
                       ref=TI_REF)
    return mol, mf


def test_ti_d2_tilts_descend_radial_rises(ti_scan_ref):
    """Fact 21 (ticket 12): the d2 3F anatomy inverts the p2 story the
    SAME way Sc 3d1 2D does, but ~2x weaker:

    - SINGLE-member tilts DESCEND (measured -163 uEh at delta = 0.1 vs
      Sc's ~-300 uEh) with <L^2> ~ 12 barely moved -- a Jahn-Teller-like
      instability of the orbitally degenerate 3F term, so free single-root
      optimizers wander and sa_radial_pur stops at its flat-gauge floor
      (no restoring force tightens the member gauge);
    - the m-SYMMETRIC radial rotation of the d shell RISES steeply
      (+28 mEh at theta = 0.05): radial diffusion of the contracted 3d
      only loses, exactly like Sc (and unlike C/Si, whose p-shell radial
      breathing recovers the full CASSCF gain);
    - the {0,1,3,4} tilt pattern (everything but dz^2 = equal tilts of
      both +-m pairs) stays Lz-BLIND (|[H, Lz]| ~ 1e-13) while [H, L^2]
      fires -- the fact-18c blind-spot structure on a 5-member shell.

    Deterministic: no optimizer anywhere in the loop."""
    mol, mf = ti_scan_ref
    d3, d4 = select_shell_ao(mol, "3d"), select_shell_ao(mol, "4d")
    c4s = np.eye(mol.nao_nr())[:, select_shell_ao(mol, "4s")].reshape(-1, 1)
    ncore, na, nb = 9, 3, 1

    def state(thetas):
        c_act = np.hstack([scan_window(mol, mf, thetas, d3, d4), c4s])
        return scan_state(mol, mf, c_act, ncore, na, nb, with_l2=True)

    p0 = state([0.0] * 5)
    assert p0.ground_l2 == pytest.approx(12.0, abs=1e-8)  # 3F ground

    # m-symmetric radial rotation: harmful (the Sc class, not C/Si;
    # measured +9.8 mEh at theta = 0.05, rising steeply)
    p_rad = state([0.05] * 5)
    assert p0.e_active + 0.005 < p_rad.e_active < p0.e_active + 0.02
    assert p_rad.comm_h_lz < 1e-10          # but exactly closed
    assert p_rad.ground_l2 == pytest.approx(12.0, abs=1e-6)

    # single-member tilt: DESCENDS (JT class), breaks L
    p_t1 = state([0.10] * 1 + [0.0] * 4)
    assert -8e-4 < p_t1.e_active - p0.e_active < -5e-5
    assert p_t1.comm_h_lz > 0.1
    assert 11.9 < p_t1.ground_l2 < 12.0     # term identity nearly kept

    # two members: roughly twice the descent
    p_t2 = state([0.10, 0.10, 0.0, 0.0, 0.0])
    assert p_t2.e_active < p_t1.e_active < p0.e_active

    # the Lz-blind pattern: equal tilts of both +-m pairs (no dz^2)
    p_bl = state([0.10, 0.10, 0.0, 0.10, 0.10])
    assert p_bl.comm_h_lz < 1e-10           # Lz cannot see it
    assert p_bl.comm_h_l2 > 0.1             # L^2 fires
