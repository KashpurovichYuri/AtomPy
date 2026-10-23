"""U ticket 23: CASSCF/SA started from NO windows -- relaxation vs growth.

The question (ticket 23): can orbital RELAXATION recover what window
GROWTH delivers -- and does the NO start (exactly closed, term-averaged,
compressed) tame the free-CASSCF escape basins? The window axis is the
measurement: every rung runs on the NO-group ladder nk = 1..4 and is
read against the frozen NO-CASCI curve (results/u_14e/u_run_no.log):

    nk=1   5 cols      5 dets   -49.660672   <L2>=6    (d4, 5D block)
    nk=2  12 cols    495 dets   -49.889731   <L2>=72   (f2d2, 5L block)
    nk=3  17 cols   2380 dets   -49.890855   <L2>=72
    nk=4  24 cols  10626 dets   -49.894287   (= campaign full window)

If CASSCF(nk) lands near frozen(nk'), relaxation substituted for growth;
if the CASSCF curve is flat above frozen(4), growth is irreducible.
Matched-cost shell pairs: 1f-tier shell window (7f+5d = 12 cols, 495
dets -- the same cost as NO-12) and the 2f2d shell window (24 cols).

Stages (run one at a time, inspect between; cheap first):
    prep    ROHF (cached chk) + pur core + full-window CASCI -> term
            average -> NOs -> occupation groups; npz-cached.
    casscf  free single-root CASSCF: NO starts nk=1,2,3, then the
            matched-cost pairs (shell-12 vs NO-12), then the heavy
            24-col pair (shell vs NO start, same span -- start-gauge
            check). Per run: ortho_err, E, <L2>, <S2>, ||[H,Lz]|| (own
            fold), span certificate, per-column f/d weights.
    sa      17-root SA over the measured 5L M_L block (weights uniform)
            at nk=2 and nk=3; per-root table from the built-in CASCI
            pass. The block is MEASURED per window (nk=1 is a 5-root 5D
            block -- anchor on the ground's own <L2>, never hardcode).
    relax   core-relaxation rungs on the NO substrate (span frozen),
            nk=1,2,3: the honest symmetric CASSCF rungs.
    six13   CAS(6,13) 7s2-as-valence probe (16,731 dets, exact eigh),
            intent header printed first (the manifest precursor).

Usage: python scripts/run_u_sa.py <stage> [<stage> ...]
Threads: OPENBLAS_NUM_THREADS=64 OMP_NUM_THREADS=64 PYSCF_NUM_THREADS=64.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_u import CFG, CHK, NCORE, REF, build_window
from run_u_no import BANDS, CAMPAIGN_TOTAL, closed_shell_core_energy
from atomterm import drivers
from atomterm.detbasis import commutator, frobenius_norm
from atomterm.natural import ci_rdm1, natural_orbitals, occupation_groups
from atomterm.operators import (
    build_operators,
    window_closure,
)

NLOW = 48
CACHE = Path("results/u_14e/u_no_cache.npz")
ANCHOR = {5: -49.660672, 12: -49.889731, 17: -49.890855, 24: -49.894287}
SHELL12_BANDS = {"f": (1.0, 2.0), "d": (0.1, 0.25)}   # the campaign 1f tier
NA, NB = 4, 0   # maximal M_S = 2 sector; the 5L ground lives here
ORTHOGUARD = 1e-8


# ---------------------------------------------------------------- setup
def setup():
    mol, mf = drivers.run_rohf(CFG, CHK, ref=REF)
    closed = list(np.flatnonzero(mf.mo_occ == 2.0))
    assert len(closed) == NCORE, (len(closed), NCORE)
    c_pur, _ = drivers.purified_core_orbitals(mol, mf, NCORE)
    e_core = closed_shell_core_energy(mf, c_pur)
    return mol, mf, closed, c_pur, e_core


def build_nos(mol, mf, c_pur, e_core):
    """Full-window CASCI -> term-average -> NOs; npz-cached (the eigh is
    the ~6-min piece; CASSCF stages must not pay it again)."""
    if CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        return z["c_no"], z["occ"], list(z["groups"]), float(z["e_full"])
    closed = list(np.flatnonzero(mf.mo_occ == 2.0))
    _, _, win = build_window(mol, BANDS)
    mo = drivers.build_full_mo(mf, closed, win, c_core=c_pur)
    c_act = mo[:, NCORE:NCORE + len(win)]
    ops = build_operators(mol, mf, c_act, NCORE, NA, NB, c_core=c_pur)
    ev, vec = np.linalg.eigh(ops["H"].real)
    l2 = np.array([vec[:, r] @ ops["L2"].real @ vec[:, r] for r in range(NLOW)])
    block = [r for r in range(NLOW)
             if abs(l2[r] - l2[0]) < 1e-6 and ev[r] - ev[0] < 1e-6]
    assert len(block) == 17 and abs(l2[0] - 72.0) < 1e-4, (len(block), l2[0])
    w = np.zeros(NLOW)
    w[block] = 1.0 / len(block)
    gamma = ci_rdm1(ops["dets"], vec[:, :NLOW], len(win), weights=w)
    c_no, occ = natural_orbitals(c_act, gamma)
    groups = occupation_groups(occ)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, c_no=c_no, occ=occ, groups=np.array(groups, dtype=object),
             e_full=ev[0] + e_core)
    print(f"NO cache written: groups = {[len(g) for g in groups]}, "
          f"E_full = {ev[0] + e_core:.6f}")
    return c_no, occ, list(groups), float(ev[0] + e_core)


def no_start(mf, c_pur, c_no, groups, nk):
    """Block-safe full MO matrix: [pur core | first nk NO groups | complement]."""
    nact = sum(len(groups[i]) for i in range(nk))
    cols = np.hstack([c_pur, c_no[:, :nact]])
    mo = drivers.complete_orthonormal(cols, mf.get_ovlp())
    err = float(np.abs(mo.conj().T @ mf.get_ovlp() @ mo - np.eye(mo.shape[1])).max())
    assert err < ORTHOGUARD, (nk, err)
    return mo, err


# ----------------------------------------------------------- measurement
def col_weights(mol, c_win, f_rows, d_rows):
    """Per-column f/d S-weight (l-content of optimized columns)."""
    s = mol.intor_symmetric("int1e_ovlp")
    out = []
    for k in range(c_win.shape[1]):
        c = c_win[:, k]
        wf = float(c[f_rows] @ (s @ c)[f_rows])
        wd = float(c[d_rows] @ (s @ c)[d_rows])
        out.append((wf, wd))
    return out


def report(mol, mf, c_pur, tag, cas, anchor):
    """Basin forensics for one finished single-root run."""
    mo = cas.mo_coeff
    ncas = cas.mc.ncas   # NEVER mo.shape[1] - NCORE: free CASSCF returns
    assert ncas in ANCHOR, ncas   # the FULL matrix (window + virtuals)
    ortho = float(np.abs(mo.conj().T @ mf.get_ovlp() @ mo
                         - np.eye(mo.shape[1])).max())
    c_act, c_own = mo[:, NCORE:NCORE + ncas], mo[:, :NCORE]
    ops_pur = build_operators(mol, mf, c_act, NCORE, NA, NB, c_core=c_pur)
    ops_own = build_operators(mol, mf, c_act, NCORE, NA, NB, c_core=c_own)
    # the honest labels: exact diagonalization in the converged window.
    # mc.ci can be a Davidson-straggled mixture inside degenerate
    # clusters (the run_u.py disease) -- never trust it for <L^2>.
    ev_x, vec_x = np.linalg.eigh(ops_own["H"].real)
    ci = np.asarray(cas.ci[0]).ravel()
    l2_ci = float(ci @ ops_pur["L2"].real @ ci)
    l2 = float(vec_x[:, 0] @ ops_pur["L2"].real @ vec_x[:, 0])
    s2 = float(vec_x[:, 0] @ ops_pur["S2"].real @ vec_x[:, 0])
    clz = frobenius_norm(commutator(ops_own["H"], ops_own["Lz"]))
    wc = window_closure(mol, c_act)
    f_rows, d_rows, _ = build_window(mol, BANDS)
    wts = col_weights(mol, c_act, f_rows, d_rows)
    big = sorted(range(ncas), key=lambda k: -max(wts[k]))[:6]
    desc = "  ".join(
        f"col{k}:{'f' if wts[k][0] >= wts[k][1] else 'd'}"
        f"{max(wts[k]):.2f}" for k in big)
    print(f"\n[{tag}] converged={cas.mc.converged}  ortho_err = {ortho:.2e}")
    print(f"  E_tot   = {cas.energies[0]:.6f}  "
          f"(vs frozen {ncas}-col anchor {ANCHOR[ncas]:.6f}: "
          f"{(cas.energies[0] - ANCHOR[ncas]) * 1e3:+.2f} mEh,"
          f"  vs campaign {CAMPAIGN_TOTAL:.6f}: "
          f"{(cas.energies[0] - CAMPAIGN_TOTAL) * 1e3:+.2f} mEh)")
    print(f"  <L2> = {l2:.4f} (exact eigh ground; mc.ci: {l2_ci:.4f})   "
          f"<S2> = {s2:.4f}   ||[H,Lz]|| (own fold) = {clz:.2e}")
    print(f"  window span: |[P,Lz]| = {wc['lz']:.2e}  |[P,Lperp]| = {wc['lperp']:.2e}")
    print(f"  top columns: {desc}")
    return {"tag": tag, "ncas": ncas, "E": cas.energies[0], "L2": l2,
            "S2": s2, "clz": clz, "lz": wc["lz"], "lperp": wc["lperp"],
            "ortho": ortho, "conv": bool(cas.mc.converged)}


def free_casscf(mol, mf, closed, c_pur, tag, ncas, mo_full, name):
    cas = drivers.constrained_casscf(mf, ncas, (NA, NB), CHK,
                                     mo_full=mo_full, name=name)
    return report(mol, mf, c_pur, tag, cas, ANCHOR[ncas])


def ground_block(mol, mf, c_pur, c_act):
    """The frozen ground-term block of a window, MEASURED: the exact
    eigenspace at ev[0], with its own <L2>. Returns (nroots, l2, weights)."""
    ops = build_operators(mol, mf, c_act, NCORE, NA, NB, c_core=c_pur)
    ev, vec = np.linalg.eigh(ops["H"].real)
    nchk = min(60, len(ev))
    l2 = np.array([vec[:, r] @ ops["L2"].real @ vec[:, r] for r in range(nchk)])
    block = [r for r in range(nchk)
             if abs(l2[r] - l2[0]) < 1e-6 and ev[r] - ev[0] < 1e-6]
    w = np.zeros(len(ev))
    w[block] = 1.0 / len(block)
    return ops, ev, vec, len(block), float(l2[0]), w


# ---------------------------------------------------------------- stages
def stage_prep():
    mol, mf, closed, c_pur, e_core = setup()
    print(f"pur core: {NCORE} cols, E_core = {e_core:.8f} Eh")
    c_no, occ, groups, e_full = build_nos(mol, mf, c_pur, e_core)
    print(f"groups: {[len(g) for g in groups]}   "
          f"occ: {[round(float(occ[g[0]]), 6) for g in groups]}")
    print(f"full-window anchor: {e_full:.6f} (campaign {CAMPAIGN_TOTAL})")


def stage_casscf():
    mol, mf, closed, c_pur, e_core = setup()
    c_no, occ, groups, _ = build_nos(mol, mf, c_pur, e_core)
    rows = []
    # (1) the NO ladder, cheap first: nk = 1, 2, 3
    for nk in (1, 2, 3):
        ncas = sum(len(groups[i]) for i in range(nk))
        mo, _ = no_start(mf, c_pur, c_no, groups, nk)
        rows.append(free_casscf(mol, mf, closed, c_pur, f"free/no-{ncas}col",
                                ncas, mo, f"casscf_no{ncas}"))
    # (2) matched-cost shell pair at 12 cols: 1f-tier shell vs NO groups
    _, _, win12 = build_window(mol, SHELL12_BANDS)
    mo = drivers.build_full_mo(mf, closed, win12, c_core=c_pur)
    rows.append(free_casscf(mol, mf, closed, c_pur, "free/shell-12col",
                            len(win12), mo, "casscf_shell12"))
    # (3) the heavy 24-col pair, same span from two starts: shell (the
    # ticket-0 run) and NO -- the start-gauge check at fixed span
    _, _, win24 = build_window(mol, BANDS)
    mo = drivers.build_full_mo(mf, closed, win24, c_core=c_pur)
    rows.append(free_casscf(mol, mf, closed, c_pur, "free/shell-24col",
                            len(win24), mo, "casscf_shell24"))
    mo, _ = no_start(mf, c_pur, c_no, groups, 4)
    rows.append(free_casscf(mol, mf, closed, c_pur, "free/no-24col",
                            24, mo, "casscf_no24"))
    print("\nbasin table (frozen curve: 5/12/17/24 cols = "
          "-49.6607 / -49.8897 / -49.8909 / -49.8943):")
    for r in rows:
        print(f"  {r['tag']:>18}: E = {r['E']:.6f}  <L2> = {r['L2']:8.4f}"
              f"  ||[H,Lz]|| = {r['clz']:.1e}  |[P,Lp]| = {r['lperp']:.1e}"
              f"  ortho = {r['ortho']:.1e}  conv = {r['conv']}")


def stage_sa(nk=2):
    mol, mf, closed, c_pur, e_core = setup()
    c_no, occ, groups, _ = build_nos(mol, mf, c_pur, e_core)
    ncas = sum(len(groups[i]) for i in range(nk))
    mo, _ = no_start(mf, c_pur, c_no, groups, nk)
    # the SA block is MEASURED in the start window (nk=1 would be a
    # 5-root 5D block; nk>=2 the 17-root 5L block)
    ops0, ev0, vec0, nroots, l2_0, w = ground_block(
        mol, mf, c_pur, mo[:, NCORE:NCORE + ncas])
    print(f"[sa/no-{ncas}col] measured ground block: {nroots} roots, "
          f"<L2> = {l2_0:.4f}, E_frozen = {ev0[0] + e_core:.6f}")
    cas = drivers.run_sa_casscf(mf, ncas, (NA, NB), mo,
                                [1.0 / nroots] * nroots, CHK)
    print(f"[sa/no-{ncas}col] ensemble E = {cas.mc.e_tot:.6f}  "
          f"converged={cas.mc.converged}")
    c_act = cas.mo_coeff[:, NCORE:NCORE + cas.mc.ncas]
    ops = build_operators(mol, mf, c_act, NCORE, NA, NB,
                          c_core=cas.mo_coeff[:, :NCORE])
    wc = window_closure(mol, c_act)
    print(f"  span: |[P,Lz]| = {wc['lz']:.2e}  |[P,Lperp]| = {wc['lperp']:.2e}")
    print(f"  root          E, Eh      <L2>      <S2>    dE vs frozen (mEh)")
    # window-identity meter: total f/d S-weight over the converged
    # columns (start: 7.000 f / 5.000 d for the nk=2 group window) --
    # drift means the optimizer traded window columns with core/virtuals
    f_rows, d_rows, _ = build_window(mol, BANDS)
    wts = col_weights(mol, c_act, f_rows, d_rows)
    print(f"  window identity: f = {sum(w[0] for w in wts):.3f}  "
          f"d = {sum(w[1] for w in wts):.3f}  (start 7.000/5.000)")
    # which material entered/left: S-weight of the converged window in
    # each l of the AO space
    s_full = mol.intor_symmetric("int1e_ovlp")
    sect = drivers.lm_sectors(mol)
    tot: dict = {}
    for (l, _m), rows in sorted(sect.items()):
        c = c_act[rows]
        tot[l] = tot.get(l, 0.0) + float(
            np.trace(c.conj().T @ s_full[np.ix_(rows, rows)] @ c))
    print("  l-content of the converged window (S-weight): "
          + "  ".join(f"l{l}={v:.3f}" for l, v in sorted(tot.items())))
    for r in range(nroots):
        ci = np.asarray(cas.ci[r]).ravel()
        l2 = float(ci @ ops["L2"].real @ ci)
        s2 = float(ci @ ops["S2"].real @ ci)
        d = (cas.energies[r] - (ev0[0] + e_core)) * 1e3
        print(f"  {r:3d}  {cas.energies[r]:.6f}  {l2:9.4f}  {s2:7.4f}  {d:+8.2f}")
    # the honest labels: exact eigh in the converged window (the roots
    # pass is Davidson -- straggles in degenerate clusters, the run_u.py
    # disease; the exact ground decides whether the SA window still
    # hosts a quantized term)
    from atomterm.detbasis import commutator as _comm, frobenius_norm as _fn
    clz = _fn(_comm(ops["H"], ops["Lz"]))
    ev_x, vec_x = np.linalg.eigh(ops["H"].real)
    l2_x = float(vec_x[:, 0] @ ops["L2"].real @ vec_x[:, 0])
    s2_x = float(vec_x[:, 0] @ ops["S2"].real @ vec_x[:, 0])
    print(f"  exact eigh ground: E = {ev_x[0]:.6f}  <L2> = {l2_x:.4f}  "
          f"<S2> = {s2_x:.4f}   ||[H,Lz]|| (SA fold) = {clz:.2e}")


def stage_relax():
    mol, mf, closed, c_pur, e_core = setup()
    c_no, occ, groups, _ = build_nos(mol, mf, c_pur, e_core)
    for nk in (1, 2, 3):
        ncas = sum(len(groups[i]) for i in range(nk))
        mo, _ = no_start(mf, c_pur, c_no, groups, nk)
        cas = drivers.run_casscf_core_relax(mf, ncas, (NA, NB), mo, CHK)
        report(mol, mf, c_pur, f"relax/no-{ncas}col", cas, ANCHOR[ncas])


def stage_six13():
    mol, mf, closed, c_pur, e_core = setup()
    c_no, occ, groups, _ = build_nos(mol, mf, c_pur, e_core)

    # intent header (the manifest precursor): declare before computing
    print("""INTENT -- CAS(6,13), sector (5,1), 16731 dets, exact eigh:
  valence    : 5f^3 6d^1 7s^2 (7s promoted out of the pur core)
  window     : NO groups {d5, f7} (ticket-22) + the 7s pur column
  reachable  : f^k d^m s^n, k+m+n = 6, k<=3, m<=2, n<=2 -- incl. f2d2s2,
               f3d1s2, f4s2; tests whether the 7s2 pair, made active,
               moves the 5L ground verdict or the f2d2/f3d1 ordering
  fold       : pur core minus the 7s column (still m-scalar: an s column
               leaves, the 6p trio stays m-complete)
  bookkeeping: E_core recomputed on the 4-column core; differences exact
""")

    # identify the 7s column: pur columns are l-pure [s, s, p, p, p];
    # the 7s = the s column whose dominant OCCUPIED-block canonical MO
    # (occupied_core_idx -- canonical order is meaningless on this
    # basis, never mo_energy[:ncore]!) has the highest energy
    s = mf.get_ovlp()
    sect = drivers.lm_sectors(mol)
    occ_idx = drivers.occupied_core_idx(mf, NCORE)
    canon = mf.mo_coeff[:, occ_idx]
    ov = np.abs(c_pur.conj().T @ s @ canon)
    dom = ov.argmax(axis=1)              # pur col -> occupied canonical idx
    eps_occ = mf.mo_energy[occ_idx]
    s_cols = []
    for k in range(NCORE):
        c = c_pur[:, k]
        w0 = sum(float(c[rows] @ (s @ c)[rows])
                 for (l, _m), rows in sect.items() if l == 0)
        if w0 > 0.99:
            s_cols.append(k)
    assert len(s_cols) == 2, (s_cols, NCORE)
    k7s = max(s_cols, key=lambda k: eps_occ[dom[k]])
    print(f"pur s columns {s_cols}; 7s = col {k7s} "
          f"(canonical eps = {eps_occ[dom[k7s]]:.4f}; "
          f"other: {eps_occ[dom[[k for k in s_cols if k != k7s][0]]]:.4f})")
    c_core4 = np.delete(c_pur, k7s, axis=1)
    # window layout: [d5 (0-4), f7 (5-11), 7s (12)] -- the 7s column LAST
    cols = np.hstack([c_core4,
                      c_no[:, :sum(len(groups[i]) for i in range(2))],
                      c_pur[:, [k7s]]])
    mo = drivers.complete_orthonormal(cols, s)
    err = float(np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max())
    assert err < ORTHOGUARD, err
    ncore, ncas = 4, 13
    e_core4 = closed_shell_core_energy(mf, c_core4)
    print(f"core4 E = {e_core4:.6f} (drop = {e_core - e_core4:.6f})")

    ops = build_operators(mol, mf, mo[:, ncore:ncore + ncas], ncore, 5, 1,
                          c_core=c_core4)
    dets = ops["dets"]
    assert len(dets) == 16731, len(dets)
    n7s = np.array([sum(1 for so_ in det if so_ % ncas == 12)
                    for det in dets], dtype=float)
    nf = np.array([sum(1 for so_ in det if 5 <= so_ % ncas < 12)
                   for det in dets], dtype=float)
    nd = np.array([sum(1 for so_ in det if so_ % ncas < 5)
                   for det in dets], dtype=float)
    assert n7s.min() >= 0 and n7s.max() <= 2
    ev, vec = np.linalg.eigh(ops["H"].real)
    e0 = ev[0] + e_core4
    nr = 12
    norm0 = float(vec[:, 0] @ vec[:, 0])
    assert abs(norm0 - 1.0) < 1e-10, norm0
    l2 = np.array([vec[:, r] @ ops["L2"].real @ vec[:, r] for r in range(nr)])
    s2 = np.array([vec[:, r] @ ops["S2"].real @ vec[:, r] for r in range(nr)])
    # diagonal-operator expectations are |c_i|^2-weighted, NEVER the
    # signed sum (the signed dot produced "negative occupations" once)
    p2 = np.abs(vec) ** 2
    n7s_e = np.array([p2[:, r] @ n7s for r in range(nr)])
    nf_e = np.array([p2[:, r] @ nf for r in range(nr)])
    nd_e = np.array([p2[:, r] @ nd for r in range(nr)])
    # bulletproof consistency: every det carries exactly 6 electrons,
    # so ANY unit vector must give <N_total> = 6 exactly
    ntot = n7s_e + nf_e + nd_e
    assert np.allclose(ntot, 6.0, atol=1e-9), ntot
    wc = window_closure(mol, mo[:, ncore:ncore + ncas])
    print(f"\nCAS(6,13): E_0 = {e0:.6f}   <L2> = {l2[0]:.4f}  <S2> = {s2[0]:.4f}"
          f"   |[P,Lz]| = {wc['lz']:.2e}  |[P,Lp]| = {wc['lperp']:.2e}")
    print(f"  vs campaign 5L (-49.894287): {(e0 - CAMPAIGN_TOTAL) * 1e3:+.2f} mEh"
          "   (both totals carry their own E_core; compare totals as printed)")
    p0 = np.abs(vec[:, 0]) ** 2
    print("  dominant configurations of the ground root:")
    for i in np.argsort(-p0)[:5]:
        print(f"    weight {p0[i]:.3f}:  Nf = {nf[i]:.0f}  Nd = {nd[i]:.0f}"
              f"  N7s = {n7s[i]:.0f}")
    print("  low roots:")
    for r in range(nr):
        print(f"   {r:2d}  E = {ev[r] + e_core4:.6f}  <L2> = {l2[r]:8.4f}"
              f"  <S2> = {s2[r]:6.3f}  <N7s> = {n7s_e[r]:5.3f}"
              f"  <Nf> = {nf_e[r]:6.3f}  <Nd> = {nd_e[r]:6.3f}")


def stage_diag():
    """One representative free run with the optimizer trace visible and a
    long cycle budget: is the free-run endpoint a PLATEAU (basin) or a
    TRUNCATED DESCENT (runaway)? Decides how the casscf-stage table reads."""
    mol, mf, closed, c_pur, e_core = setup()
    c_no, occ, groups, _ = build_nos(mol, mf, c_pur, e_core)
    ncas = sum(len(groups[i]) for i in range(2))
    mo, _ = no_start(mf, c_pur, c_no, groups, 2)
    cas = drivers.constrained_casscf(
        mf, ncas, (NA, NB), CHK, mo_full=mo, name="casscf_no12_diag",
        mc_opts={"max_cycle_macro": 100}, verbose=4)
    print(f"\n[diag] converged = {cas.mc.converged}  "
          f"E = {cas.energies[0]:.6f}")


STAGES = {"prep": stage_prep, "casscf": stage_casscf, "sa": stage_sa,
          "relax": stage_relax, "six13": stage_six13, "diag": stage_diag}

if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "sa":
        stage_sa(int(args[1]) if len(args) > 1 else 2)
    else:
        for name in args or ["prep"]:
            print(f"\n===== stage {name} =====")
            STAGES[name]()
