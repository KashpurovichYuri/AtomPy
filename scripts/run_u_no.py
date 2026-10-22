"""U NO-window pilot (ticket 22, methodology D): compress the 2f+2d
CASCI through its term-averaged density.

The 2f+2d window (24 columns, 10,626 dets, exact eigh) is the lightest
window holding BOTH configuration blocks (f2d2 ground, f1d3 excited).
Its 5L ground-term roots are averaged into the 1-RDM -- a
rotation-invariant root set, hence an m-scalar density, hence NOs in
complete occupation-degenerate m-sets (the ticket-22 claim, measured
here on real data). Occupation-GROUP windows then re-solve the sector
at fractional determinant cost; every NO window carries the
one-particle closure certificate. Under test: 5L quality at a fraction
of the cost, and the f2d2/f1d3 ordering preserved.

Bookkeeping: E_core is computed directly from the pur core (raw
mf field) and VALIDATED against the campaign total; all rows share it,
so energy DIFFERENCES are exact regardless.

Usage: python scripts/run_u_no.py
Threads: OPENBLAS_NUM_THREADS=128 OMP_NUM_THREADS=128 PYSCF_NUM_THREADS=128.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_u import CFG, CHK, NCORE, REF, build_window
from atomterm import drivers
from atomterm.natural import (
    ci_rdm1,
    natural_orbitals,
    occupation_groups,
    select_no_window,
)
from atomterm.operators import (
    build_operators,
    population_operator,
    window_closure,
)

BANDS = {"f": (0.5, 2.0), "d": (0.1, 0.5)}   # the campaign 2f2d tier
CAMPAIGN_TOTAL = -49.894287                  # tier2f2d_ms4 5L ground (log)
NLOW = 48                                    # roots kept for the average


def closed_shell_core_energy(mf, c_core):
    """E of the doubly-occupied pur core in the raw mf field."""
    from pyscf import ao2mo

    hc = c_core.T @ mf.get_hcore() @ c_core
    n = c_core.shape[1]
    eri = ao2mo.restore(1, ao2mo.kernel(mf.mol, c_core), n)
    e = 2.0 * float(np.trace(hc))
    for i in range(n):
        e += float(np.sum(2.0 * np.diag(eri[i, i])
                          - eri[i, :, :, i].diagonal()))
    return e


def main():
    mol, mf = drivers.run_rohf(CFG, CHK, ref=REF)
    closed = list(np.flatnonzero(mf.mo_occ == 2.0))
    assert len(closed) == NCORE, (len(closed), NCORE)
    C_PUR, _ = drivers.purified_core_orbitals(mol, mf, NCORE)
    e_core = closed_shell_core_energy(mf, C_PUR)
    print(f"pur core: {NCORE} orbitals, E_core = {e_core:.8f} Eh")

    f_rows, d_rows, win = build_window(mol, BANDS)
    ncas = len(win)
    print(f"window 2f+2d: {ncas} columns")
    mo = drivers.build_full_mo(mf, closed, win, c_core=C_PUR)
    c_act = mo[:, NCORE:NCORE + ncas]

    ops = build_operators(mol, mf, c_act, NCORE, 4, 0, c_core=C_PUR)
    dets = ops["dets"]
    ev, vec = np.linalg.eigh(ops["H"].real)
    e_full = ev[0] + e_core
    print(f"full window: E_5L = {e_full:.6f} Eh  (campaign "
          f"{CAMPAIGN_TOTAL:.6f}, dE = {(e_full - CAMPAIGN_TOTAL) * 1e3:+.3f} mEh"
          "  -- validates the E_core bookkeeping)")

    # root classification: the ground term = the degenerate eigenspace
    # at ev[0] (rotationally invariant H -> all M_L components exactly
    # degenerate); anchor on the ground's own <L^2>, never a hardcoded
    # term value
    exp_l2 = np.array([vec[:, r] @ ops["L2"].real @ vec[:, r]
                       for r in range(NLOW)])
    exp_s2 = np.array([vec[:, r] @ ops["S2"].real @ vec[:, r]
                       for r in range(NLOW)])
    hist = {}
    for r in range(NLOW):
        key = round(float(exp_l2[r]), 4)
        hist[key] = hist.get(key, 0) + 1
    print("low-root <L^2> histogram:", dict(sorted(hist.items())))
    l2_ground = float(exp_l2[0])
    block_l = [r for r in range(NLOW)
               if abs(exp_l2[r] - l2_ground) < 1e-6
               and ev[r] - ev[0] < 1e-6]
    block_rest = [r for r in range(NLOW) if r not in block_l]
    l = 0.5 * (math.sqrt(4.0 * l2_ground + 1.0) - 1.0)
    print(f"ground term: L = {l:.0f} (<L^2> = {l2_ground:.4f}), "
          f"{len(block_l)} roots; rest: {len(block_rest)}")

    # ---- the term-averaged density -> NOs
    w = np.zeros(NLOW)
    w[block_l] = 1.0 / len(block_l)
    gamma = ci_rdm1(dets, vec[:, :NLOW], ncas, weights=w)
    c_no, occ = natural_orbitals(c_act, gamma)
    groups = occupation_groups(occ)
    s = mol.intor_symmetric("int1e_ovlp")
    print("\nNO occupation groups (the ticket-22 claim: degenerate "
          "m-sets, sizes 7/7/5/5):")
    for gi, g in enumerate(groups):
        c_g = c_no[:, g]
        cf, cd = c_g[f_rows], c_g[d_rows]
        w_f = float(np.sum(cf.T @ s[np.ix_(f_rows, f_rows)] @ cf))
        w_d = float(np.sum(cd.T @ s[np.ix_(d_rows, d_rows)] @ cd))
        print(f"  group {gi}: n = {len(g):2d}  occ = {occ[g[0]]:.6f}"
              f"  f-weight {w_f:.3f}  d-weight {w_d:.3f}")

    # ---- group windows: certificate + re-solve
    print(f"\n{'groups':>7} {'cols':>4} {'dets':>6}   {'E_5L, Eh':>12}"
          f"  {'dE, mEh':>8}  {'<L2>':>8}  {'<S2>':>6}  {'<Nf>':>6}"
          f"  {'|[P,Lz]|':>9} {'|[P,Lperp]|':>11}")
    nf_op = None
    for nk in (1, 2, 3, 4):
        n_max = sum(len(groups[i]) for i in range(nk))
        c_win, kept = select_no_window(c_no, occ, n_max)
        k = c_win.shape[1]
        ndet = math.comb(k, 4)
        wc = window_closure(mol, c_win)
        ops_k = build_operators(mol, mf, c_win, NCORE, 4, 0,
                                c_core=C_PUR)
        ev_k, vec_k = np.linalg.eigh(ops_k["H"].real)
        l2_0 = float(vec_k[:, 0] @ ops_k["L2"].real @ vec_k[:, 0])
        s2_0 = float(vec_k[:, 0] @ ops_k["S2"].real @ vec_k[:, 0])
        nf = population_operator(mol, c_win, f_rows, 4, 0)
        nf0 = float(vec_k[:, 0] @ nf.real @ vec_k[:, 0])
        print(f"{nk:>7} {k:>4} {ndet:>6}   {ev_k[0] + e_core:>12.6f}"
              f"  {(ev_k[0] - ev[0]) * 1e3:>+8.2f}  {l2_0:>8.4f}"
              f"  {s2_0:>6.3f}  {nf0:>6.3f}"
              f"  {wc['lz']:>9.2e} {wc['lperp']:>11.2e}")

    # ---- ordering check on the 19-column window (3 groups):
    # the low spectrum must keep f2d2 below f1d3 (campaign: -49.8943 /
    # -49.8763)
    c_win, _ = select_no_window(c_no, occ,
                                sum(len(groups[i]) for i in range(3)))
    ops_k = build_operators(mol, mf, c_win, NCORE, 4, 0, c_core=C_PUR)
    nf = population_operator(mol, c_win, f_rows, 4, 0)
    ev_k, vec_k = np.linalg.eigh(ops_k["H"].real)
    print("\nordering check (19 columns): lowest 20 roots "
          "(f2d2: <Nf> = 2, f1d3: <Nf> = 1)")
    for r in range(20):
        l2_r = float(vec_k[:, r] @ ops_k["L2"].real @ vec_k[:, r])
        nf_r = float(vec_k[:, r] @ nf.real @ vec_k[:, r])
        print(f"  root {r}: E = {ev_k[r] + e_core:.6f}  "
              f"<L2> = {l2_r:8.4f}  <Nf> = {nf_r:.4f}")


if __name__ == "__main__":
    main()
