"""Sc/Ti NO-window pilot (ticket 22): the NO column of the matched
tables. Same construction as the U pilot (run_u_no.py): the {3d,4d,4s,
4p} window CASCI's ground-term roots -> term-averaged 1-RDM -> NOs in
occupation groups -> group windows -> energies + certificates.

Expected structure: l-pure groups of 5/5/1/3 (3d/4d/4s/4p radials);
the leading groups should reproduce the campaign CAS(k) shell windows
(one radial per l = the same span), the later groups add the
correlating radials in occupation order.

Usage: python scripts/run_tm_no.py [sc|ti]
Threads: OPENBLAS_NUM_THREADS=128 OMP_NUM_THREADS=128 PYSCF_NUM_THREADS=128.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_tm_methb import CASES, RES
from atomterm import drivers
from atomterm.natural import (
    ci_rdm1,
    natural_orbitals,
    occupation_groups,
    select_no_window,
)
from atomterm.operators import build_operators, window_closure
from atomterm.system import AtomConfig

PATTERNS = ["3d", "4d", "4s", "4p"]
NLOW = 24


def _ti_ref():
    from atomterm.drivers import ReferenceSpec
    return ReferenceSpec(
        occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
        occ_open=("3d",), occ_e_ref=-848.40648875)


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


def run_atom(name):
    kw = CASES[name]
    sym, spin = kw["symbol"], kw["spin"]
    na, nb = kw["sector"]
    ncore = kw["ncore"]
    print(f"\n{'#' * 70}\n# {sym} NO windows (ticket 22)\n{'#' * 70}")
    cfg = AtomConfig(sym, "ccpvdz", spin=spin)
    chk = str(RES / f"{sym}_ccpvdz_s{spin}_rohf.chk")
    mol, mf = drivers.run_rohf(
        cfg, chk, ref=(kw["ref_factory"]() if kw.get("ref_factory")
                       else None))
    c_pur, contam = drivers.purified_core_orbitals(mol, mf, ncore)
    print(f"pur core: max contamination {max(contam):.2e}")
    closed = list(np.flatnonzero(mf.mo_occ == 2.0))[:ncore]
    e_core = closed_shell_core_energy(mf, c_pur)

    rows = drivers.select_valence_ao(mol, PATTERNS)
    ncas = len(rows)
    mo = drivers.build_full_mo(mf, closed, rows, c_core=c_pur)
    c_act = mo[:, ncore:ncore + ncas]

    ops = build_operators(mol, mf, c_act, ncore, na, nb, c_core=c_pur)
    dets = ops["dets"]
    ev, vec = np.linalg.eigh(ops["H"].real)
    e_full = ev[0] + e_core
    print(f"window {PATTERNS}: {ncas} cols, E_ground = {e_full:.6f} Eh")

    nlow = min(NLOW, len(dets))
    exp_l2 = np.array([vec[:, r] @ ops["L2"].real @ vec[:, r]
                       for r in range(nlow)])
    l2_ground = float(exp_l2[0])
    block = [r for r in range(nlow)
             if abs(exp_l2[r] - l2_ground) < 1e-6
             and ev[r] - ev[0] < 1e-6]
    l = 0.5 * (math.sqrt(4.0 * l2_ground + 1.0) - 1.0)
    print(f"ground term: L = {l:.0f} (<L^2> = {l2_ground:.4f}), "
          f"{len(block)} roots averaged")

    w = np.zeros(nlow)
    w[block] = 1.0 / len(block)
    gamma = ci_rdm1(dets, vec[:, :nlow], ncas, weights=w)
    c_no, occ = natural_orbitals(c_act, gamma)
    groups = occupation_groups(occ)
    s = mol.intor_symmetric("int1e_ovlp")
    shells = {pat: drivers.select_valence_ao(mol, [pat])
              for pat in PATTERNS}
    print("NO occupation groups:")
    for gi, g in enumerate(groups):
        c_g = c_no[:, g]
        wts = {}
        for pat, srows in shells.items():
            sub = c_g[srows]
            wts[pat] = float(np.sum(
                sub.T @ s[np.ix_(srows, srows)] @ sub))
        tag = " ".join(f"{p}={v:.2f}" for p, v in wts.items() if v > 0.01)
        print(f"  group {gi}: n = {len(g):2d}  occ = {occ[g[0]]:.6f}"
              f"  [{tag}]")

    print(f"\n{'groups':>7} {'cols':>4} {'dets':>6}   {'E, Eh':>12}"
          f"  {'dE, mEh':>8}  {'<L2>':>8}  {'<S2>':>6}"
          f"  {'|[P,Lz]|':>9} {'|[P,Lperp]|':>11}")
    for nk in range(1, len(groups) + 1):
        n_max = sum(len(groups[i]) for i in range(nk))
        c_win, kept = select_no_window(c_no, occ, n_max)
        k = c_win.shape[1]
        if k < max(na, nb):
            print(f"{nk:>7} {k:>4}     0   (fewer columns than "
                  f"electrons; smallest sector opens at k = {max(na, nb)})")
            continue
        ndet = math.comb(k, na) * math.comb(k, nb)
        if ndet > 12000:
            print(f"{nk:>7} {k:>4} {ndet:>6}   (above the exact cap)")
            continue
        wc = window_closure(mol, c_win)
        ops_k = build_operators(mol, mf, c_win, ncore, na, nb,
                                c_core=c_pur)
        ev_k, vec_k = np.linalg.eigh(ops_k["H"].real)
        l2_0 = float(vec_k[:, 0] @ ops_k["L2"].real @ vec_k[:, 0])
        s2_0 = float(vec_k[:, 0] @ ops_k["S2"].real @ vec_k[:, 0])
        print(f"{nk:>7} {k:>4} {ndet:>6}   {ev_k[0] + e_core:>12.6f}"
              f"  {(ev_k[0] - ev[0]) * 1e3:>+8.2f}  {l2_0:>8.4f}"
              f"  {s2_0:>6.3f}"
              f"  {wc['lz']:>9.2e} {wc['lperp']:>11.2e}")


def main():
    wanted = sys.argv[1:] or ["sc", "ti"]
    for name in wanted:
        run_atom(name)


if __name__ == "__main__":
    main()
