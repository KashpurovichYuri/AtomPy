#!/usr/bin/env python
"""Basin map: does the CASSCF end answer depend on the start?

Runs free CASSCF from N perturbed starts of the SAME window, then
classifies every endpoint: E, <L^2>, commutator, and the pairwise span
overlap between converged windows. Within a cluster (same energy) the
differences are gauge; across clusters they are the basins of the
Motivation section.

Usage:
    python scripts/basin_map.py --atom C --starts 5
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scipy.linalg import expm

from atomterm.angular import l_matrices  # noqa: F401  (libcint init)
from atomterm.casscf import constrained_casscf
from atomterm.drivers import occupied_core_idx
from atomterm.operators import build_operators
from atomterm.reference import run_rohf
from atomterm.system import AtomConfig, build_mol, ncore_of
from atomterm.windows import build_full_mo, select_valence_ao, span_overlap

CFG = {
    "C": (AtomConfig("C", "ccpvdz", spin=2), ["2p"], 2,
          "elements/C/rungs/rohf.chk"),
    "Ti": (AtomConfig("Ti", "ccpvdz", spin=2), ["3d"], 2,
           "elements/Ti/rungs/rohf.chk"),
}


def perturb(mo, rng, scale):
    """Small random unitary mixing of the MO columns (orthonormality kept)."""
    n = mo.shape[1]
    a = rng.standard_normal((n, n))
    kappa = scale * (a - a.T) / np.linalg.norm(a - a.T)
    return mo @ expm(kappa)


def endpoint_stats(mol, mf, c_win, ncore, na, nb, tag):
    ops = build_operators(mol, mf, c_win, ncore, na, nb)
    res = analyze_tag(ops, tag)
    ground = res.states[0]
    return res, ground


def analyze_tag(ops, tag):
    from atomterm.analysis import analyze
    return analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag,
                   sz=ops["Sz"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atom", default="C", choices=sorted(CFG))
    ap.add_argument("--starts", type=int, default=5)
    ap.add_argument("--scale", type=float, default=0.3)
    args = ap.parse_args()

    cfg, patterns, nelec, rohf_chk = CFG[args.atom]
    rng = np.random.default_rng(2026)

    mol = build_mol(cfg)
    mf = run_rohf(cfg, rohf_chk, ref=None)[1]
    ncore = ncore_of(mol.nelectron, nelec)
    act_idx = select_valence_ao(mol, patterns)
    ncas = len(act_idx)
    na = (nelec + cfg.spin) // 2
    nb = nelec - na
    mo_base = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
    s = mf.get_ovlp()

    rows = []
    wins = []
    for k in range(args.starts):
        start = mo_base if k == 0 else perturb(mo_base, rng, args.scale)
        tag = f"basin_start{k}"
        res = constrained_casscf(
            mf, ncas, nelec, f"elements/{args.atom}/basin_start{k}.chk",
            mo_full=start, name=tag)
        c_win = res.mo_coeff[:, ncore:ncore + ncas]
        wins.append(c_win)
        ops = build_operators(mol, mf, c_win, ncore, na, nb)
        an = analyze_tag(ops, tag)
        rows.append((k, res.energies[0], an.states[0].l2, an.comm_norm_h_lz))
        print(f"start {k}: E = {res.energies[0]:.8f} Eh   <L^2> = "
              f"{an.states[0].l2:.4f}   |[H,Lz]| = {an.comm_norm_h_lz:.2e}",
              flush=True)

    es = [r[1] for r in rows]
    print(f"\nenergy spread = {(max(es) - min(es)) * 1e3:.3f} mEh")

    n = len(wins)
    ov = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            ov[i, j] = ov[j, i] = span_overlap(wins[i], wins[j], s)
    print("\nspan overlap between endpoints:")
    print("      " + "".join(f"   start{j:d}" for j in range(n)))
    for i in range(n):
        print(f"start{i} " + "".join(f"  {ov[i, j]:8.5f}" for j in range(n)))


if __name__ == "__main__":
    main()
