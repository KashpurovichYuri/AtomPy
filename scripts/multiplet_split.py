#!/usr/bin/env python
"""Multiplet split: the SAME window read through the raw and the averaged
folds. Both folds' Hamiltonians are diagonalized in our determinant basis
and labeled by <L_z> (a good quantum number in both folds): the averaged
fold is degenerate per term; the raw fold splits it by |M_L| -- the
Motivation branching as a table.

Usage: python scripts/multiplet_split.py [--nroots 9]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.drivers import occupied_core_idx
from atomterm.operators import build_operators
from atomterm.reference import run_rohf
from atomterm.system import AtomConfig, build_mol, ncore_of
from atomterm.windows import build_full_mo, select_valence_ao


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atom", default="Ti")
    ap.add_argument("--window", default="3d4s")
    ap.add_argument("--nroots", type=int, default=9)
    args = ap.parse_args()

    cfg = AtomConfig(args.atom, "ccpvdz", spin=2)
    mol = build_mol(cfg)
    mf = run_rohf(cfg, f"elements/{args.atom}/rungs/rohf.chk",
                  ref=None)[1]
    patterns = {"3d4s": ["3d", "4s"], "3d": ["3d"]}[args.window]
    nelec = {"3d4s": 4, "3d": 2}[args.window]
    ncore = ncore_of(mol.nelectron, nelec)
    act_idx = select_valence_ao(mol, patterns)
    ncas = len(act_idx)
    na = (nelec + cfg.spin) // 2
    nb = nelec - na

    c_win = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
    c_win = c_win[:, ncore:ncore + ncas]

    h_raw = np.real(build_operators(mol, mf, c_win, ncore, na, nb,
                                    spherical_core=False)["H"])
    h_avg = np.real(build_operators(mol, mf, c_win, ncore, na, nb,
                                    spherical_core=True)["H"])

    from atomterm.detbasis import determinant_basis
    from atomterm.windows import lm_sectors
    dets = determinant_basis(ncas, na, nb)
    sectors = lm_sectors(mol)
    local = {r: k for k, r in enumerate(act_idx)}
    m_of_local = {}
    for (l, member), rows in sectors.items():
        for r in rows:
            if r in local:
                m_of_local[local[r]] = member - l
    ml_of_det = [sum(m_of_local[so % ncas] for so in det) for det in dets]

    def blocks(h):
        out = {}
        for ml in sorted(set(ml_of_det)):
            idx = [k for k, m in enumerate(ml_of_det) if m == ml]
            sub = h[np.ix_(idx, idx)]
            out[ml] = np.linalg.eigvalsh(sub)
        return out

    raw_blocks = blocks(h_raw)
    avg_blocks = blocks(h_avg)

    print("\n== raw fold: M_L blocks ==")
    print("| M_L | N states | E levels, Eh |")
    print("|---|---|---|")
    for ml in sorted(raw_blocks):
        levels = ", ".join(f"{e:.6f}" for e in raw_blocks[ml])
        print(f"| {ml:+d} | {len(raw_blocks[ml])} | {levels} |")

    print("\n== averaged fold: M_L blocks ==")
    print("| M_L | N states | E levels, Eh |")
    print("|---|---|---|")
    for ml in sorted(avg_blocks):
        levels = ", ".join(f"{e:.6f}" for e in avg_blocks[ml])
        print(f"| {ml:+d} | {len(avg_blocks[ml])} | {levels} |")

    for ml in sorted(raw_blocks):
        if len(raw_blocks[ml]) > 1:
            sp = (raw_blocks[ml][-1] - raw_blocks[ml][0]) * 1e3
            print(f"raw fold, M_L = {ml:+d}: split {sp:.2f} mEh")


if __name__ == "__main__":
    main()
