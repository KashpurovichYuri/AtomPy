#!/usr/bin/env python
"""Multiplet split: the SAME window read through the raw / averaged /
fractional folds, with exact M_L labels.

Real shell AOs are +-|m| combinations, NOT m eigenstates (finding 5),
so determinant-sum labeling over real AOs is invalid. Instead the
per-shell unitary U^(l) = angular.real_to_complex(l) is embedded in
the window-orbital basis and lifted to the Fock space through its
one-body Slater-Condon matrix W (the second-quantized transform).
H_complex = W^dag H W lives in the |l,m> basis where Lz is diagonal:
each fold's roots are labeled by their resolved M_L (degenerate
clusters diagonalized in Lz) and their <Lz^2>.

Ladder validation (--validate): on an isotropic fold the CAS(2,5)
Ms=1 sector is 3F(7) + 3P(3); the resolved <Lz^2> ladder must read
0,1,1,4,4,9,9 | 0,1,1 -- the transform's end-to-end test.

Usage:
  python scripts/multiplet_split.py --atom Ti --window 3d4s
  python scripts/multiplet_split.py --atom C --window 3d --validate
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.angular import real_to_complex
from atomterm.detbasis import determinant_basis
from atomterm.drivers import occupied_core_idx
from atomterm.operators import build_operators
from atomterm.reference import run_frac_occ, run_rohf
from atomterm.system import AtomConfig, build_mol, ncore_of
from atomterm.windows import build_full_mo, select_valence_ao

TOL = 1e-8  # degenerate-cluster resolution threshold (Eh)


def _ao_to_bas(mol) -> dict:
    out = {}
    loc = mol.ao_loc_nr()  # contracted shells may hold nctr x (2l+1) AOs
    for b in range(mol.nbas):
        for ao in range(loc[b], loc[b + 1]):
            out[ao] = b
    return out


def shell_u(mol, act_idx: list[int]) -> np.ndarray:
    """The real->complex unitary embedded in the window-orbital basis.

    The active AOs are grouped by angular momentum l (each shell's AOs
    contiguous, in PySCF real order); each group is transformed by
    U^(l) = angular.real_to_complex(l).
    """
    n = len(act_idx)
    u = np.zeros((n, n), dtype=complex)
    bas_of = _ao_to_bas(mol)
    l_vals = [int(mol._bas[bas_of[ao], 1]) for ao in act_idx]
    i = 0
    while i < n:
        l = l_vals[i]
        j = i
        while j < n and l_vals[j] == l:
            j += 1
        u[i:j, i:j] = real_to_complex(l)
        i = j
    assert np.allclose(u.conj().T @ u, np.eye(n), atol=1e-12), "U not unitary"
    return u


def fock_transform(u_spatial: np.ndarray, dets) -> np.ndarray:
    """The Fock-space basis change W of the one-body unitary u (spatial):
    the WEDGE representation, W_IJ = det(u_so[I_occ, J_occ]) over the two
    determinants' occupied spin-orbitals. (The Slater-Condon one-body
    matrix SUM u a^+a is the operator, not the basis change.)"""
    n = u_spatial.shape[0]
    u_so = np.block([[u_spatial, np.zeros((n, n))],
                     [np.zeros((n, n)), u_spatial]]).astype(complex)
    ndet = len(dets)
    w = np.zeros((ndet, ndet), dtype=complex)
    for a, I in enumerate(dets):
        rows = np.asarray(I)
        for b, J in enumerate(dets):
            w[a, b] = np.linalg.det(u_so[np.ix_(rows, np.asarray(J))])
    return w


def resolve(h: np.ndarray, lz: np.ndarray, w: np.ndarray, dets,
            tol: float = TOL):
    """Diagonalize h in the complex (|l,m>) basis; inside each degenerate
    cluster resolve M_L by diagonalizing Lz (H and Lz commute on every
    m-diagonal fold), then read <Lz^2> off the Lz eigenvector. Returns
    per-root (energy, M_L, <Lz^2>)."""
    lzc = w.conj().T @ lz @ w
    lzc2 = lzc @ lzc
    evals, evecs = np.linalg.eigh(w.conj().T @ h @ w)
    roots = []
    start = 0
    for i in range(1, len(evals) + 1):
        if i == len(evals) or evals[i] - evals[start] > tol:
            cluster = evecs[:, start:i]
            m_vals, u_vecs = np.linalg.eigh(
                cluster.conj().T @ lzc @ cluster)
            for k in range(i - start):
                v = cluster @ u_vecs[:, k]
                lz2 = float((v.conj() @ lzc2 @ v).real)
                roots.append((float(evals[start + k]), float(m_vals[k]),
                              lz2))
            start = i
    return roots


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atom", default="Ti")
    ap.add_argument("--window", default="3d4s")
    ap.add_argument("--nroots", type=int, default=9)
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()

    cfg = AtomConfig(args.atom, "ccpvdz", spin=2)
    mol = build_mol(cfg)
    mf = run_rohf(cfg, f"elements/{args.atom}/rungs/rohf.chk", ref=None)[1]
    windows = {"2p": (["2p"], 2), "3p": (["3p"], 2),
               "3d": (["3d"], 2), "3d4s": (["3d", "4s"], 4)}
    patterns, nelec = windows[args.window]
    ncore = ncore_of(mol.nelectron, nelec)
    act_idx = select_valence_ao(mol, patterns)
    ncas = len(act_idx)
    na = (nelec + cfg.spin) // 2
    nb = nelec - na

    c_win = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
    c_win = c_win[:, ncore:ncore + ncas]

    ops_avg = build_operators(mol, mf, c_win, ncore, na, nb,
                              spherical_core=True)
    ops_raw = build_operators(mol, mf, c_win, ncore, na, nb,
                              spherical_core=False)
    dets = determinant_basis(ncas, na, nb)
    lz = ops_avg["Lz"]

    u = shell_u(mol, act_idx)
    w = fock_transform(u, dets)

    folds = {"raw": np.real(ops_raw["H"]),
             "averaged": np.real(ops_avg["H"])}
    if args.window == "3d4s":
        _, mf_frac = run_frac_occ(
            cfg, f"elements/{args.atom}/rungs/fracocc.chk", nelec, ncas,
            open_idx=act_idx)
        folds["fractional"] = np.real(build_operators(
            mol, mf_frac, c_win, ncore, na, nb, spherical_core=True)["H"])

    if args.validate:
        # p^2 Ms=1 is 3P only (L=1: 0,1,1); d^2 Ms=1 is 3F(7)+3P(3):
        # L=3 -> 0,1,1,4,4,9,9 plus L=1 -> 0,1,1
        want_map = {("C", "2p"): [0, 1, 1],
                    ("Ti", "3d"): [0, 1, 1, 4, 4, 9, 9, 0, 1, 1]}
        want = sorted(want_map[(args.atom, args.window)])
        roots = resolve(folds["averaged"], lz, w, dets)
        got = sorted(round(r[2], 6) for r in roots)
        ok = np.allclose(got, want, atol=1e-6)
        print(f"# VALIDATE {args.atom} CAS({nelec},{ncas}) averaged: "
              f"ladder {'OK' if ok else 'MISMATCH'}")
        print("# got     :", got)
        print("# expected:", want)
        sys.exit(0 if ok else 1)

    # consistency print: the round trip must restore the real H
    hc = w.conj().T @ folds["averaged"] @ w
    rt = np.abs(w @ hc @ w.conj().T - folds["averaged"]).max()
    print(f"# round-trip max|W Hc W^dag - H| = {rt:.2e}")

    for name, h in folds.items():
        roots = resolve(h, lz, w, dets)[:args.nroots]
        print(f"\n== {name} fold: per-root (resolved) ==")
        print("| # | E, Eh | M_L | <Lz^2> |")
        print("|---|---|---|---|")
        for k, (e, ml, lz2) in enumerate(roots):
            ml_s = f"{ml:+.1f}" if abs(ml - round(ml)) < 1e-6 \
                else f"{ml:+.3f}"
            print(f"| {k} | {e:.6f} | {ml_s} | {lz2:.4f} |")


if __name__ == "__main__":
    main()
