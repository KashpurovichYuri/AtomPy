"""The D2h hard-proof probe for the U14e wrapper mol: is PySCF's
symm_orb assignment CORRECT?  Independent of PySCF's own machinery,
label-free, transcription-free:

  1. each D2h operation acts on a central atom's AO functions as a
     SIGN (the real spherical harmonics' nodal planes align with the
     coordinate planes, so every real harmonic is an eigenfunction of
     every diagonal operation) -- the signs are computed NUMERICALLY:
     evaluate each AO on a generic grid, apply the operation to the
     points, read the constant ratio;
  2. the signed operation matrix R_ao (diagonal, nao x nao) is then
     checked against EVERY irrep's symm_orb block: a correct
     assignment means R_ao @ D_k == chi_k(R) * D_k, where chi_k(R)
     is the textbook D2h character table (hardcoded, standard
     ordering);
  3. the character table's own orthogonality
     (sum_R chi_k chi_l = 8 delta_kl) is checked as well.

A pass certifies the D2h assignment end-to-end: the block structure,
the block-to-irrep naming, and the group algebra, all against
numerics + the textbook table, independent of PySCF.

Run from the repo root:
    OMP_NUM_THREADS=8 python scripts/probes/probe_u14e_d2h_verify.py \
        | tee logs/probe_u14e_d2h_verify.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np  # noqa: E402

# the textbook D2h character table, standard ordering:
# operations (E, C2z, C2y, C2x, i, s_xy, s_xz, s_yz);
# the axial system: z = the principal axis
OPS = ("E", "C2z", "C2y", "C2x", "i", "s_xy", "s_xz", "s_yz")
CHAR = {
    "Ag":  (1, 1, 1, 1, 1, 1, 1, 1),
    "B1g": (1, 1, -1, -1, 1, 1, -1, -1),
    "B2g": (1, -1, 1, -1, 1, -1, 1, -1),
    "B3g": (1, -1, -1, 1, 1, -1, -1, 1),
    "Au":  (1, 1, 1, 1, -1, -1, -1, -1),
    "B1u": (1, 1, -1, -1, -1, -1, 1, 1),
    "B2u": (1, -1, 1, -1, -1, 1, -1, 1),
    "B3u": (1, -1, -1, 1, -1, 1, 1, -1),
}
# the operations as the signed coordinate maps (x, y, z) -> ...
MAPS = {
    "E":    (1, 1, 1),
    "C2z":  (-1, -1, 1),
    "C2y":  (-1, 1, -1),
    "C2x":  (1, -1, -1),
    "i":    (-1, -1, -1),
    "s_xy": (1, 1, -1),
    "s_xz": (1, -1, 1),
    "s_yz": (-1, 1, 1),
}


def main() -> int:
    t00 = time.time()
    from pyscf import gto
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    mol = gto.M(atom=f"{spec.symbol} 0 0 0", basis=spec.basis,
                spin=spec.spin, ecp=spec.ecp, symmetry="D2h", verbose=0)
    nao = mol.nao
    print(f"  [mol] the D2h U atom: nao {nao}, the irreps "
          f"{mol.irrep_name}", flush=True)

    # ---- 1. the numeric AO signs per operation -------------------------
    rng = np.random.default_rng(20261124)
    pts = rng.uniform(-0.7, 0.7, size=(9, 3))   # generic, off-plane
    vals = mol.eval_gto("GTOval_sph", pts)      # (npts, nao)
    signs = {}
    for op, (sx, sy, sz) in MAPS.items():
        pts_op = pts * np.array([sx, sy, sz])
        vals_op = mol.eval_gto("GTOval_sph", pts_op)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = vals_op / vals
        per_ao = np.zeros(nao)
        ok = 0
        for i in range(nao):
            r = ratio[:, i]
            r = r[np.abs(r) > 1e-8]
            if len(r) == 0:
                continue
            med = np.median(r)
            if np.max(np.abs(r - med)) < 1e-7 and abs(abs(med) - 1) < 1e-7:
                per_ao[i] = med
                ok += 1
        signs[op] = per_ao
        print(f"  [{op}] the numeric signs determined for {ok}/{nao} "
              f"AOs (the remainder vanish on the grid)", flush=True)

    # ---- 2. R_ao @ D_k == chi_k(R) * D_k for all ops x irreps ----------
    fails = 0
    checks = 0
    for k, name in enumerate(mol.irrep_name):
        d = np.asarray(mol.symm_orb[k])
        for opi, op in enumerate(OPS):
            r_ao = np.diag(signs[op])
            chi = CHAR[name][opi]
            dev = float(np.max(np.abs(r_ao @ d - chi * d)))
            checks += 1
            if dev > 1e-10:
                fails += 1
                print(f"  [FAIL] {name} x {op}: deviation {dev:.2e}",
                      flush=True)
    print(f"  [verdict] {checks} op x irrep checks, {fails} failures "
          f"-- the D2h assignment "
          f"{'CERTIFIED' if fails == 0 else 'REFUTED'}", flush=True)

    # ---- 3. the character table's own orthogonality --------------------
    names = list(mol.irrep_name)
    ortho_ok = all(
        abs(sum(CHAR[a][i] * CHAR[b][i] for i in range(8))
            - (8 if a == b else 0)) == 0
        for a in names for b in names)
    print(f"  [table] the character-table orthogonality: "
          f"{'PASS' if ortho_ok else 'FAIL'}", flush=True)

    # ---- 4. the subduction counts (the group-theory cross-check) -------
    print("  [note] the conflation is group theory: s and dz2/dx2-y2 "
          "are ALL invariant under every D2h operation -- no axis "
          "relabeling separates them; the 7 f functions MUST "
          "distribute over the four u-irreps (2 each + f0 in Au)",
          flush=True)

    good = fails == 0 and ortho_ok
    print(f"\n{'CERTIFIED' if good else 'REFUTED'} "
          f"({time.time() - t00:.0f} s total)", flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
