"""U methodology-B scan: raw-ROHF-orbital windows — the broken-symmetry
size-vs-extension curve (ticket 24, light tier).

Window B = the pur core + the 4 open ROHF MOs + k virtuals ranked by
S-weight on the f+d AO rows (character selection — canonical order is
junk on this basis, ticket-21 finding 3). Unlike methodology A (m-complete
primitive lattices) these windows are m-INCOMPLETE: <L^2> goes off-grid
and the commutators / window_closure measure the closure defect directly.
Every step is small (C(4+k,4) dets <= 495), exact eigh, seconds each.

Per admitted virtual k the curve reads:
    E_ground(k)      the broken-symmetry optimum at this cost
    E(k) - E(0)      fold-consistent descent (k=0 IS the bare ROHF
                     determinant: 4 open MOs, 1 determinant)
    E(k) - E_ROHF    the physical crossing statement (carries the
                     pur-fold bookkeeping offset -- measured at k=0)
    <L^2> ground     the breaking meter (A gives exactly 72 for 5L-deg)
    ||[H,Lz]|| etc   the many-body closure defect of the MO window
    ||[P,Lz]|| etc   the one-particle span certificate (fact 19)

The k=8 point (12 MOs, 495 dets) is MATCHED-COST to methodology A's
1f+1d lattice point (12 AOs, 495 dets): "does exactness cost energy at
matched cost" on U. Sectors are stored in the campaign chk under
atomterm/methb_k<k>/ (reproducibility convention).

Usage: python scripts/run_u_methb.py [k ...]   (default: 0 1 2 3 4 6 8)
Threads: OPENBLAS_NUM_THREADS=128 OMP_NUM_THREADS=128 PYSCF_NUM_THREADS=128.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_u import CFG, CHK, NCORE, REF
from atomterm import drivers
from atomterm.detbasis import commutator, frobenius_norm
from atomterm.operators import (
    build_operators,
    parity_operator,
    population_operator,
    window_closure,
)

K_STEPS = [0, 1, 2, 3, 4, 6, 8]
NTRACK = 24  # roots scanned for the 5L-deg-like (<L^2> ~ 72) member


def f_rows(mol):
    """The f-chain AO rows (for <N_f>: the honest all-chain f population)."""
    return [r for _, rs in drivers.primitive_chain(mol, 3) for r in rs]


def valence_rows(mol):
    """The VALENCE-BAND f+d rows (the same bands as the 4f2d lattice
    tier): the character metric for virtual admission. All-chain rows
    would admit the near-core tight f MOs (measured: canonical mo123,
    eps = +38.8 Eh, f-character 1.000 on the tight primitives) whose
    singles sit tens of Eh up and never mix -- the curve would be flat."""
    rows = []
    for l, band in [(3, (0.2, 4.0)), (2, (0.1, 0.5))]:
        rows.extend(drivers.select_primitives(mol, l, exponent_band=band))
    return rows


def main():
    wanted = [int(a) for a in sys.argv[1:]] or K_STEPS
    mol, mf = drivers.run_rohf(CFG, CHK, ref=REF)
    e_rohf = mf.e_tot
    print(f"reference: E = {e_rohf:.8f} Eh  (cached gated Newton)")

    closed = list(np.flatnonzero(mf.mo_occ == 2.0))
    assert len(closed) == NCORE, (len(closed), NCORE)
    c_pur, contam = drivers.purified_core_orbitals(mol, mf, NCORE)
    print(f"pur core: max contamination {max(contam):.2e}")

    open_cols = list(np.flatnonzero(mf.mo_occ == 1.0))
    assert len(open_cols) == 4, open_cols

    # the measurement: S-weight of every MO on the valence-band rows
    rows_v = valence_rows(mol)
    rows_f = f_rows(mol)
    s = mol.intor_symmetric("int1e_ovlp")

    def s_weight(col, rr):
        c = mf.mo_coeff[:, col]
        return float(c[rr] @ (s @ c)[rr])

    virt = [k for k in range(mf.mo_coeff.shape[1]) if mf.mo_occ[k] == 0.0]
    char = {k: s_weight(k, rows_v) for k in virt}
    rank = sorted(virt, key=lambda k: -char[k])
    print("open MOs valence-band character:",
          " ".join(f"{s_weight(c, rows_v):.3f}" for c in open_cols))
    print("virtual admission order (rank:mo:eps:valence S-weight):")
    print("  " + "  ".join(
        f"{i}:mo{r}:{mf.mo_energy[r]:+.3f}:{char[r]:.3f}"
        for i, r in enumerate(rank[:8])))

    results = {}
    for k in wanted:
        cols = open_cols + rank[:k]
        c_win = mf.mo_coeff[:, cols]
        # canonical MOs are exactly orthonormal and S-orthogonal to the
        # core span; the pur core spans the same space -> no re-rotation
        mo = np.hstack([c_pur, c_win])
        ncas = len(cols)
        ndet = math.comb(ncas, 4)

        ops = build_operators(mol, mf, mo[:, NCORE:], NCORE, 4, 0,
                              c_core=mo[:, :NCORE])
        h = ops["H"].real
        ev, vec = np.linalg.eigh(h)

        # PySCF side + core energy (same convention as the lattice run):
        # k=0 is a 1-det CASCI (exact trivially)
        cas = drivers.run_casci(mf, ncas, (4, 0), mo, CHK,
                                nroots=min(NTRACK, ndet),
                                name=f"methb_k{k}")
        core = cas.core_energy
        d_e = float(np.max(np.abs(ev[:1] + core - cas.energies[:1])))

        # the reference determinant inside the CI space: its diagonal
        # (Brillouin-decoupled from singles; the mixing gain reads as
        # ground < det) and the ground CI vector's top components
        from atomterm.detbasis import determinant_basis
        dets = determinant_basis(ncas, 4, 0)
        ref = dets.index(tuple(range(4)))
        e_det = h[ref, ref] + core
        top = np.argsort(-np.abs(vec[:, 0]))[:3]
        gtxt = " ".join(f"{dets[i]}:{vec[i, 0]:+.3f}" for i in top)

        nroots = min(NTRACK, ndet)
        l2s = np.einsum("ik,ij,jk->k", vec[:, :nroots].conj(), ops["L2"],
                        vec[:, :nroots]).real
        near = int(np.argmin(np.abs(l2s - 72.0)))
        nf_op = population_operator(mol, mo[:, NCORE:], rows_f, 4, 0)
        nf0 = float(np.einsum("ik,ij,jk->", vec[:, :1].conj(), nf_op,
                              vec[:, :1]).real)
        par = float(np.diag(parity_operator(mol, mo[:, NCORE:], 4, 0))[0])
        h_norm = max(frobenius_norm(h), 1e-300)
        clz = frobenius_norm(commutator(h, ops["Lz"]))
        cl2 = frobenius_norm(commutator(h, ops["L2"]))
        wc = window_closure(mol, c_win)

        e_g = ev[0] + core
        results[k] = dict(ncas=ncas, ndet=ndet, e_g=e_g, e_det=e_det,
                          l2=l2s[0], near=near, l2near=l2s[near],
                          enear=ev[near] + core, clz=clz, cl2=cl2,
                          hn=h_norm, wc=wc, de=d_e, nf=nf0, par=par,
                          gtxt=gtxt)
        print(f"\n[k={k}] norb={ncas} ndet={ndet}")
        print(f"  E_ground      = {e_g:.6f} Eh   <N_f> = {nf0:.3f}   "
              f"<L^2> = {l2s[0]:.4f}   parity = {par:+.0f}")
        print(f"  E_ref det     = {e_det:.6f}   ground-det = "
              f"{(e_g - e_det) * 1e3:+.3f} mEh mixing gain")
        print(f"  ground vec: {gtxt}")
        print(f"  5L-like root #{near}: <L^2> = {l2s[near]:.4f}   "
              f"E = {ev[near] + core:.6f}")
        print(f"  ||[H,Lz]|| = {clz:.2e} (rel {clz / h_norm:.2e})   "
              f"||[H,L2]|| = {cl2:.2e} (rel {cl2 / h_norm:.2e})")
        print(f"  span: ||[P,Lz]|| = {wc['lz']:.2e}   "
              f"||[P,Lperp]|| = {wc['lperp']:.2e}")
        print(f"  max|dE| vs PySCF (ground) = {d_e:.2e}")

    if 0 in results:
        print(f"\n{'=' * 72}\n"
              "broken-symmetry descent curve (U, sector (4,0), pur fold)\n"
              f"{'=' * 72}")
        e0 = results[0]["e_det"]
        print(f"reference determinant (k=0) vs ROHF total: "
              f"{(e0 - e_rohf) * 1e3:+.2f} mEh (pur-fold bookkeeping)")
        print(f"{'k':>3} {'norb':>4} {'ndet':>5} {'E_ground':>13} "
              f"{'vs ref':>9} {'vs ROHF':>9} {'mix(mEh)':>9} {'<L2>g':>8} "
              f"{'<L2>5L':>8} {'|[H,Lz]|':>9}")
        for k in wanted:
            r = results[k]
            print(f"{k:>3} {r['ncas']:>4} {r['ndet']:>5} {r['e_g']:>13.6f} "
                  f"{(r['e_g'] - e0) * 1e3:>+9.2f} "
                  f"{(r['e_g'] - e_rohf) * 1e3:>+9.2f} "
                  f"{(r['e_g'] - r['e_det']) * 1e3:>+9.3f} "
                  f"{r['l2']:>8.4f} {r['l2near']:>8.4f} "
                  f"{r['clz']:>9.1e}")


if __name__ == "__main__":
    main()
