"""Sc/Ti methodology B at matched det counts (ticket 24, user-approved).

Window B = pur core + the valence occupied MOs (4s^2 closed pair + 3d^1
or 3d^2 open) + k virtuals ranked by S-weight on the valence shell rows
({3d,4s,4d,4p}), sized so the sector determinant count MATCHES the
contracted A-lattice points:

    Sc (2,1): w=6 -> C(6,2)*6   =   90  (A: {3d,4s}     90)
              w=11 -> C(11,2)*11 =  605  (A: {3d,4s,4d}  605)
              w=14 -> C(14,2)*14 = 1274  (A: {3d,4d,4s,4p} 1274)
    Ti (3,1): w=6 -> C(6,3)*6   =   120  (A: {3d,4s}     120)
              w=11 -> C(11,3)*11 = 1815  (A: {3d,4s,4d} 1815)
              w=14 -> C(14,3)*14 = 5096  (A: {3d,4d,4s,4p} 5096)

Question answered: does exactness cost energy at matched cost, per atom
class? All sectors <= EXACT_MAX: both sides exact.

Usage: python scripts/run_tm_methb.py [sc|ti]
Threads: OPENBLAS_NUM_THREADS=128 OMP_NUM_THREADS=128 PYSCF_NUM_THREADS=128.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_tm_lattice import RES, anchors
from atomterm import drivers
from atomterm.detbasis import commutator, frobenius_norm
from atomterm.operators import build_operators
from atomterm.system import AtomConfig

def _ti_ref():
    from atomterm.drivers import ReferenceSpec
    return ReferenceSpec(
        occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
        occ_open=("3d",), occ_e_ref=-848.40648875)


CASES = {
    "sc": dict(symbol="Sc", spin=1, sector=(2, 1), ncore=9, ref=None,
               w_list=[6, 11, 14], a_tags=["ds", "dsd", "dsdp"]),
    "ti": dict(symbol="Ti", spin=2, sector=(3, 1), ncore=9,
               ref_factory=_ti_ref, w_list=[6, 11, 14],
               a_tags=["ds", "dsd", "dsdp"]),
}


def main():
    wanted = sys.argv[1:] or list(CASES)
    for name in wanted:
        kw = CASES[name]
        sym, spin = kw["symbol"], kw["spin"]
        na, nb = kw["sector"]
        ncore = kw["ncore"]
        print(f"\n{'#' * 70}\n# {sym} methodology B (matched cost)\n{'#' * 70}")
        cfg = AtomConfig(sym, "ccpvdz", spin=spin)
        chk = str(RES / f"{sym}_ccpvdz_s{spin}_rohf.chk")
        mol, mf = drivers.run_rohf(
            cfg, chk, ref=(kw["ref_factory"]() if kw.get("ref_factory") else None))
        c_pur, contam = drivers.purified_core_orbitals(mol, mf, ncore)
        print(f"pur core: max contamination {max(contam):.2e}")

        closed = drivers.occupied_core_idx(mf, ncore)
        occ2 = list(np.flatnonzero(mf.mo_occ == 2.0))
        val_occ = occ2[ncore:] + list(np.flatnonzero(mf.mo_occ == 1.0))
        assert len(val_occ) == 1 + int(cfg.spin), val_occ
        rows_v = drivers.select_valence_ao(mol, ["3d", "4s", "4d", "4p"])
        s = mol.intor_symmetric("int1e_ovlp")
        virt = [k for k in range(mf.mo_coeff.shape[1]) if mf.mo_occ[k] == 0.0]
        char = {k: float(mf.mo_coeff[:, k][rows_v]
                         @ s[np.ix_(rows_v, rows_v)]
                         @ mf.mo_coeff[:, k][rows_v]) for k in virt}
        rank = sorted(virt, key=lambda k: -char[k])
        print("virtual admission (mo:weight):",
              "  ".join(f"{r}:{char[r]:.3f}" for r in rank[:12]))

        anc = anchors(sym, spin)
        print(f"{'w':>3} {'norb':>4} {'ndet':>6} {'E_ground':>16} "
              f"{'<L^2>':>8} {'<S^2>':>7} {'|[H,Lz]|':>9} {'dEvsPy':>8}  "
              f"vs A-anchor")
        for w, a_tag in zip(kw["w_list"], kw["a_tags"]):
            k = w - len(val_occ)
            cols = val_occ + rank[:k]
            mo = np.hstack([c_pur, mf.mo_coeff[:, cols]])
            ncas = len(cols)
            ndet = math.comb(ncas, na) * math.comb(ncas, nb)
            cas = drivers.run_casci(mf, ncas, (na, nb), mo, chk, nroots=1,
                                    name=f"methb_w{w}",
                                    fci_opts={"pspace_size": ndet})
            ops = build_operators(mol, mf, mo[:, ncore:ncore + ncas],
                                  ncore, na, nb, c_core=mo[:, :ncore])
            ev, vec = np.linalg.eigh(ops["H"].real)
            e_g = ev[0] + cas.core_energy
            l2 = float(vec[:, 0] @ ops["L2"] @ vec[:, 0])
            s2 = float(vec[:, 0] @ ops["S2"] @ vec[:, 0])
            clz = frobenius_norm(commutator(ops["H"], ops["Lz"]))
            d_e = abs(e_g - cas.energies[0])
            a_e = anc.get(a_tag)
            extra = (f"  A({a_tag}) = {a_e:.6f}: "
                     f"{(e_g - a_e) * 1e3:+9.2f} mEh" if a_e else "")
            print(f"{w:>3} {ncas:>4} {ndet:>6} {e_g:>16.8f} {l2:>8.4f} "
                  f"{s2:>7.4f} {clz:>9.2e} {d_e:>8.1e}{extra}")


if __name__ == "__main__":
    main()
