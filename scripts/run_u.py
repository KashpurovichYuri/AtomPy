"""U campaign (u_14e): Stuttgart RLC ECP78 + cut-Dyall basis.

Serial stages, each inspectable: gated Newton reference (cached) ->
window selection BY MEASUREMENT (the 5f lives on mid-tight primitives,
kin 3.0/6.5; the diffuse end is empty) -> CASCI term tables with parity
labels in the Ms=2/1/0 sectors -> closure + commutator fingerprints.

Windows (f primitives by open-MO S-weight, d by the same):
  tier 1a: {f #5, d #3}              12 AOs   (structure sanity)
  tier 1b: {f #4, #5, d #3}          19 AOs   (the real tier-1 table)
7s^2 is FROZEN in the core (ncore = 5: 6s 6p x3 7s): 7s^2 is 1S-closed,
so the ground term of 5f^3 6d^1 7s^2 IS the term of 5f^3 6d^1, and the
f^4 manifold lives in the same window.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
from pyscf import gto

from atomterm.system import AtomConfig, ncore_of
from atomterm import drivers
from atomterm.analysis import analyze, format_report, term_symbol
from atomterm.drivers import ReferenceSpec
from atomterm.operators import build_operators, parity_operator, window_closure

HERE = Path(__file__).resolve().parent
ECP = (HERE / ".." / "data" / "u_stuttgart_rlc_ecp78.txt").read_text()
BAS = {"U": gto.parse((HERE / ".." / "data" / "u_dyall_cut_basis.txt")
                      .read_text())}

CFG = AtomConfig(symbol="U", basis=BAS, ecp=ECP, spin=4,
                 tag="U_stuttgart78_dyall")
REF = ReferenceSpec(occ_closed=("l0", "l1"), occ_open=("l2", "l3"),
                    occ_e_ref=-50.38185233, scf_method="newton")
OUT = Path("results/u_14e")
OUT.mkdir(parents=True, exist_ok=True)
NELEC_ACT, NCORE = 4, 5  # 5f^3 6d^1 active; 6s 6p^3 7s core
CLOSED_IDX = None  # set in main(): the OCCUPIED closed block (canonical
# order is meaningless on near-degenerate uncontracted bases)


def prim_rows(mol, l, kin_lo, kin_hi):
    """AO rows of the primitives of angular momentum l with kinetic
    diagonal in [kin_lo, kin_hi] (m-complete by construction)."""
    kin = mol.intor("int1e_kin").diagonal()
    rows = [r for (ll, m), rws in drivers.lm_sectors(mol).items()
            if ll == l for r in rws]
    return sorted(r for r in rows if kin_lo <= kin[r] <= kin_hi)


def run_sector(mol, mf, win, na, nb, nroots, tag):
    ncas = len(win)
    ndet = math.comb(ncas, na) * math.comb(ncas, nb)
    mo = drivers.build_full_mo(mf, CLOSED_IDX, win)
    # analyze the ORTHONORMALIZED window (build_full_mo Loewdin-rotates
    # overlapping same-l primitives): span- and closure-preserving, and
    # the basis the CASCI actually solved in. RAW overlapping columns
    # break every product operator (measured: L^2 eigenvalues 172.98,
    # non-quantized -- fact 4's metric trap, first multi-radial window)
    c_act = mo[:, NCORE:NCORE + ncas]
    cas = drivers.run_casci(mf, ncas, (na, nb), mo,
                            str(OUT / "U_stuttgart78_dyall_rohf.chk"),
                            nroots=min(nroots, ndet), name=tag)
    ops = build_operators(mol, mf, c_act, NCORE, na, nb,
                          c_core=mo[:, :NCORE])  # the OCCUPIED core block
    # (build_full_mo's columns) -- the DEFAULT c_core is mf.mo_coeff[:,
    # :ncore] = canonical order, which is junk on this near-degenerate
    # uncontracted basis (the same trap occupied_core_idx guards)
    res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag)
    par = np.diag(parity_operator(mol, c_act, na, nb))
    ev, vec = np.linalg.eigh(ops["H"])
    p_exp = np.einsum("ik,i,jk->k", vec.conj(), par, vec).real
    ours = np.sort(ev)[: len(cas.energies)] + cas.core_energy
    max_de = float(np.max(np.abs(ours - np.sort(cas.energies))))
    title = f"[{tag}] ndet={ndet}  max|dE| vs PySCF = {max_de:.2e}"
    text = format_report(res, title)
    text += "\n  #      parity\n"
    for k in range(min(30, ndet)):
        text += (f"  {k:2d}   {'odd' if p_exp[k] < 0 else 'even'}   "
                 f"{term_symbol(res.states[k].s2, res.states[k].l2)}"
                 f"{'°' if p_exp[k] < 0 else ''}\n")
    print(text)
    return res


def main():
    mol, mf = drivers.run_rohf(CFG, str(OUT / "U_stuttgart78_dyall_rohf.chk"),
                               ref=REF)
    print(f"reference: E = {mf.e_tot:.8f} Eh  (5f^3 6d^1 7s^2, Newton-gated)")
    closed_idx = list(np.flatnonzero(mf.mo_occ == 2.0))
    assert len(closed_idx) == NCORE, (len(closed_idx), NCORE)
    # canonical order is MEANINGLESS on this basis (near-degenerate junk
    # virtuals first): the core is the OCCUPIED closed block, or the fold
    # silently uses five garbage columns (measured: max|dE| = 14.5 Eh)
    global CLOSED_IDX
    CLOSED_IDX = closed_idx

    tiers = [
        ("1a", prim_rows(mol, 3, 5.0, 8.0) + prim_rows(mol, 2, 0.5, 0.7)),
        ("1b", prim_rows(mol, 3, 2.5, 8.0) + prim_rows(mol, 2, 0.5, 0.7)),
    ]
    for tier, win in tiers:
        cl = window_closure(mol, np.eye(mol.nao_nr())[:, win])
        print(f"\n{'=' * 68}\ntier {tier}: {len(win)} AOs  "
              f"(f: {sum(1 for r in win if r in prim_rows(mol, 3, 0, 1e9))}, "
              f"rest d)  raw-window closure lz={cl['lz']:.1e} "
              f"lperp={cl['lperp']:.1e}\n{'=' * 68}")
        for na, nb, nroots in [(4, 0, 30), (3, 1, 40), (2, 2, 40)]:
            run_sector(mol, mf, win, na, nb, nroots,
                       f"tier{tier}_ms{na - nb}")


if __name__ == "__main__":
    main()
