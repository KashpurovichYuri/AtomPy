"""Sc/Ti window lattices (ticket 24): methodology A on the transition
metals, both basis kinds, exact symmetry everywhere.

The experiment: grow the valence radial span at m-complete closure and
watch (a) where the energy goes relative to the CASSCF escape basin
(fact 15d) and the radial ceiling (fact 16b: +24.4/+25.1 mEh on Sc) --
the user's full-space intuition says the exact-symmetry lattice MUST
close that gap once the window includes the escape directions
m-completely; (b) the contracted-vs-decontracted granularity question
(shell = one converged radial vs primitive chains).

Two kinds:
  contracted   shells as radials (select_valence_ao pattern lists); the
               campaign ROHF chks (trusted basins, incl. Ti's gated
               fact-21a reference) are reused.
  decontracted pyscf decontract_basis; fresh ROHF with whole-l-span
               forced occupations (the label-based Ti spec does not
               transfer); windows = primitive bands around the MEASURED
               open-MO weights (the U protocol: never labels).

Every sector is tiny (<= ~1.2k dets contracted, <= ~1k decontracted):
exact eigh + PySCF pspace-exact cross-check, seconds-to-minutes each.
Anchors: the campaign ladder json (casci_avg_pur = the same-substrate
CASCI the lattice's {3d,4s}/{3d,4s} point must reproduce; casscf_pyscf =
the escaped single-root basin; casscf_radial = the radial ceiling).

Usage: python scripts/run_tm_lattice.py [sc|ti ...] [--kind contracted|decontracted]
Threads: OPENBLAS_NUM_THREADS=128 OMP_NUM_THREADS=128 PYSCF_NUM_THREADS=128.
"""

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_u import prim_weights
from atomterm import drivers
from atomterm.analysis import term_symbol
from atomterm.detbasis import commutator, frobenius_norm
from atomterm.drivers import ReferenceSpec
from atomterm.operators import build_operators, window_closure
from atomterm.system import AtomConfig

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
EXACT_MAX = 12000

SHELL_TIERS = [
    ("d", ["3d"]),
    ("ds", ["3d", "4s"]),
    ("dd", ["3d", "4d"]),
    ("dsd", ["3d", "4s", "4d"]),
    ("dsp", ["3d", "4s", "4p"]),
    ("dsdp", ["3d", "4d", "4s", "4p"]),
    # the crossing candidates (ticket 24): one more radial on the axes
    ("dsdpd5", ["3d", "4d", "5d", "4s", "4p"]),
    ("dsdpf", ["3d", "4d", "4s", "4p", "4f"]),
]

# whole-l-span occupations for the decontracted references: [Ar] (s+p)
# closed incl. the 4s pair, all d open (Sc 3d1, Ti 3d2 by electron count)
REF_UNC = ReferenceSpec(occ_closed=("l0", "l1"), occ_open=("l2",))

CASES = {
    "sc": dict(symbol="Sc", spin=1, nelec=3, sector=(2, 1), ncore=9),
    "ti": dict(symbol="Ti", spin=2, nelec=4, sector=(3, 1), ncore=9,
               ref_contract=ReferenceSpec(
                   occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
                   occ_open=("3d",), occ_e_ref=-848.40648875)),
}


def smearing(mol, mf, win_rows, c_core):
    """fact-21a basin detector: max S-overlap of the window AO span with
    the core columns BEFORE projection (good basin ~1e-8, metastable
    ~1.6e-2 on Ti contracted)."""
    s = mol.intor_symmetric("int1e_ovlp")
    worst = 0.0
    for r in win_rows:
        v = np.zeros(mol.nao_nr())
        v[r] = 1.0
        worst = max(worst, float(np.max(np.abs(v @ s @ c_core))))
    return worst


def anchors(atom, spin):
    """Campaign ladder rung energies (Eh) + names, for the delta print."""
    p = RES / f"{atom}_ccpvdz_s{spin}_ladder.json"
    if not p.exists():
        return {}
    d = json.load(open(p))
    return {k: v["e_ground"] for k, v in d.items()
            if isinstance(v, dict) and v.get("e_ground") is not None}


def run_tier(mol, mf, tag, act_rows, ncore, na, nb, chk):
    closed = drivers.occupied_core_idx(mf, ncore)
    c_pur, contam = drivers.purified_core_orbitals(mol, mf, ncore)
    mo = drivers.build_full_mo(mf, closed, act_rows, c_core=c_pur)
    ncas = len(act_rows)
    ndet = math.comb(ncas, na) * math.comb(ncas, nb)
    sm = smearing(mol, mf, act_rows, mo[:, :ncore])
    fci_opts = ({"pspace_size": ndet} if ndet <= EXACT_MAX
                else {"conv_tol": 1e-12})
    cas = drivers.run_casci(mf, ncas, (na, nb), mo, chk, nroots=1,
                            name=tag, fci_opts=fci_opts)
    # build_full_mo returns [core | window | complement] (fact 7): the
    # analysis space is the WINDOW slice only
    c_act = mo[:, ncore:ncore + ncas]
    ops = build_operators(mol, mf, c_act, ncore, na, nb,
                          c_core=mo[:, :ncore])
    ev, vec = np.linalg.eigh(ops["H"].real)
    e_g = ev[0] + cas.core_energy
    l20 = float(vec[:, 0] @ ops["L2"] @ vec[:, 0])
    s20 = float(vec[:, 0] @ ops["S2"] @ vec[:, 0])
    h_norm = max(frobenius_norm(ops["H"]), 1e-300)
    clz = frobenius_norm(commutator(ops["H"], ops["Lz"]))
    cl2 = frobenius_norm(commutator(ops["H"], ops["L2"]))
    wc = window_closure(mol, np.eye(mol.nao_nr())[:, act_rows])
    d_pyscf = abs(e_g - cas.energies[0])
    return dict(tag=tag, ncas=ncas, ndet=ndet, e=e_g, l2=l20, s2=s20,
                clz=clz, cl2=cl2, hn=h_norm, wc=wc, d=d_pyscf, sm=sm,
                contam=max(contam), term=term_symbol(s20, l20))


def unchained(mol):
    """(l -> list of (exp, rows)) for the DECONTRACTED mol: one entry per
    contracted radial, contraction-outer / m-inner rows (fact 16e).
    Official accessors only (the _bas column layout is not a contract)."""
    out = {}
    loc = mol.ao_loc_nr()
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        nctr = mol.bas_nctr(i)
        exps = mol.bas_exp(i)
        for c in range(nctr):
            rows = [loc[i] + c * (2 * l + 1) + m for m in range(2 * l + 1)]
            out.setdefault(l, []).append((float(exps[0]), rows))
    return out


def main():
    args = sys.argv[1:]
    kind = ("decontracted" if "--kind" in args
            and args[args.index("--kind") + 1] == "decontracted"
            else "contracted")
    drop = {i for i, a in enumerate(args) if a.startswith("--")}
    drop |= {i + 1 for i, a in enumerate(args) if a == "--kind"}
    tiers = None
    if "--tiers" in args:
        i = args.index("--tiers")
        tiers = set(args[i + 1].split(","))
        drop |= {i, i + 1}
    wanted = [a for i, a in enumerate(args) if i not in drop] or list(CASES)
    for name in wanted:
        cfg_kw = CASES[name]
        sym, spin = cfg_kw["symbol"], cfg_kw["spin"]
        na, nb = cfg_kw["sector"]
        ncore = cfg_kw["ncore"]
        print(f"\n{'#' * 72}\n# {sym} ({kind})  nelec={cfg_kw['nelec']}"
              f"  sector M_S={abs(na - nb) / 2}\n{'#' * 72}")

        if kind == "contracted":
            cfg = AtomConfig(sym, "ccpvdz", spin=spin)
            chk = str(RES / f"{sym}_ccpvdz_s{spin}_rohf.chk")
            mol, mf = drivers.run_rohf(cfg, chk,
                                       ref=cfg_kw.get("ref_contract"))
            tier_spec = list(SHELL_TIERS)  # (tag, patterns); lazy resolve
            anc = anchors(sym, spin)
        else:
            from pyscf import gto
            from pyscf.gto import mole
            m0 = gto.M(atom=f"{sym} 0 0 0", basis="ccpvdz", spin=spin,
                       verbose=0)
            # decontract_basis splits every shell into primitives EXCEPT
            # one fused core-s block (measured: 19 exps x 1 ctr survives)
            # -- primitive_chain refuses fused slots by design, so the
            # chains below are read straight from _bas: one entry per
            # contracted radial (nprim=1 entries are pure primitives;
            # the fused block is ONE converged radial, honestly kept).
            mu, _ = mole.decontract_basis(m0)
            cfg = AtomConfig(sym, basis=mu._basis, spin=spin,
                             tag=f"{sym}_ccpvdz_unc")
            chk = str(RES / f"{sym}_ccpvdz_unc_s{spin}_rohf.chk")
            # PLAIN ROHF (no whole-l forcing: closing the ENTIRE
            # uncontracted s span smears the core over 20 primitives --
            # measured E = -737.2 vs -759.74; the fact-15c runaway was a
            # FRAC-OCC effect, and the contracted plain ROHF is clean)
            mol, mf = drivers.run_rohf(cfg, chk)
            e_ref_contracted = {"Sc": -759.7362551, "Ti": -848.4064888}[sym]
            assert mf.e_tot < e_ref_contracted + 1e-4, (
                "decontracted ROHF failed to reproduce the mean field: "
                f"{mf.e_tot} vs contracted {e_ref_contracted}")
            print(f"ROHF: E = {mf.e_tot:.8f}  "
                  f"(closed {int((mf.mo_occ == 2).sum())}, "
                  f"open {int((mf.mo_occ == 1).sum())})")
            s = mol.intor_symmetric("int1e_ovlp")
            open_cols = np.flatnonzero(mf.mo_occ == 1.0)
            # valence columns for band measurement: the open MOs PLUS the
            # valence s pair (the highest closed MO = 4s^2 here) -- the
            # open-MO-only weights would leave the s band empty by
            # construction. p has NO occupied valence carrier: its band
            # comes from the VIRTUAL character (the 4p polarization
            # radial lives in virtuals).
            val_cols = list(open_cols) + \
                [int(np.flatnonzero(mf.mo_occ == 2.0)[-1])]
            vir_cols = np.flatnonzero(mf.mo_occ == 0.0)
            chains = unchained(mol)
            band = {}
            for l, nm, cols in [(0, "s", val_cols), (2, "d", val_cols),
                                (1, "p", vir_cols)]:
                meas = []
                for exp, rows in chains.get(l, []):
                    srr = s[np.ix_(rows, rows)]
                    tot = max(sum(float(mf.mo_coeff[:, k][rows] @ srr
                                       @ mf.mo_coeff[:, k][rows])
                                  for k in cols), 0.0)
                    meas.append((exp, rows, tot))
                meas.sort(key=lambda t: -t[2])
                if meas:
                    print(f"{nm} radials by {'open+4s' if l != 1 else 'virtual'}"
                          " S-weight:",
                          " ".join(f"{e:.3f}->{w:.2f}" for e, _, w in meas))
                # RELATIVE threshold: the absolute one keeps the
                # near-linearly-dependent needle radials (fact 16e; the
                # diffuse p primitives crash build_full_mo's Loewdin at
                # a zero eigenvalue)
                wmax = max((w for _, _, w in meas), default=0.0)
                band[l] = [rows for _, rows, w in meas
                           if w >= 0.05 * wmax]
            bd, bs, bp = band.get(2, []), band.get(0, []), band.get(1, [])
            print(f"bands (w >= 0.05): d x{len(bd)}  s x{len(bs)}  "
                  f"p x{len(bp)}")
            tiers = []
            for k in range(1, len(bd) + 1):
                tiers.append((f"d{k}", [r for rs in bd[:k] for r in rs]))
            tiers.append(("dalls",
                          [r for rs in bd for r in rs]
                          + [r for rs in bs for r in rs]))
            tiers.append(("dallsp",
                          [r for rs in bd for r in rs]
                          + [r for rs in bs for r in rs]
                          + [r for rs in bp for r in rs]))
            anc = anchors(sym, spin)

        if anc:
            print("anchors (campaign ladder, Eh):",
                  "  ".join(f"{k}:{v:.6f}" for k, v in sorted(anc.items())))
        base = anc.get("casci_avg_pur")
        print(f"{'tier':>7} {'ncas':>4} {'ndet':>6} {'E_ground':>16} "
              f"{'<L^2>':>8} {'<S^2>':>7} {'term':>6} {'|[H,Lz]|':>9} "
              f"{'|[H,L2]|':>9} {'smear':>8} {'dEvsPy':>8}")
        for item in tier_spec:
            if kind == "contracted":
                t, pats = item
                if tiers and t not in tiers:
                    continue
                try:
                    rows = drivers.select_valence_ao(mol, pats)
                except ValueError as exc:
                    print(f"{t:>7}: skipped ({exc})")
                    continue
            else:
                t, rows = item
                if tiers and t not in tiers:
                    continue
            r = run_tier(mol, mf, f"lat_{kind}_{t}", rows, ncore, na, nb,
                         chk)
            extra = ""
            if base:
                extra = f"  vs casci_avg_pur: {(r['e'] - base) * 1e3:+9.2f} mEh"
                for key, nm in [("casscf_pyscf", "escape"),
                                ("casscf_radial", "ceil"),
                                ("casscf_relax_pur", "relaxP")]:
                    if key in anc:
                        extra += f"  {nm}:{(r['e'] - anc[key]) * 1e3:+8.2f}"
            print(f"{t:>7} {r['ncas']:>4} {r['ndet']:>6} {r['e']:>16.8f} "
                  f"{r['l2']:>8.4f} {r['s2']:>7.4f} {r['term']:>6} "
                  f"{r['clz']:>9.2e} {r['cl2']:>9.2e} {r['sm']:>8.1e} "
                  f"{r['d']:>8.1e}{extra}")


if __name__ == "__main__":
    main()
