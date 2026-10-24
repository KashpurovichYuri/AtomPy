"""Ticket 28: the completion runs, executed AS DECLARED STUDIES through
the study runner (atomterm.study) -- declared product, per-cell budgets,
quarantine, exact-eigh labels (diagnostics="full"), results rendered
into the study markers of elements/<Atom>/TABLES.md.

The studies (ticket 28, approved stop rules):

    sc / ti  NO-group windows (from the term-averaged density of the
             {3d,4d,4s,4p} window CASCI ground term) x
             {casci, relax, sa-casscf} on the pur substrate; read
             against the shell-lattice anchors. All cells are cheap
             (<= 5096 dets).
    u24      U NO windows g1..g2 (12 col) and g1..g4 (24 col = the full
             window) x {sa-casscf, relax}, sector (4,0). The SA-24 cell
             is the long one (17 roots x 10,626 dets per macro) --
             budget declared: max_cycle_macro = 60. Saturation check of
             the SA curve (SA-12 = -50.423810 exact 5L) and the
             relax-curve 24-col point (relax-12 = -49.899729,
             relax-17 = -49.900366).
    six13    U CAS(6,13), sector (5,1): the 7s pair promoted out of the
             pur core into a {d5, f7, 7s} NO window -- does the 7s^2
             pair, made active, move the 5L verdict (CASCI said inert,
             -0.58 mEh)? SA-CASSCF, 17 roots expected, 16,731 dets.

Stop rules (approved 2026-10-24): (a) certificate stop -- a bigger
window whose optimizer gauge breaks m-completeness is reported, not
pushed; (b) saturation stop -- window-independent term identity at mEh;
(c) cost stop -- the declared budgets above.

Usage: python scripts/run_completion.py <stage> [...]   # sc ti u24 six13
Threads: OPENBLAS_NUM_THREADS=64 OMP_NUM_THREADS=64 PYSCF_NUM_THREADS=64.
"""

import math
import re
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_tm_methb import CASES, RES as TM_RES
from run_u import CFG, CHK, NCORE as U_NCORE, REF
from run_u_sa import (NA, NB, build_nos, setup as u_setup)
from atomterm import drivers
from atomterm.study import StudySpec, plan_study, render, run_study
from atomterm.system import AtomConfig

REPO = Path(__file__).resolve().parent.parent
U_RES = Path("results/u_14e")
TI_REF = dict(occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
              occ_open=("3d",), occ_e_ref=-848.40648875)


# --------------------------------------------------------------- plumbing
def inject(atom: str, block: str) -> None:
    """Replace the content between the study markers of the element file."""
    p = REPO / f"elements/{atom}/TABLES.md"
    src = p.read_text()
    pat = re.compile(r"(<!-- study:begin -->[^\n]*\n)(.*?)(<!-- study:end -->)",
                     re.S)
    assert pat.search(src), f"no study markers in {p}"
    p.write_text(pat.sub(lambda m: m.group(1) + block + m.group(3), src),
                 encoding="utf-8")
    print(f"[inject] {p} updated")


def execute(spec: StudySpec, mol, mf, out: Path, c_core, atom: str,
            extra: str = "") -> dict:
    print(plan_study(spec, mol))
    cells = run_study(spec, mol, mf, out, c_core=c_core)
    ok = {k: c for k, c in cells.items() if c["status"] == "ok"}
    for k, c in ok.items():
        print(f"  {k}: E = {c['energy']:.6f}  <L2> = {c.get('l2')}"
              f"  span = ({c['span_lz']:.1e}, {c['span_lperp']:.1e})")
    for k, c in cells.items():
        if c["status"] != "ok":
            print(f"  {k}: ERROR {c['error'][:90]}")
    inject(atom, extra + render(cells))
    return cells


# ------------------------------------------------------------- Sc / Ti
def stage_tm(name: str) -> None:
    from atomterm.natural import (ci_rdm1, natural_orbitals,
                                  occupation_groups, select_no_window)
    kw = CASES[name]
    sym, spin = kw["symbol"], kw["spin"]
    na, nb = kw["sector"]
    ncore = kw["ncore"]
    print(f"\n{'#' * 66}\n# {sym}: NO-group study (ticket 28)\n{'#' * 66}")
    cfg = AtomConfig(sym, "ccpvdz", spin=spin)
    chk = str(TM_RES / f"{sym}_ccpvdz_s{spin}_rohf.chk")
    mol, mf = drivers.run_rohf(
        cfg, chk, ref=(kw["ref_factory"]() if kw.get("ref_factory") else None))
    c_pur, _ = drivers.purified_core_orbitals(mol, mf, ncore)
    closed = list(np.flatnonzero(mf.mo_occ == 2.0))[:ncore]

    # the source window + its term-averaged density (the run_tm_no recipe)
    patterns = ["3d", "4d", "4s", "4p"]
    rows = drivers.select_valence_ao(mol, patterns)
    mo = drivers.build_full_mo(mf, closed, rows, c_core=c_pur)
    c_act = mo[:, ncore:ncore + len(rows)]
    ops = drivers.build_operators(mol, mf, c_act, ncore, na, nb,
                                  c_core=c_pur)
    ev, vec = np.linalg.eigh(ops["H"].real)
    dets = ops["dets"]
    nlow = min(48, len(dets))
    l2r = np.array([vec[:, r] @ ops["L2"].real @ vec[:, r]
                    for r in range(nlow)])
    block = [r for r in range(nlow)
             if abs(l2r[r] - l2r[0]) < 1e-6 and ev[r] - ev[0] < 1e-6]
    print(f"source window ground: {len(block)} roots, <L2> = {l2r[0]:.4f}")
    w = np.zeros(nlow)
    w[block] = 1.0 / len(block)
    gamma = ci_rdm1(dets, vec[:, :nlow], len(rows), weights=w)
    c_no, occ = natural_orbitals(c_act, gamma)
    groups = occupation_groups(occ)
    print(f"NO groups: {[len(g) for g in groups]}")

    windows = []
    for nk in range(1, len(groups) + 1):
        n_max = sum(len(groups[i]) for i in range(nk))
        c_win, _ = select_no_window(c_no, occ, n_max)
        ndet = math.comb(c_win.shape[1], na) * math.comb(c_win.shape[1], nb)
        assert ndet <= 12000, (nk, ndet)   # cost stop (approved cap)
        windows.append(_no_win(f"no{nk}", c_win))
    spec = StudySpec(label=f"{sym.lower()}_no", ncore=ncore, na=na,
                     nb=nb,
                     windows=tuple(windows),
                     treatments=("casci", "relax", "sa-casscf"),
                     budgets={"sa-casscf": {"max_cycle_macro": 60}},
                     diagnostics="full")
    execute(spec, mol, mf, TM_RES, c_pur, sym)


def _no_win(label: str, c_win) -> object:
    from atomterm.windows import Window
    return Window(label, None, c_win, True)


# ----------------------------------------------------------------- U
def _u_nos(mol, mf, c_pur, e_core):
    c_no, occ, groups, _ = build_nos(mol, mf, c_pur, e_core)
    return c_no, groups


def stage_u24() -> None:
    print(f"\n{'#' * 66}\n# U: NO-24 SA saturation + relax point (ticket 28)\n"
          f"{'#' * 66}")
    mol, mf, closed, c_pur, e_core = u_setup()
    c_no, groups = _u_nos(mol, mf, c_pur, e_core)
    windows = []
    for nk in (2, 4):
        nact = sum(len(groups[i]) for i in range(nk))
        print(f"NO window g1..g{nk}: {nact} cols")
        windows.append(_no_win(f"no{nk}", c_no[:, :nact]))
    spec = StudySpec(label="u_no", ncore=U_NCORE, na=NA, nb=NB,
                     windows=tuple(windows),
                     treatments=("sa-casscf", "relax"),
                     budgets={"sa-casscf": {"max_cycle_macro": 60}},
                     diagnostics="full")
    extra = (f"\nAnchors: SA-12 = -50.423810 (exact 5L); "
             f"relax-12 = -49.899729, relax-17 = -49.900366; "
             f"frozen NO-24 CASCI = -49.894287.\n")
    execute(spec, mol, mf, U_RES, c_pur, "U", extra)


def stage_six13() -> None:
    print(f"\n{'#' * 66}\n# U: CAS(6,13) 7s^2 SA (ticket 28)\n{'#' * 66}")
    from run_u_no import closed_shell_core_energy
    mol, mf, closed, c_pur, e_core = u_setup()
    c_no, groups = _u_nos(mol, mf, c_pur, e_core)

    # the 7s pur column: the l=0 column living ABOVE the 6p rows
    s = mf.get_ovlp()
    sect = drivers.lm_sectors(mol)
    l1max = max(r for (l, _m), rows in sect.items() if l == 1 for r in rows)
    s7rows = [r for (l, _m), rows in sect.items() if l == 0
              for r in rows if r > l1max]
    w7 = np.real(np.einsum("ri,rj,ij->j", c_pur[s7rows].conj(),
                           s[np.ix_(s7rows, s7rows)], c_pur[s7rows]))
    i7 = int(np.argmax(w7))
    assert w7[i7] > 0.8, (i7, w7)
    print(f"7s pur column: index {i7}, S-weight on the 7s rows "
          f"= {w7[i7]:.4f}")
    c_core6 = np.delete(c_pur, i7, axis=1)
    e_core6 = closed_shell_core_energy(mf, c_core6)
    print(f"E_core: 5-col pur {e_core:.6f} -> 4-col (7s out) {e_core6:.6f}"
          f"  (dE_bookkeeping = {(e_core - e_core6) * 1e3:+.1f} mEh)")

    c_win = np.hstack([c_no[:, :12], c_pur[:, [i7]]])
    print(f"six13 window: {c_win.shape[1]} cols")
    win = _no_win("d5f7+7s", c_win)
    spec = StudySpec(label="u_six13_sa", ncore=U_NCORE - 1, na=5, nb=1,
                     windows=(win,), treatments=("sa-casscf",),
                     budgets={"sa-casscf": {"max_cycle_macro": 60}},
                     diagnostics="full")
    extra = (f"\nCAS(6,13) = NO {{d5, f7}} + the 7s pur column; core = "
             f"pur minus 7s; E_core(4col) = {e_core6:.6f}; "
             f"CASCI verdict was inert (-0.58 mEh).\n")
    execute(spec, mol, mf, U_RES, c_core6, "U", extra)


STAGES = {"sc": lambda: stage_tm("sc"), "ti": lambda: stage_tm("ti"),
          "u24": stage_u24, "six13": stage_six13}

if __name__ == "__main__":
    for st in sys.argv[1:]:
        STAGES[st]()
