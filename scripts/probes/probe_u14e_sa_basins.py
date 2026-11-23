"""The (6,13)@r1 SA/CASSCF basin probe -- the valley map, printed
only (NO store writes; the user evaluates before any table/leg
change).

THE QUESTION (the user's review of the sa-casscf (free) block): the
single-root CASSCF descends the QUINTET valley (-50.4623, <L2> 72.08,
<S2> 6) while the free SA descends the SEPTET valley (-50.4529,
<L2> 90, <S2> 12 -- the exchange-maximal six-parallel class).  Which
wins -- the ensemble OBJECTIVE or the SEED?  Five runs:

  A  CASSCF from the input window      the quintet valley record
  B  free SA x17 from the input window the septet landing + the
                                       occupation FINGERPRINT of the
                                       L=9 root (the prediction:
                                       5f3(m=3,2,1) 6d2(m=2,1) 7s1(m=0),
                                       Sum m = 9; the alternative f3d3
                                       reads the s column ~0)
  C  free SA from A's orbitals         objective vs seed
  D  CASSCF from B's orbitals          does the single root escape to
                                       the quintet valley?
  E1 the pinned SA (the term's own     the 5L-block orbitals
     wfnsym-anchored ensemble)
  E2 free SA from E1's orbitals        the gauge-anchoring test: does
                                       the free SA STAY on the 5L
                                       block when started from the
                                       term's orbitals?

Every optimized run prints: E, <L2>/<S2> of the ground root, the
converged-window CLOSURE meters (the (a) measurement: the input
windows close at 1e-13; the defect is created by the optimization --
here measured, not inferred), the macro-cycle budget, and the active
OCCUPATIONS (the class fingerprint).

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/probes/probe_u14e_sa_basins.py \
        | tee logs/probe_u14e_sa_basins.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402
from pyscf.mcscf import addons as _mcaddons  # noqa: E402

from _vv_common import (core_part_of, free_kernel, load_reference,  # noqa: E402
                        make_blocks, meters, meters_all, pinned_sa_run)

from atomterm.operators import window_closure  # noqa: E402

SECTOR = (5, 1)
NCORE = 4
ROOTS = 17
BUDGET = {"max_cycle_macro": 200}


def _occ_print(tag, cas, ncas, root=0):
    """The active occupations of one root -- the class fingerprint."""
    ci = np.asarray(cas.ci[root] if isinstance(cas.ci, list)
                    else cas.ci).ravel()
    rdm = cas.fcisolver.make_rdm1(ci, ncas, SECTOR)
    occ = np.diag(rdm).real
    print(f"  [{tag}] occupations:", flush=True)
    print("    " + "  ".join(f"{i}:{o:.2f}" for i, o in enumerate(occ)),
          flush=True)
    return occ


def _ground_report(tag, mol, mo_win, ncas, cas, root=0):
    """<L2>/<S2> of one root + the converged-window closure."""
    expect = meters(mol, mo_win, *SECTOR)
    ci = np.asarray(cas.ci[root] if isinstance(cas.ci, list)
                    else cas.ci).ravel()
    l2 = expect("l2", ci)
    s2 = expect("s2", ci)
    cl = window_closure(mol, mo_win)
    ss, mult = _mcaddons.spin_square(cas, ci=cas.ci[root]
                                     if isinstance(cas.ci, list)
                                     else cas.ci)
    print(f"  [{tag}] ground root: <L2> {l2:.6f}, <S2> {s2:.6f} "
          f"(spin_square: {ss:.6f}, {mult:.1f}); the converged "
          f"window's closure: ||[P,Lz]|| {cl['lz']:.2e}, "
          f"||[P,Lperp]|| {cl['lperp']:.2e}", flush=True)
    return l2, s2, cl


def main() -> int:
    t00 = time.time()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    import u14e_vv_6_13_r1 as r1
    mo_full, ncas, cl_in, j_s = r1.production_window_r1(mol, mf, spec)
    print(f"  [window] the input closure ||[P,Lz]|| {cl_in['lz']:.2e}, "
          f"||[P,Lperp]|| {cl_in['lperp']:.2e}; the s column at "
          f"{j_s}", flush=True)
    weights = [1.0 / ROOTS] * ROOTS

    # ---- A: CASSCF from the input window ------------------------------
    print("===== A: CASSCF, input window =====", flush=True)
    t0 = time.time()
    mc_a, cas_a = free_kernel(mf, ncas, SECTOR, mo_full,
                              track_iters=True)
    print(f"  [A] E {mc_a.e_tot:.8f}, converged {mc_a.converged} "
          f"({mc_a.n_macro} macro, |dE| {mc_a.final_de:.1e}, "
          f"{time.time() - t0:.0f} s)", flush=True)
    _ground_report("A", mol, mc_a.mo_coeff[:, NCORE:NCORE + ncas],
                   ncas, cas_a)
    _occ_print("A", cas_a, ncas)

    # ---- B: free SA x17 from the input window -------------------------
    print("===== B: free SA x17, input window =====", flush=True)
    t0 = time.time()
    mc_b, cas_b = free_kernel(mf, ncas, SECTOR, mo_full,
                              weights=weights, nroots=ROOTS,
                              track_iters=True)
    core_e, e_cas = core_part_of(cas_b)
    print(f"  [B] ensemble {mc_b.e_tot:.8f}, converged "
          f"{mc_b.converged} ({mc_b.n_macro} macro, |dE| "
          f"{mc_b.final_de:.1e}, {time.time() - t0:.0f} s); root 0 "
          f"e_ref {e_cas[0] + core_e:.8f}", flush=True)
    _ground_report("B", mol, mc_b.mo_coeff[:, NCORE:NCORE + ncas],
                   ncas, cas_b)
    occ_b = _occ_print("B", cas_b, ncas)
    print(f"  [B] the fingerprint read: s column occ "
          f"{occ_b[j_s]:.3f}; f block sum {occ_b[:7].sum():.2f}; "
          f"d block sum {occ_b[7:12].sum():.2f}", flush=True)

    # ---- C: free SA from A's orbitals ---------------------------------
    print("===== C: free SA x17, started from A's orbitals =====",
          flush=True)
    t0 = time.time()
    mc_c, cas_c = free_kernel(mf, ncas, SECTOR, mc_a.mo_coeff,
                              weights=weights, nroots=ROOTS,
                              track_iters=True)
    print(f"  [C] ensemble {mc_c.e_tot:.8f}, converged "
          f"{mc_c.converged} ({mc_c.n_macro} macro, "
          f"{time.time() - t0:.0f} s)", flush=True)
    _ground_report("C", mol, mc_c.mo_coeff[:, NCORE:NCORE + ncas],
                   ncas, cas_c)

    # ---- D: CASSCF from B's orbitals ----------------------------------
    print("===== D: CASSCF, started from B's orbitals =====", flush=True)
    t0 = time.time()
    mc_d, cas_d = free_kernel(mf, ncas, SECTOR, mc_b.mo_coeff,
                              track_iters=True)
    print(f"  [D] E {mc_d.e_tot:.8f}, converged {mc_d.converged} "
          f"({mc_d.n_macro} macro, {time.time() - t0:.0f} s)", flush=True)
    _ground_report("D", mol, mc_d.mo_coeff[:, NCORE:NCORE + ncas],
                   ncas, cas_d)

    # ---- E1: the pinned SA (the term's own ensemble) -------------------
    print("===== E1: the pinned SA on the @r1 window =====", flush=True)
    import u14e_vv_unification as uni
    mol_d2h, mf_d2h = uni.wrapper_substrate(mol, mf, spec)
    act_s, ncas_s, purity = uni.act_s_mo_of(mol, mf, mol_d2h, spec)
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import CellKey, key_core
    key = CellKey.from_address(r1.KEY_613_R1)
    _nc, _cc, c_pur = key_core(key, mol, mf, spec.shell_bands)
    mo_d2h = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, _nc),
                           act_s, c_core=c_pur)
    t0 = time.time()
    out_e1, _ = pinned_sa_run(mol_d2h, mf_d2h, mo_d2h, _nc, ncas_s,
                              SECTOR, 8, "u", c_pur,
                              c_read=np.asarray(
                                  mo_d2h)[:, _nc:_nc + ncas_s],
                              mc_opts=BUDGET)
    blocks_e1, expect_e1, _ = make_blocks(mol_d2h, mf_d2h,
                                          np.asarray(out_e1["mo"]),
                                          _nc, ncas_s, SECTOR, 8, "u",
                                          c_pur, out_e1["roots"],
                                          c_read=np.asarray(
                                              out_e1["mo"])[:,
                                                            _nc:_nc + ncas_s])
    l2_e1, s2_e1 = meters_all(blocks_e1, expect_e1)
    print(f"  [E1] ensemble {out_e1['e_tot']:.8f}, converged "
          f"{out_e1['converged']} (cycles {out_e1['cycles']}, "
          f"{time.time() - t0:.0f} s); <L2> [{min(l2_e1):.6f}, "
          f"{max(l2_e1):.6f}], <S2> max {max(s2_e1):.6f}", flush=True)

    # ---- E2: free SA from the pinned SA's orbitals ---------------------
    print("===== E2: free SA x17, started from E1's orbitals =====",
          flush=True)
    t0 = time.time()
    mc_e2, cas_e2 = free_kernel(mf, ncas, SECTOR, out_e1["mo"],
                                weights=weights, nroots=ROOTS,
                                track_iters=True)
    print(f"  [E2] ensemble {mc_e2.e_tot:.8f}, converged "
          f"{mc_e2.converged} ({mc_e2.n_macro} macro, "
          f"{time.time() - t0:.0f} s)", flush=True)
    _ground_report("E2", mol, mc_e2.mo_coeff[:, NCORE:NCORE + ncas],
                   ncas, cas_e2)
    _occ_print("E2", cas_e2, ncas)

    # ---- the summary ---------------------------------------------------
    print("\n===== the valley map =====", flush=True)
    print(f"  A  casscf   (input seed)      E {mc_a.e_tot:.8f}",
          flush=True)
    print(f"  D  casscf   (B's seed)        E {mc_d.e_tot:.8f}",
          flush=True)
    print(f"  B  sa x17   (input seed)      ensemble {mc_b.e_tot:.8f}",
          flush=True)
    print(f"  C  sa x17   (A's seed)        ensemble {mc_c.e_tot:.8f}",
          flush=True)
    print(f"  E1 pinned sa                  ensemble "
          f"{out_e1['e_tot']:.8f}", flush=True)
    print(f"  E2 sa x17   (E1's seed)       ensemble "
          f"{mc_e2.e_tot:.8f}", flush=True)
    print(f"({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
