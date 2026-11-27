"""The CAS(6,13) V&V leg -- the 7s-promotion leg of the flagship
CAS-size axis, sector (5,1).  THE DISCLOSURE LEG: everything is
printed; only the certified-class check gates.

THE WINDOW, IN THE SCRIPT (the user's directive): the production
recipe's exact calls, written out here -- ``key_core`` +
``window_columns`` (the distilled span: the 10f/13d bands are the
identity; the 7s radial is the production distillation across its
three s bands) + ``build_full_mo`` -- followed by the leg's own
DIAGNOSTIC that exposes the assembly's behavior on this window: each
distilled column's S-norm after removing the pur core span (the 7s
is core-lean), and where the assembled window actually landed.

THE FIVE ROUTES: free CASCI / free CASSCF / free SA on the production
window; pinned CI / pinned SA on the lane window (the lane
construction).  THE WINDOW IS HEALTHY (the 2026-11-22 correction, two
sources): the assembled span IS the projected distilled span exactly
(principal cosines all 1 against (I-P)act), the closure meters read
~1e-13 (the campaign cell's win_lz/win_lperp agree) -- the earlier
"L-mixed 18.0085" was an ARTIFACT of this leg's shifted meter slice
(NCORE=5 vs the key's ncore=4; the kernels themselves always ran the
true window: PySCF derives ncore=4 from the sector).  THE REAL
DISCLOSURE: the frozen ground is the 7F (L=3, S=3 definite -- the
d5 high-spin class + one f, the 7s column UNOCCUPIED at 3e-03): the
as-built window's variational ground abandons the intended
configuration -- the weak distillation captures (10f 0.82, 7s
0.27-0.71) drive it.  THE CERTIFIED-CLASS CHECK: the pinned
components read the L=9 block ground (16 at ~90 + one at ~72) -- the
declared 5L is not this block's ground; the free-SA meters are
draw-dependent on this flat near-degenerate manifold and are printed,
not asserted.  THE SPIN CROSS-CHECK: the <S2> = 12 roots are verified
as S = 3 by PySCF's own spin_square -- an independent spin
construction agreeing with atomterm's meters (the electron count is
right; the septet is the all-parallel d5+f1 class, not a 7s-broken
pair).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_6_13.py \
        [--ground-only]  |  tee logs/u14e_vv_6_13.log
"""

import sys
import time

import numpy as np
from pyscf import mcscf
from pyscf.mcscf import addons as _mcaddons

from _vv_common import (BAND_CLASS, L_L2, Verdict, conv_of, core_part_of,
                        drop_row, fmt_range, free_kernel, load_reference,
                        make_blocks, meters, meters_all, nevpt_rows,
                        pinned_ci_run, pinned_sa_run, price, std_row)

KEY_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"
KEY_613_LANE = "u14e_5f6d7s@6s6p_5-1_sa-casscf_ao_-_pur_avg"
SECTOR = (5, 1)
# the core PAIR count of this key: 14 electrons - (5+1) active = 8 in
# 4 pairs.  MUST equal key_core's ncore (asserted in production_window)
# -- the 2026-11-22 incident: hardcoded 5 here shifted every leg meter
# slice by one (the kernels were unaffected; PySCF derives ncore=4
# itself), and the shifted-slice <L2> read the artifact 18.0085.
NCORE = 4
L, PARITY = 8, "u"
ROOTS = 17


def production_window(mol, mf, spec):
    """The production recipe, written out (``run_cell``'s exact three
    calls) + the leg's preamble: the key's core-pair identity (the
    assert that guards the slice incident), each distilled column's
    S-norm after removing the pur core span, the assembled window's
    placement against act, and the CLOSURE METERS of the assembled
    window (rotation closure is the frozen ground's definite-L
    certificate -- measured every run)."""
    from atomterm.cells import CellKey, key_core, window_columns
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.operators import window_closure

    key = CellKey.from_address(KEY_613)
    nc, c_core, c_pur = key_core(key, mol, mf, spec.shell_bands)
    assert nc == NCORE, (f"the key's core pairs {nc} != the leg's "
                         f"NCORE {NCORE} -- the meter slices would "
                         f"shift (the 2026-11-22 incident)")
    act, ncas, _rep = window_columns(key, mol, mf, spec.shell_bands,
                                     spec.shell_radials, ncol=nc)
    S = mf.get_ovlp()
    # the diagnostic: the distilled columns against the pur core span
    norms = np.sqrt(np.diag(act.T @ S @ act)
                    - np.diag(act.T @ S @ c_pur @ np.linalg.solve(
                        c_pur.T @ S @ c_pur, c_pur.T @ S @ act)))
    lean = [(lab.split()[-1], float(n))
            for lab, n in zip(_labels_of(act, mol, S), norms)]
    print(f"  [diagnostic] distilled columns' S-norm after removing "
          f"the pur core span: "
          f"{', '.join(f'{n:.4f}' for _, n in lean)}", flush=True)
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, NCORE), act,
                            c_core=c_pur)
    mo_win = mo_prod[:, NCORE:NCORE + ncas]
    cl = window_closure(mol, mo_win)
    print(f"  [closure] assembled window: ||[P,Lz]|| "
          f"{cl['lz']:.2e}, ||[P,Lperp]|| {cl['lperp']:.2e} "
          f"(rotation closure -- the frozen ground's definite-L "
          f"certificate)", flush=True)
    O = np.abs(mo_win.T @ S @ act)
    place = int(np.argmax(O.max(axis=0) > 0.95)) if \
        (O.max(axis=0) > 0.95).any() else -1
    print(f"  [diagnostic] the assembled window block "
          f"[{NCORE}:{NCORE + ncas}] matches act's columns from "
          f"index {place} onward (the placement the production "
          f"assembly chose)", flush=True)
    if norms.min() < 1e-1:
        print(f"  [diagnostic] FINDING: a distilled direction is "
              f"core-lean (min norm {norms.min():.4f})", flush=True)
    return mo_prod, ncas


def _labels_of(act, mol, S):
    """Each act column's dominant AO label (the diagnostic's naming)."""
    labs = []
    for j in range(act.shape[1]):
        w = np.abs(act[:, j])
        i = int(np.argmax(w))
        labs.append(mol.ao_labels()[i].split()[-1])
    return labs


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pt2_roots = [0] if ground_only else None
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    mo_full, ncas = production_window(mol, mf, spec)
    weights = [1.0 / 17] * 17
    expect_in = meters(mol, mo_full[:, NCORE:NCORE + ncas], *SECTOR)

    # ---- free CASCI (the frozen input window; the pricing class) -----
    t0 = time.time()
    cas_f = mcscf.CASCI(mf, ncas, SECTOR)
    cas_f.canonicalization = False
    cas_f.kernel(mo_full)
    r_f = price(cas_f, expect_in, roots=[0])[0]
    l_f = r_f["l2"]
    l_def = min(abs(l_f - lx * (lx + 1)) for lx in range(16))
    print(f"  [casci] e_tot {cas_f.e_tot:.8f}, <L2> {l_f:.6f} "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("free CASCI (the frozen input window)",
           f"<L2> {l_f:.6f} -- closest L(L+1) within {l_def:.2e} "
           f"(the window is rotation-closed, so the ground has a "
           f"definite L; WHICH L is the disclosure), E(2) "
           f"{r_f['e2'] * 1e3:+.4f} mEh")
    drop_row(std_row("casci", "free@pur", "(6,13)", float(cas_f.e_tot),
                     r_f["e2"] * 1e3,
                     float(cas_f.e_tot) + r_f["e2"],
                     f"{l_f:.4f}", f"{r_f['s2']:.4f}",
                     "definite-L: the 7F, S=3 (the d5 high-spin "
                     "ground; the 7s column unoccupied)",
                     conv=conv_of(cas_f.fcisolver.conv_tol, None)))

    # ---- free CASSCF (single root; no assertion) ----------------------
    t0 = time.time()
    mc_c, cas_c = free_kernel(mf, ncas, SECTOR, mo_full, track_iters=True)
    r_c = price(cas_c, expect_in, roots=[0])[0]
    print(f"  [casscf] e_tot {mc_c.e_tot:.8f}, converged "
          f"{mc_c.converged}, <L2> {r_c['l2']:.6f} "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("free CASSCF (the single-root draw class)",
           f"<L2> {r_c['l2']:.6f}, E(2) {r_c['e2'] * 1e3:+.4f} mEh")
    drop_row(std_row("casscf", "free@pur", "(6,13)", float(mc_c.e_tot),
                     r_c["e2"] * 1e3,
                     float(mc_c.e_tot) + r_c["e2"],
                     f"{r_c['l2']:.4f}", f"{r_c['s2']:.4f}",
                     "escape draw",
                     conv=conv_of(mc_c.final_de, mc_c.n_macro,
                                  bool(mc_c.converged))))

    # ---- free SA (the draw-dependent meters, printed) -----------------
    t0 = time.time()
    mc_sa, cas_sa = free_kernel(mf, ncas, SECTOR, mo_full,
                                weights=weights, nroots=ROOTS,
                                track_iters=True)
    core_e, e_cas = core_part_of(cas_sa)
    expect_sa = meters(mol, cas_sa.mo_coeff[:, NCORE:NCORE + ncas],
                       *SECTOR)
    rows = price(cas_sa, expect_sa,
                 roots=range(1) if ground_only else range(ROOTS))
    for r in rows:
        r["e_ref"] = float(e_cas[r["root"]] + core_e)
        print(f"  root {r['root']}: e_ref {r['e_ref']:.8f}, E(2) "
              f"{r['e2'] * 1e3:+.4f} mEh, <L2> {r['l2']:.6f}, "
              f"<S2> {r['s2']:.6f}", flush=True)
    r0 = rows[0]
    l2s_sa = [r["l2"] for r in rows]
    at_l9 = max(abs(x - 90.0) for x in l2s_sa) <= BAND_CLASS
    l9_class = "the L=9 class (90), S=3" if at_l9 else "mixed classes"
    v.note("free SA (the draw-dependent meters)",
           f"<L2> in [{min(l2s_sa):.4f}, {max(l2s_sa):.4f}], <S2> up "
           f"to {max(r['s2'] for r in rows):.4f} -- the free-SA roots "
           f"partition this flat near-degenerate manifold "
           f"draw-dependently (this draw: {l9_class}); E(2) root 0 "
           f"{r0['e2'] * 1e3:+.4f} mEh, total "
           f"{r0['e_ref'] + r0['e2']:.4f} (the prior committed "
           f"run's rows read -50.732 -- git history)")
    print(f"  [sa] ensemble {mc_sa.e_tot:.8f}, converged "
          f"{mc_sa.converged} ({time.time() - t0:.1f} s)", flush=True)
    drop_row(std_row("sa-casscf", "free@pur", "(6,13)",
                     float(e_cas[0] + core_e), r0["e2"] * 1e3,
                     float(r0["e_ref"] + r0["e2"]),
                     fmt_range(l2s_sa),
                     f"{max(r['s2'] for r in rows):.4f}",
                     l9_class,
                     conv=conv_of(mc_sa.final_de, mc_sa.n_macro,
                                  bool(mc_sa.converged))))

    # ---- the <S2> = 12 class: PySCF's OWN spin meter ------------------
    # the (5,1) sector hosts S = 3 (6 electrons all parallel).  The
    # frozen occupations MEASURED (2026-11-22): the d5 high-spin class
    # (five d columns at ~0.99) + one f (10f+0 at 0.998) -- the 7s
    # column UNOCCUPIED (3e-03): the window's variational ground
    # abandons the intended 5f3 6d1 7s2 configuration, the weak
    # distillation captures (7s 0.27-0.71) driving it.  spin_square is
    # an independent construction of S^2 (fci.spin_op algebra vs
    # atomterm's operator application); its second return is 2S+1, so
    # S = 3 reads (12, 7)
    k12 = max(range(len(rows)), key=lambda k: rows[k]["s2"])
    ss, mult = _mcaddons.spin_square(cas_sa, ci=cas_sa.ci[k12])
    s_of = (np.sqrt(4.0 * ss + 1.0) - 1.0) / 2.0
    s_at = (np.sqrt(4.0 * rows[k12]["s2"] + 1.0) - 1.0) / 2.0
    v.check("the free-SA <S2>=12 root is S = 3 by PySCF's own "
            "spin_square (the two independent spin constructions "
            "agree)",
            abs(ss - rows[k12]["s2"]) <= 1e-6 and abs(s_of - 3.0) <= 1e-6,
            f"root {k12}: spin_square (S^2, 2S+1) = ({ss:.8f}, "
            f"{mult:.4f}) -> S = {s_of:.6f}; atomterm's meters read "
            f"<S2> {rows[k12]['s2']:.8f} -> S = {s_at:.6f} -- the "
            f"all-parallel d5+f1 septet class, the electron count "
            f"RIGHT")

    # ---- the pinned route (re-derives the certified class) ------------
    from _vv_common import lane_route
    lmol, lmf, c_core_l, nf, mo_l, ncas_l, _act_l = lane_route(
        spec, KEY_613_LANE, SECTOR)

    # pinned CI
    t0 = time.time()
    roots_ci, _ = pinned_ci_run(lmol, lmf, mo_l, nf, ncas_l, SECTOR,
                                L, PARITY, c_core_l)
    blocks, expect_p, _ = make_blocks(lmol, lmf, mo_l, nf, ncas_l,
                                      SECTOR, L, PARITY, c_core_l,
                                      roots_ci)
    l2s, s2s = meters_all(blocks, expect_p)
    e2s = [row.e_corr for _, mc in blocks
           for row in nevpt_rows(mc, pt2_roots)]
    term_means = []
    for _, mc in blocks:
        ec = np.atleast_1d(mc.e_cas).ravel()
        cb = float(mc.e_tot) - float(ec[0])
        term_means.append(float(np.mean(ec + cb)))
    e_var_ci = float(np.mean(term_means))
    print(f"  [pinned casci] term e_var {e_var_ci:.8f}, <L2> "
          f"[{min(l2s):.6f}, {max(l2s):.6f}] ({time.time() - t0:.1f} s)",
          flush=True)
    v.note("pinned casci", f"<L2> in [{min(l2s):.4f}, {max(l2s):.4f}], "
           f"<S2> max {max(s2s):.4f}; E(2) mean "
           f"{np.mean(e2s) * 1e3:+.4f} mEh")
    # the lane-core route is outside the table's keep-set (the
    # pinned@pur rows are the table's pinned route): disclosed only
    v.note("the lane-pinned casci (print-only)",
           f"E_method {e_var_ci:.8f}, E(2) mean "
           f"{np.mean(e2s) * 1e3:+.4f} mEh")

    # pinned SA + THE CERTIFIED-CLASS CHECK
    t0 = time.time()
    out_sa, _ = pinned_sa_run(lmol, lmf, mo_l, nf, ncas_l, SECTOR,
                              L, PARITY, c_core_l)
    blocks_s, expect_s, _ = make_blocks(lmol, lmf,
                                        np.asarray(out_sa["mo"]),
                                        nf, ncas_l, SECTOR, L, PARITY,
                                        c_core_l, out_sa["roots"])
    l2s_s, s2s_s = meters_all(blocks_s, expect_s)
    e2s_s = [row.e_corr for _, mc in blocks_s
             for row in nevpt_rows(mc, pt2_roots)]
    near90 = sum(abs(x - (90.0)) <= BAND_CLASS for x in l2s_s)
    near72 = sum(abs(x - L_L2) <= BAND_CLASS for x in l2s_s)
    print(f"  [pinned sa] ensemble {out_sa['e_tot']:.8f}, <L2> "
          f"[{min(l2s_s):.6f}, {max(l2s_s):.6f}] "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.check("the pinned (6,13) components read the certified L=9 "
            "block ground (16 at ~90 + one at ~72)",
            near90 == 16 and near72 == 1,
            f"{len(l2s_s)} components: {near90} in 90 +- "
            f"{BAND_CLASS:g}, {near72} in 72 +- {BAND_CLASS:g} "
            f"(the declared 5L is not this block's ground -- the "
            f"route-wave record's classification, re-derived)")
    # the L=9 block ground's spin class, by PySCF's own meter
    # (the flat component index -> its block + component, the same
    # walk meters_all uses)
    k_s = max(range(len(s2s_s)), key=lambda k: s2s_s[k])
    acc, bmc, bci = 0, None, None
    for _, mc in blocks_s:
        n = len(np.atleast_1d(mc.e_cas))
        if acc + n > k_s:
            bmc, bci = mc, k_s - acc
            break
        acc += n
    ss_p, mult_p = _mcaddons.spin_square(bmc, ci=bmc.ci[bci])
    v.note("the pinned block ground's spin class (PySCF spin_square)",
           f"max-<S2> component: (S^2, 2S+1) = ({ss_p:.8f}, "
           f"{mult_p:.4f}) vs atomterm's meters {max(s2s_s):.8f}")
    v.note("the lane-pinned SA (print-only)",
           f"ensemble {out_sa['e_tot']:.8f}, E(2) mean "
           f"{np.mean(e2s_s) * 1e3:+.4f} mEh")

    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
