"""The CAS(6,13) V&V leg -- the 7s-fix leg of the flagship CAS-size
axis, sector (5,1).  The as-built window (``10f13d7s``, the incumbent
distillation from the reference's occupied block) abandons the
intended 5f3 6d1 7s2 configuration at the frozen level (the 7F
ground; the 7s column unoccupied at 3e-03 -- the companion leg
``u14e_vv_6_13.py``'s disclosure).  THIS leg carries the fix
candidate, MIXED-SOURCE: the f/d channels from the PRIOR-AXIS ``@r1``
density-block distillation (``basis.series_window`` against the rung
itself; measured captures 5f 0.93 / 6d 0.98 -- better than the
incumbent's 0.82 / 0.92), and the 7s channel as THE REFERENCE'S OWN
PAIR MO, PROJECTED ONTO ITS PURE-s SECTOR.

The s-channel mechanism, measured (2026-11-23, this leg's probes):
the promoted pair is the reference's occupied column 143 (occ 2.0 --
the 5th closed pair by character, exactly the (6,13) promotion), and
the declared 7s realization (the 18s/19s/20s bands) does NOT contain
it: the pair MO spreads over >= 6 s primitives (18s 0.524, 19s
0.506, 15s 0.366, 12s 0.327, 14s 0.303, 10s 0.255 raw coefficients
-- s-diffuse across the basis), the distilled top radial of the
declared span captures cos^2 = 0.7094 of it, and the single-band
pick fallback would capture less (the best band alone ~0.27).  No
band-limited construction can host the pair faithfully; the honest s
column is the pair's own orbital.  It is already core-orthogonal
(weight on the pur core span 6.5e-07 -- canonical S-orthonormality
against the folded pairs).  ONE correction is required: the C1
gauge's pair MO carries an m-mixed d tail (the l-sector weights
measured s 0.9873 / d 0.0127 -- fact 12's residue); the raw MO in
the window BREAKS the rotation closure (measured ||[P,Lperp]||
8.4e-02, the definite-L certificate lost at the tail's scale), so
the window takes the pair MO's PURE-s SECTOR (cos^2 0.9873 of the
pair kept) -- an l = 0 direction is L-invariant, the assembled span
stays rotation-closed.

THE ACCEPTANCE CONTROLS, all on the frozen CASCI: the 7s pair forms
(the s-column occupation ~2 -- the as-built read 3e-03); the ground
reads a definite L on a rotation-closed window (the closure meters
~1e-13); the frozen total sits near the reference (dE << +1065 mEh,
the as-built's frozen dE).  Then the free CASSCF / free SA (the
campaign budgets: max_cycle_macro 200 -- the optimizer-budget
directive) with PT2 -> the ``...|(6,13)@r1+MOs`` table rows.  The
manifold physics to expect: the (6,13) quasi-degenerate block (the
certified lane record: the declared 5L is not the block ground --
16 x ~90 + one x ~72); the frozen ground's class and pair
occupation are the disclosure.  NO pinned route here (the
unification ladder is the pinned business; step C).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_6_13_r1.py \
        [--ground-only]  |  tee logs/u14e_vv_6_13_r1.log
"""

import sys
import time

import numpy as np
from pyscf import mcscf
from pyscf.mcscf import addons as _mcaddons

from _vv_common import (Verdict, conv_of, core_part_of, drop_row,
                        fmt_range, free_kernel, load_reference, meters,
                        price, std_row)

KEY_613_R1 = "u14e_5f6d7s@r1_5-1_casci_ao_-_pur_avg"
SECTOR = (5, 1)
# the core PAIR count of this key: 14 electrons - (5+1) active = 8 in
# 4 pairs.  MUST equal key_core's ncore (asserted in production_window)
# -- the slice-incident guard (see the companion leg).
NCORE = 4
ROOTS = 17


def _reference_s_mo(mol, mf) -> tuple[int, np.ndarray]:
    """The reference's own 7s PAIR MO: among the occupied columns
    OUTSIDE the 4-pair core set, the one dominated by an s band --
    asserted unique.  It is the 5th closed pair by character (the
    (6,13) promotion promotes exactly this pair out of the core):
    occ 2.0, and the only occupied s-dominated column outside the
    folded prefix (the open columns are f/d-dominated)."""
    from atomterm.windows import occupied_core_idx
    core = set(occupied_core_idx(mf, NCORE))
    occ = [j for j in range(mf.mo_coeff.shape[1])
           if mf.mo_occ[j] > 0 and j not in core]
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    hits = []
    for j in occ:
        v = np.abs(np.asarray(mf.mo_coeff)[:, j])
        if labels[int(np.argmax(v))].endswith("s"):
            hits.append(j)
    assert len(hits) == 1, (f"the valence s column is not unique: "
                            f"occupied non-core columns {hits} are "
                            f"s-dominated -- the identification needs "
                            f"the character spans, not the labels")
    return hits[0], np.asarray(mf.mo_coeff)[:, hits[0]]


def production_window_r1(mol, mf, spec):
    """The mixed-source window, written out: ``key_core`` +
    ``window_columns`` (the ``@r1`` tag routes inside to
    ``basis.series_window`` with ``series_mfs=[mf]``, the rung as the
    head -- the f/d radials) + the s channel REPLACED by the
    reference's own 7s pair MO, PROJECTED ONTO ITS PURE-s SECTOR (the
    mechanism in the module docstring) + ``build_full_mo``.  The
    preamble: the key's core-pair identity, the s replacement's
    verification (the l-sector decomposition -- the m-mixed C1 d
    tail measured 1.27e-02 and removed; the cos^2 against the full
    MO ~0.987), the distilled columns' S-norms after removing the
    pur core span, and the CLOSURE METERS of the assembled window
    (the definite-L certificate).

    The shell layer here is the case's ``shell_bands`` (the f/d
    radials' three-band splits; ``prior_bands`` declares f/d alone
    but its realization coincides for these shells), so
    ``window_columns`` is called without ``prior_bands``."""
    import re as _re
    from atomterm.cells import CellKey, key_core, window_columns
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.operators import window_closure

    key = CellKey.from_address(KEY_613_R1)
    nc, c_core, c_pur = key_core(key, mol, mf, spec.shell_bands)
    assert nc == NCORE, (f"the key's core pairs {nc} != the leg's "
                         f"NCORE {NCORE} -- the meter slices would "
                         f"shift (the 2026-11-22 incident)")
    act, ncas, rep = window_columns(key, mol, mf, spec.shell_bands,
                                    spec.shell_radials, ncol=nc)
    S = mf.get_ovlp()
    # the composition disclosure: per shell, the distilled k of the
    # band chain, the discarded spectrum mass, the capture of the
    # head's significant columns
    if rep is not None:
        for sh, r in rep.items():
            caps = [c for c in r["capture"] if c]
            cap_s = ("capture n/a (no significant source columns)"
                     if not caps else
                     f"capture {min(min(c) for c in caps):.4f}"
                     f"-{max(max(c) for c in caps):.4f}")
            print(f"  [window] {sh}: k={r['k']} of {len(r['bands'])} "
                  f"bands, discarded {r['discarded']:.4f}, {cap_s}",
                  flush=True)
    norms = np.sqrt(np.diag(act.T @ S @ act)
                    - np.diag(act.T @ S @ c_pur @ np.linalg.solve(
                        c_pur.T @ S @ c_pur, c_pur.T @ S @ act)))
    lean = [(lab.split()[-1], float(n))
            for lab, n in zip(_labels_of(act, mol, S), norms)]
    print(f"  [diagnostic] distilled columns' S-norm after removing "
          f"the pur core span: "
          f"{', '.join(f'{n:.4f}' for _, n in lean)}", flush=True)

    # ---- the s channel: the reference's own 7s pair MO, s-pure -----
    j_s, mo_s = _reference_s_mo(mol, mf)
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    s_rows = [i for i, lab in enumerate(labels)
              if _re.match(r"\d+s$", lab)]
    d_rows = [i for i, lab in enumerate(labels)
              if _re.match(r"\d+d", lab)]
    v = mo_s
    w_s = float(v[s_rows] @ S[np.ix_(s_rows, s_rows)] @ v[s_rows])
    w_d = float(v[d_rows] @ S[np.ix_(d_rows, d_rows)] @ v[d_rows])
    j_raw = _s_column(act, mol, S)
    cos2_old = float((v @ S @ act[:, j_raw]) ** 2
                     / ((v @ S @ v)
                        * (act[:, j_raw] @ S @ act[:, j_raw])))
    mo_pure = np.zeros_like(v)
    mo_pure[s_rows] = v[s_rows]
    mo_pure = mo_pure / np.sqrt(float(mo_pure @ S @ mo_pure))
    cos2_pure = float((mo_pure @ S @ v) ** 2)
    print(f"  [s channel] the reference's 7s pair MO (occupied "
          f"column {j_s}, occ {mf.mo_occ[j_s]:.1f}) replaces the "
          f"distilled radial at {j_raw}: its l-sector weights s "
          f"{w_s:.6f} / d {w_d:.6f} (the m-mixed C1 tail -- fact "
          f"12's residue; PROJECTED OUT), the s-pure column keeps "
          f"cos^2 {cos2_pure:.6f} of the pair MO, whose cos^2 "
          f"against the distilled radial was {cos2_old:.4f} (the "
          f"declared 3-band span misses it)", flush=True)
    act = np.array(act, dtype=float)
    act[:, j_raw] = mo_pure

    mo_prod = build_full_mo(mf, occupied_core_idx(mf, NCORE), act,
                            c_core=c_pur)
    mo_win = mo_prod[:, NCORE:NCORE + ncas]
    cl = window_closure(mol, mo_win)
    ortho = float(np.max(np.abs(
        mo_win.T @ S @ mo_win - np.eye(ncas))))
    print(f"  [closure] assembled window: ||[P,Lz]|| "
          f"{cl['lz']:.2e}, ||[P,Lperp]|| {cl['lperp']:.2e} "
          f"(rotation closure -- the definite-L certificate)",
          flush=True)
    print(f"  [s channel] the assembled block's orthonormality "
          f"deviation {ortho:.2e}", flush=True)
    return mo_prod, ncas, cl, j_raw


def _labels_of(act, mol, S):
    """Each act column's dominant AO label (the diagnostic's naming)."""
    labs = []
    for j in range(act.shape[1]):
        w = np.abs(act[:, j])
        i = int(np.argmax(w))
        labs.append(mol.ao_labels()[i].split()[-1])
    return labs


def _s_column(act, mol, S):
    """The act block's s-column index (the 7s slot): the column whose
    dominant AO label is an s band -- never a hardcoded position."""
    labs = _labels_of(act, mol, S)
    hits = [j for j, lab in enumerate(labs) if lab.endswith("s")]
    assert len(hits) == 1, (f"the window's s sector is {len(hits)} "
                            f"columns ({[labs[j] for j in hits]}) -- "
                            f"the (6,13) anatomy is one s radial")
    return hits[0]


def _occupations(cas) -> np.ndarray:
    """The frozen root's active-space occupations (the CAS 1-RDM
    diagonal in the window basis -- PySCF's own construction)."""
    return np.diag(cas.make_rdm1())


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pt2_roots = [0] if ground_only else None
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    mo_full, ncas, cl, j_s = production_window_r1(mol, mf, spec)
    weights = [1.0 / 17] * 17
    expect_in = meters(mol, mo_full[:, NCORE:NCORE + ncas], *SECTOR)
    print(f"  [anatomy] the window's s column at index {j_s} of "
          f"{ncas}", flush=True)

    # ---- free CASCI (the frozen input window; the pricing class) -----
    t0 = time.time()
    cas_f = mcscf.CASCI(mf, ncas, SECTOR)
    cas_f.canonicalization = False
    cas_f.kernel(mo_full)
    r_f = price(cas_f, expect_in, roots=[0])[0]
    l_f = r_f["l2"]
    l_def = min(abs(l_f - lx * (lx + 1)) for lx in range(16))
    occ_s = float(_occupations(cas_f)[j_s])
    de_meh = (float(cas_f.e_tot) - mf.e_tot) * 1e3
    print(f"  [casci] e_tot {cas_f.e_tot:.8f}, <L2> {l_f:.6f}, "
          f"occ(7s) {occ_s:.4f}, dE {de_meh:+.1f} mEh "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.check("the @r1+MOs window is rotation-closed (the definite-L "
            "certificate)",
            cl["lz"] <= 1e-10 and cl["lperp"] <= 1e-10,
            f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
            f"{cl['lperp']:.2e} (the band: the as-built's observed "
            f"~1e-13)")
    v.check("the frozen ground reads a definite L",
            l_def <= 1e-4,
            f"<L2> {l_f:.6f} -- closest L(L+1) within {l_def:.2e} "
            f"(WHICH L is the disclosure)")
    v.check("the 7s pair forms (the fix's acceptance: the s "
            "occupation ~2; the as-built read 3e-03)",
            occ_s >= 1.5,
            f"occ(7s) = {occ_s:.4f}")
    v.check("the frozen total sits near the reference (dE << +1065, "
            "the as-built's frozen dE)",
            de_meh < 400.0,
            f"dE {de_meh:+.1f} mEh (the @r1 (4,12) class sits at "
            f"+149.3; the frozen (6,13) has one more shell's "
            f"variational freedom)")
    v.note("free CASCI (the frozen input window)",
           f"E(2) {r_f['e2'] * 1e3:+.4f} mEh, <S2> "
           f"{r_f['s2']:.4f}")
    drop_row(std_row("casci", "free", "(6,13)@r1+MOs",
                     float(cas_f.e_tot),
                     r_f["e2"] * 1e3,
                     float(cas_f.e_tot) + r_f["e2"],
                     f"{l_f:.4f}", f"{r_f['s2']:.4f}",
                     "definite-L (the 7s pair = the reference's own "
                     "MO column)",
                     conv=conv_of(cas_f.fcisolver.conv_tol, None)))

    # ---- free CASSCF (single root; no assertion) ----------------------
    t0 = time.time()
    mc_c, cas_c = free_kernel(mf, ncas, SECTOR, mo_full, track_iters=True)
    r_c = price(cas_c, expect_in, roots=[0])[0]
    occ_c = float(_occupations(cas_c)[j_s])
    print(f"  [casscf] e_tot {mc_c.e_tot:.8f}, converged "
          f"{mc_c.converged}, <L2> {r_c['l2']:.6f}, occ(7s) "
          f"{occ_c:.4f} ({time.time() - t0:.1f} s)", flush=True)
    v.note("free CASSCF (the single-root draw class)",
           f"<L2> {r_c['l2']:.6f}, occ(7s) {occ_c:.4f}, E(2) "
           f"{r_c['e2'] * 1e3:+.4f} mEh")
    drop_row(std_row("casscf", "free", "(6,13)@r1+MOs",
                     float(mc_c.e_tot),
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
    v.note("free SA (the draw-dependent meters)",
           f"<L2> in [{min(l2s_sa):.4f}, {max(l2s_sa):.4f}], <S2> up "
           f"to {max(r['s2'] for r in rows):.4f} -- the free-SA roots "
           f"partition the manifold draw-dependently; E(2) root 0 "
           f"{r0['e2'] * 1e3:+.4f} mEh, total "
           f"{r0['e_ref'] + r0['e2']:.4f}")
    # the max-<S2> root's spin class, by PySCF's own meter (whatever
    # class this manifold's high-spin member is -- no pre-judged
    # class here)
    k12 = max(range(len(rows)), key=lambda k: rows[k]["s2"])
    ss, mult = _mcaddons.spin_square(cas_sa, ci=cas_sa.ci[k12])
    s_of = (np.sqrt(4.0 * ss + 1.0) - 1.0) / 2.0
    s_at = (np.sqrt(4.0 * rows[k12]["s2"] + 1.0) - 1.0) / 2.0
    v.check("the free-SA max-<S2> root's spin class agrees between "
            "PySCF's spin_square and atomterm's meters (the two "
            "independent spin constructions)",
            abs(ss - rows[k12]["s2"]) <= 1e-6
            and abs(s_of - s_at) <= 1e-6,
            f"root {k12}: spin_square (S^2, 2S+1) = ({ss:.8f}, "
            f"{mult:.4f}) -> S = {s_of:.6f}; atomterm's meters read "
            f"<S2> {rows[k12]['s2']:.8f} -> S = {s_at:.6f}")
    print(f"  [sa] ensemble {mc_sa.e_tot:.8f}, converged "
          f"{mc_sa.converged} ({time.time() - t0:.1f} s)", flush=True)
    drop_row(std_row("sa-casscf", "free", "(6,13)@r1+MOs",
                     float(e_cas[0] + core_e), r0["e2"] * 1e3,
                     float(r0["e_ref"] + r0["e2"]),
                     fmt_range(l2s_sa),
                     f"{max(r['s2'] for r in rows):.4f}",
                     "draw-dependent partition",
                     conv=conv_of(mc_sa.final_de, mc_sa.n_macro,
                                  bool(mc_sa.converged))))

    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
