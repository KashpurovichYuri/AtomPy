"""The unification ladder -- the pinned route on the PUR core through
the D2h-tagged C1 wrapper: NO D2h SCF anywhere.  The legitimacy math
(the wrapper's contract, C0): at self-consistency the density is
G-invariant, hence [F, R] = 0 for every group element R; the Fock
coupling between irrep blocks obeys F_kl = chi_k(R)* chi_l(R) F_kl,
so F_kl = 0 for k != l by character orthogonality -- the Roothaan
problem splits into per-block problems and every converged MO lives
in exactly one block.  The CI consumes AO integrals (gauge-independent)
plus given columns: the C1 rung supplies both; the D2h-tagged mol
carries the symmetry framework (the block projectors, the
``wfnsym``-pinned solvers).  A symmetry-BROKEN solution is not even
representable in the framework -- the seed only picks the basin, and
here no SCF runs at all.

THE LADDER (the approved plan):
  C0  the wrapper probe: the D2h clone's AO quantities IDENTICAL to
      the C1 mol's (measured 0.00e+00), the CASCI on the hybrid
      reproducing the C1-mol value (and the campaign's frozen @r1
      record -50.2325), the window columns already irrep-pure
      (measured 2.2e-16).
  C1  THE GATE: pin-on-pur (4,12) -- the frozen CI must read the
      certified lane record 17 x 72.000000 +- band, and the SA the
      5L ensemble.  FAIL = stop and report, the lane stays; HOLD =
      the lane-retirement decision is recorded for the user.
  C2  pin-on-pur (4,15): the 15-dim core-orthogonal irrep-pure
      window on the PUR core (the 14p survival 0.28 -- the lane
      core's span contained the slots, the pur core does not) ->
      pinned CI + pinned SA (the 200-macro-cycle budget) + PT2 ->
      rows.  New information: the class disclosure is the point.
  C3  pin-on-pur (6,13): the class re-derivation on the as-built
      distilled window (the lane's certified class: 16 x ~90 +
      one x ~72) and on step B's @r1 window (disclosed).
  C4  the pinned SA on the LEAN LANE window at (4,15) (the approved
      cheap first pinned-(4,15) SA; the lane construction has no
      purify_window refusal -- the vanish class was the PUR route's
      core treatment; PySCF-side refusals are disclosed).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_unification.py \
        [--ground-only]  |  tee logs/u14e_vv_unification.log
"""

import sys
import time

import numpy as np
from pyscf import gto, mcscf, scf

from vv.meters import meters_all
from vv.reference import load_reference
from vv.routes import BUDGET_200, lane_route, make_blocks, nevpt_rows, pinned_ci_run, pinned_sa_run
from vv.table import conv_of, drop_row, fmt_range, std_row
from vv.verdict import Verdict

ADDR_412_R1 = "u14e_5f6d@r1_4-0_casci_ao_-_pur_avg"
ADDR_415 = "u14e_10f13d14p_4-0_casci_ao_-_pur_avg"
ADDR_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"
ADDR_415_LANE = "u14e_10f13d14p@6s6p7s_4-0_sa-casscf_ao_-_pur_avg"
L, PARITY = 8, "u"
BAND_GATE = 1e-5   # the lane record's certified band (observed dev
                   # 1.9e-07 there; the wrapper's fresh read 1.1e-13)
BUDGET = BUDGET_200  # the optimizer-budget directive (from vv.routes)


def wrapper_substrate(mol_c1, mf_c1, spec):
    """The D2h-tagged wrapper: the clone mol (same basis/ECP/spin,
    ``symmetry="D2h"``) + an SCF object rebuilt around it with the C1
    rung's record spliced in -- NO SCF runs.  Returns (mol, mf)."""
    mol_d2h = gto.M(atom=f"{spec.symbol} 0 0 0", basis=spec.basis,
                    spin=spec.spin, ecp=spec.ecp, symmetry="D2h",
                    verbose=0)
    mf_d2h = scf.ROHF(mol_d2h)
    mf_d2h.__dict__.update({"mo_coeff": mf_c1.mo_coeff,
                            "mo_occ": mf_c1.mo_occ,
                            "e_tot": mf_c1.e_tot, "converged": True})
    return mol_d2h, mf_d2h


def window_of(key_addr, mol_c1, mf_c1, spec, mf_d2h, radials_ext=None):
    """One key's window, built on the C1 rung (the construction's AO
    quantities are gauge-independent), assembled on the wrapper.  The
    core: the key's own pur resolution (``key_core``).
    ``radials_ext`` -- the leg's own k declarations for shells the
    case's ``shell_radials`` does not carry (the new-window
    compositions declare their additions here, never in src).
    Returns (mo, ncas, c_pur, nf)."""
    from atomterm.cells import CellKey, key_core, window_columns
    from atomterm.windows import build_full_mo, occupied_core_idx
    key = CellKey.from_address(key_addr)
    nc, _c_core, c_pur = key_core(key, mol_c1, mf_c1, spec.shell_bands)
    radials = ({**spec.shell_radials, **radials_ext}
               if radials_ext else spec.shell_radials)
    act, ncas, _ = window_columns(key, mol_c1, mf_c1, spec.shell_bands,
                                  radials, ncol=nc)
    mo = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nc), act,
                       c_core=c_pur)
    return mo, ncas, c_pur, nc


def act_s_mo_of(mol_c1, mf_c1, mol_d2h, spec):
    """Step B's (6,13) act block (PRE-assembly): the ``@r1`` f/d
    radials with the s column REPLACED by the reference's own 7s pair
    MO, projected onto its pure-s sector (the construction of the
    ``u14e_vv_6_13_r1`` leg -- the module docstring there owns the
    mechanism), then S-orthonormalized WITHIN each dominant-irrep
    group (:func:`orthonormal_per_irrep` on the D2h WRAPPER mol --
    the projectors read ``symm_orb``, which only the symmetry clone
    carries; the AO basis and overlap are identical) -- the replaced
    column overlaps the f/d columns at ~1e-2, and the assembly's
    global Loewdin would mix the pure columns across irrep blocks;
    the per-group treatment preserves the purity the pinned solvers'
    labeling consumes.  Returns (act, ncas, purity_loss)."""
    import re as _re
    import u14e_vv_6_13_r1 as r1
    from vv.symmetry import orthonormal_per_irrep
    from atomterm.cells import CellKey, window_columns
    key = CellKey.from_address(r1.KEY_613_R1)
    act, ncas, _ = window_columns(key, mol_c1, mf_c1, spec.shell_bands,
                                  spec.shell_radials, ncol=r1.NCORE)
    _j, mo_s = r1._reference_s_mo(mol_c1, mf_c1)
    labels = [lab.split()[-1] for lab in mol_c1.ao_labels()]
    s_rows = [i for i, lab in enumerate(labels)
              if _re.match(r"\d+s$", lab)]
    mo_pure = np.zeros_like(mo_s)
    mo_pure[s_rows] = mo_s[s_rows]
    mo_pure = mo_pure / np.sqrt(float(mo_pure @ mf_c1.get_ovlp()
                                      @ mo_pure))
    act = np.array(act, dtype=float)
    act[:, r1._s_column(act, mol_c1, mf_c1.get_ovlp())] = mo_pure
    act, worst = orthonormal_per_irrep(mol_d2h, mf_c1.get_ovlp(), act)
    return act, ncas, worst


def pinned_section(v, tag, mol_d2h, mf_d2h, mo, nf, ncas, c_core,
                   sector, pt2_roots, sa=True, ci_gate=None):
    """One pinned CI (+ SA) derivation, shared by the ladder's
    steps: the block grounds + meters, the PT2 pricing, the rows.
    ``ci_gate`` -- an optional callable (l2s) -> (name, ok, detail)
    appended to the verdict after the CI meters (C1's certified-class
    check).  Returns (ci_meters, sa_result_or_None, sa_meters)."""
    t0 = time.time()
    roots_ci, _ = pinned_ci_run(mol_d2h, mf_d2h, mo, nf, ncas, sector,
                                L, PARITY, c_core)
    blocks, _expect, _ = make_blocks(mol_d2h, mf_d2h, mo, nf, ncas,
                                     sector, L, PARITY, c_core, roots_ci)
    l2s, s2s = meters_all(blocks, _expect)
    print(f"  [{tag} casci] {len(l2s)} components, <L2> in "
          f"[{min(l2s):.6f}, {max(l2s):.6f}] (max dev from 72 "
          f"{max(abs(x - 72.0) for x in l2s):.2e}), <S2> max "
          f"{max(s2s):.6f} ({time.time() - t0:.1f} s)", flush=True)
    if ci_gate is not None:
        name, ok, detail = ci_gate(l2s)
        v.check(name, ok, detail)
    # the frozen CI's PT2 + row (the term energy = the components'
    # mean plus the core bookkeeping -- the pinned casci convention)
    e2s = [row.e_corr for _, mc in blocks
           for row in nevpt_rows(mc, pt2_roots)]
    term_means = []
    for _, mc in blocks:
        ec = np.atleast_1d(mc.e_cas).ravel()
        cb = float(mc.e_tot) - float(ec[0])
        term_means.append(float(np.mean(ec + cb)))
    e_var_ci = float(np.mean(term_means))
    v.note(f"{tag} casci", f"E(2) mean {np.mean(e2s) * 1e3:+.4f} mEh "
           f"over the priced roots")
    drop_row(std_row("casci", "pinned@pur", tag, e_var_ci,
                     float(np.mean(e2s) * 1e3),
                     float(e_var_ci + np.mean(e2s)),
                     fmt_range(l2s), f"{max(s2s):.4f}",
                     "block grounds" if ci_gate is None
                     else "the certified class",
                     conv=conv_of(1e-12, None)))
    if not sa:
        return (l2s, None, None)
    t0 = time.time()
    out_sa, _ = pinned_sa_run(mol_d2h, mf_d2h, mo, nf, ncas, sector,
                              L, PARITY, c_core, mc_opts=BUDGET)
    blocks_s, _exp_s, _ = make_blocks(mol_d2h, mf_d2h,
                                      np.asarray(out_sa["mo"]),
                                      nf, ncas, sector, L, PARITY,
                                      c_core, out_sa["roots"])
    l2s_s, s2s_s = meters_all(blocks_s, _exp_s)
    print(f"  [{tag} sa] ensemble {out_sa['e_tot']:.8f}, converged "
          f"{out_sa['converged']} (cycles {out_sa['cycles']}), <L2> "
          f"[{min(l2s_s):.6f}, {max(l2s_s):.6f}] "
          f"({time.time() - t0:.1f} s)", flush=True)
    e2s_s = [row.e_corr for _, mc in blocks_s
             for row in nevpt_rows(mc, pt2_roots)]
    drop_row(std_row("sa-casscf", "pinned@pur", tag,
                     float(out_sa["e_tot"]),
                     float(np.mean(e2s_s) * 1e3),
                     float(out_sa["e_tot"] + np.mean(e2s_s)),
                     fmt_range(l2s_s), f"{max(s2s_s):.4f}",
                     "draw-dependent partition",
                     conv=conv_of(out_sa["conv_dE"], out_sa["cycles"],
                                  bool(out_sa["converged"]))))
    return (l2s, out_sa, l2s_s)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pt2_roots = [0] if ground_only else None
    v = Verdict()
    mol_c1, mf_c1 = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    mol_d2h, mf_d2h = wrapper_substrate(mol_c1, mf_c1, spec)

    # ---- C0: the wrapper probe ----------------------------------------
    S1, S2 = mf_c1.get_ovlp(), mf_d2h.get_ovlp()
    h1, h2 = mf_c1.get_hcore(), mf_d2h.get_hcore()
    dev_s = float(np.abs(S1 - S2).max())
    dev_h = float(np.abs(h1 - h2).max())
    v.check("C0: the wrapper's AO quantities are the C1 mol's "
            "(gauge-independent integrals)",
            mol_d2h.nao == mol_c1.nao and dev_s == 0.0 and dev_h == 0.0,
            f"nao {mol_c1.nao}/{mol_d2h.nao}, ovlp dev {dev_s:.2e}, "
            f"hcore dev {dev_h:.2e}")
    from atomterm.cells import CellKey, key_core, window_columns
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.irreps import orb_irreps
    from vv.symmetry import symmetrize_columns
    mo412, ncas412, c_pur412, nf412 = window_of(
        ADDR_412_R1, mol_c1, mf_c1, spec, mf_d2h)
    key = CellKey.from_address(ADDR_412_R1)
    nc_c1, _, _ = key_core(key, mol_c1, mf_c1, spec.shell_bands)
    _act, _, _rep = window_columns(key, mol_c1, mf_c1, spec.shell_bands,
                                   spec.shell_radials, ncol=nc_c1)
    mo412_c1 = build_full_mo(mf_c1, occupied_core_idx(mf_c1, nc_c1),
                             _act, c_core=c_pur412)
    cas_c1 = mcscf.CASCI(mf_c1, ncas412, (4, 0))
    cas_c1.canonicalization = False
    cas_c1.kernel(mo412_c1)
    cas_d = mcscf.CASCI(mf_d2h, ncas412, (4, 0))
    cas_d.canonicalization = False
    cas_d.kernel(mo412)
    dev_e = abs(float(cas_c1.e_tot) - float(cas_d.e_tot))
    v.check("C0: the CASCI on the hybrid reproduces the C1-mol value "
            "(and the campaign's frozen @r1 record)",
            dev_e <= 1e-6 and abs(float(cas_d.e_tot) + 50.2325075273)
            <= 1e-6,
            f"C1 {cas_c1.e_tot:.10f} | wrapper {cas_d.e_tot:.10f} "
            f"| dev {dev_e:.2e} (the campaign's frozen row: "
            f"-50.2325, +149.3 mEh)")
    c_act = np.asarray(mo412)[:, nf412:nf412 + ncas412]
    _fixed, imp = symmetrize_columns(mol_d2h, S2, c_act, tol=1.0)
    v.check("C0: the distilled window columns are irrep-pure on the "
            "wrapper (no cleanup needed; the labels resolve)",
            max(imp) <= 1e-6,
            f"max impurity {max(imp):.2e}; labels "
            f"{orb_irreps(mol_d2h, c_act).tolist()}")

    # ---- C1: THE GATE (pin-on-pur (4,12)) ------------------------------

    def _gate_412(l2s):
        n72 = sum(abs(x - 72.0) <= BAND_GATE for x in l2s)
        ok = n72 == len(l2s)
        detail = (f"{n72}/{len(l2s)} components at 72.000000 "
                  f"+- {BAND_GATE:g} (max dev "
                  f"{max(abs(x - 72.0) for x in l2s):.2e}) -- the "
                  f"certified lane record reproduced through the "
                  f"unified route")
        return ("C1 GATE: the pin-on-pur (4,12) frozen CI reads the "
                "certified lane record 17 x 72.000000", ok, detail)

    def _gate_412_sa(l2s_s):
        n72 = sum(abs(x - 72.0) <= BAND_GATE for x in l2s_s)
        ok = n72 == len(l2s_s)
        detail = (f"{n72}/{len(l2s_s)} SA components at 72.000000 "
                  f"+- {BAND_GATE:g} -- the 5L ensemble, the lane "
                  f"record's second half")
        return ("C1 GATE: the pin-on-pur (4,12) SA reads the 5L "
                "ensemble", ok, detail)

    print("  [C1] pin-on-pur (4,12): the gate", flush=True)
    l2s_ci, out_sa_412, l2s_sa = pinned_section(
        v, "(4,12)", mol_d2h, mf_d2h, mo412, nf412, ncas412,
        c_pur412, (4, 0), pt2_roots, sa=True, ci_gate=_gate_412)
    if out_sa_412 is not None:
        name, ok, detail = _gate_412_sa(l2s_sa)
        v.check(name, ok, detail)
    if not all(ok for _, ok, kind in v.rows if kind == "check"):
        print("\n===== THE GATE FAILED -- stop and report; the lane "
              "stays =====", flush=True)
        v.report()
        return 1

    # ---- C2: pin-on-pur (4,15) -----------------------------------------
    print("  [C2] pin-on-pur (4,15)", flush=True)
    mo415, ncas415, c_pur415, nf415 = window_of(
        ADDR_415, mol_c1, mf_c1, spec, mf_d2h)
    pinned_section(v, "(4,15)", mol_d2h, mf_d2h, mo415, nf415,
                   ncas415, c_pur415, (4, 0), pt2_roots, sa=True)

    # ---- C3: pin-on-pur (6,13), the class re-derivations ---------------
    # the pricing = the ground root only (the class check is the
    # point; the full per-root ladder pricing is C2's business)
    c3_roots = [0]
    print("  [C3] pin-on-pur (6,13): the as-built distilled window",
          flush=True)
    mo613, ncas613, c_pur613, nf613 = window_of(
        ADDR_613, mol_c1, mf_c1, spec, mf_d2h)
    pinned_section(v, "(6,13)", mol_d2h, mf_d2h, mo613, nf613,
                   ncas613, c_pur613, (5, 1), c3_roots, sa=False)
    print("  [C3] pin-on-pur (6,13): the @r1 window (step B's)",
          flush=True)
    act_s, ncas_s, purity_loss = act_s_mo_of(mol_c1, mf_c1, mol_d2h,
                                             spec)
    print(f"  [C3] the per-irrep orthonormalization: worst dominant-"
          f"block weight loss {purity_loss:.2e}", flush=True)
    mo_s = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nf613),
                         act_s, c_core=c_pur613)
    pinned_section(v, "(6,13)@r1", mol_d2h, mf_d2h, mo_s, nf613,
                   ncas_s, c_pur613, (5, 1), c3_roots, sa=False)

    # ---- C4: the pinned SA on the lean lane window (4,15) --------------
    # THE CERTIFIED NO-WINDOW CLASS, DISCLOSED (the handoff's "labeling
    # may refuse on the internal re-projection"): the lane window's
    # promoted 14p columns are core-span content (the 4_15 leg's
    # record; ~5e-4 core-lean), so the SA's internal re-assembly
    # (core projection -> global Loewdin) cannot deliver columns pure
    # enough for pyscf's orb labeling -- measured 3.3e-06 raw, still
    # 2.7e-07 after the disclosed cleanup.  The pinned-(4,15) SA
    # result is the PIN-ON-PUR route's (C2 above).  The refusal is
    # caught and disclosed; the lane substrate chain itself is
    # witnessed (the C1-seeded retry fired).
    print("  [C4] the lean lane window's pinned SA (4,15)", flush=True)
    t0 = time.time()
    lmol, lmf, c_core_l, nf_l, mo_l, ncas_l, _act_l = lane_route(
        spec, ADDR_415_LANE, (4, 0))
    c_read_l = np.asarray(mo_l)[:, nf_l:nf_l + ncas_l]
    from vv.symmetry import symmetrize_columns as _sym
    _clean_l, imp_l = _sym(lmol, lmf.get_ovlp(), c_read_l, tol=1.0)
    print(f"  [C4] the assembled lane window's irrep impurities: max "
          f"{max(imp_l):.2e} (the labeling seam's dirt scale)",
          flush=True)
    try:
        out_sa_l, _ = pinned_sa_run(lmol, lmf, mo_l, nf_l, ncas_l,
                                    (4, 0), L, PARITY, c_core_l,
                                    c_read=_clean_l, mc_opts=BUDGET)
        blocks_l, _exp_l, _ = make_blocks(lmol, lmf,
                                          np.asarray(out_sa_l["mo"]),
                                          nf_l, ncas_l, (4, 0), L,
                                          PARITY, c_core_l,
                                          out_sa_l["roots"],
                                          c_read=_clean_l)
        l2s_l, s2s_l = meters_all(blocks_l, _exp_l)
        print(f"  [C4 sa] ensemble {out_sa_l['e_tot']:.8f}, converged "
              f"{out_sa_l['converged']} (cycles {out_sa_l['cycles']}), "
              f"<L2> [{min(l2s_l):.6f}, {max(l2s_l):.6f}] "
              f"({time.time() - t0:.1f} s)", flush=True)
        e2s_l = [row.e_corr for _, mc in blocks_l
                 for row in nevpt_rows(mc, pt2_roots)]
        at72 = sum(abs(x - 72.0) <= 0.01 for x in l2s_l)
        v.note("C4: the lean lane window's pinned SA (4,15)",
               f"{len(l2s_l)} components, {at72} at 72 +- 0.01; E(2) "
               f"mean {np.mean(e2s_l) * 1e3:+.4f} mEh")
        # the lane-core route is outside the table's keep-set: the
        # disclosure is the note above, never a row
    except ValueError as e:
        v.note("C4: the lean lane window's pinned SA is REFUSED (the "
               "certified no-window class at this partition)",
               f"pyscf's orb labeling rejects the SA's internally "
               f"re-projected columns: {e} -- the promoted 14p "
               f"columns are core-span content (the 4_15 leg's "
               f"certified mechanism); the pinned-(4,15) SA result "
               f"is the pin-on-pur route's (C2).  The lane chain "
               f"itself fired: plain D2h diverged, the C1-seeded "
               f"retry converged")
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
