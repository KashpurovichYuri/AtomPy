"""The @r1 CAS-ladder V&V leg -- the prior-axis window at (4,12) and
(4,15), free routes, so the table's @r1 story spans the flagship
CAS-size axis at matched source (the user's comparability directive;
the (6,13)@r1 rows live in the companion leg ``u14e_vv_6_13_r1.py``,
whose s channel carries the pair-MO projection).

THE WINDOWS (the production recipe, written out per section):
  (4,12) ``5f6d@r1``     the density-block radials for the chemical
                         5f/6d shells (the declared 3-band
                         realizations; k=1 each) -- the campaign's
                         certified frozen record: 5L 72.000000 at
                         +149.3 mEh;
  (4,15) ``5f6d14p@r1``  the same f/d radials + the PROMOTED 14p
                         band at its 1:1 identity (undeclared shell
                         = the band itself -- the promotion's
                         semantics preserved; only the f/d radial
                         SIGNAL changes vs the incumbent).
  Both on the pur core (ncore=5: the 7s2 pair folded), the trusted
  rung as the @r head.  Free CASCI / CASSCF / SA (the campaign
  budgets: max_cycle_macro 200) + PT2 per ground root -> the
  ``...(4,12)@r1`` / ``...(4,15)@r1`` table rows.

THE CHECKS: the rotation-closure certificate per window (the
definite-L guarantee); the frozen ground's definite L (WHICH L is
the disclosure -- the campaign anchors (4,12)@r1 at the 5L; the
(4,15) p-promotion classes may ground elsewhere frozen, per the
unification leg's C2 finding).

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_r1_ladder.py \
        [--ground-only]  |  tee logs/u14e_vv_r1_ladder.log
"""

import sys
import time

import numpy as np
from pyscf import mcscf

from vv.meters import meters
from vv.reference import load_reference
from vv.routes import core_part_of, free_kernel, price
from vv.table import conv_of, drop_row, fmt_range, std_row
from vv.verdict import Verdict

ADDR_412_R1 = "u14e_5f6d@r1_4-0_casci_ao_-_pur_avg"
ADDR_415_R1 = "u14e_5f6d14p@r1_4-0_casci_ao_-_pur_avg"
L_L2 = 72.0
ROOTS = 17


def production_window(key_addr, ncore, sector, mol, mf, spec):
    """One key's window, written out (``key_core`` + ``window_columns``
    -- the @r1 tag routes inside to ``basis.series_window`` with the
    rung as the head) + ``build_full_mo``; the closure meters printed
    (the definite-L certificate)."""
    from atomterm.cells import CellKey, key_core, window_columns
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.operators import window_closure

    key = CellKey.from_address(key_addr)
    nc, _c_core, c_pur = key_core(key, mol, mf, spec.shell_bands)
    assert nc == ncore, (f"the key's core pairs {nc} != the leg's "
                         f"NCORE {ncore}")
    act, ncas, rep = window_columns(key, mol, mf, spec.shell_bands,
                                    spec.shell_radials, ncol=nc)
    if rep is not None:
        for sh, r in rep.items():
            caps = [c for c in r["capture"] if c]
            cap_s = ("capture n/a" if not caps else
                     f"capture {min(min(c) for c in caps):.4f}"
                     f"-{max(max(c) for c in caps):.4f}")
            print(f"  [window] {sh}: k={r['k']} of {len(r['bands'])} "
                  f"bands, discarded {r['discarded']:.4f}, {cap_s}",
                  flush=True)
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, ncore), act,
                            c_core=c_pur)
    cl = window_closure(mol, mo_prod[:, ncore:ncore + ncas])
    print(f"  [closure] ||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
          f"{cl['lperp']:.2e}", flush=True)
    return mo_prod, ncas, cl


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pt2_roots = [0] if ground_only else None
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    weights = [1.0 / 17] * 17

    for tag, addr, ncore, sector in (
            ("(4,12)@r1", ADDR_412_R1, 5, (4, 0)),
            ("(4,15)@r1", ADDR_415_R1, 5, (4, 0))):
        print(f"===== {tag} =====", flush=True)
        mo_full, ncas, cl = production_window(addr, ncore, sector,
                                              mol, mf, spec)
        expect = meters(mol, mo_full[:, ncore:ncore + ncas], *sector)
        v.check(f"{tag}: the window is rotation-closed (the "
                "definite-L certificate)",
                cl["lz"] <= 1e-10 and cl["lperp"] <= 1e-10,
                f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
                f"{cl['lperp']:.2e}")

        # ---- free CASCI (the frozen pricing class) --------------------
        t0 = time.time()
        cas_f = mcscf.CASCI(mf, ncas, sector)
        cas_f.canonicalization = False
        cas_f.kernel(mo_full)
        r_f = price(cas_f, expect, roots=[0])[0]
        l_f = r_f["l2"]
        l_def = min(abs(l_f - lx * (lx + 1)) for lx in range(16))
        de_meh = (float(cas_f.e_tot) - mf.e_tot) * 1e3
        print(f"  [casci] e_tot {cas_f.e_tot:.8f}, <L2> {l_f:.6f}, "
              f"dE {de_meh:+.1f} mEh ({time.time() - t0:.1f} s)",
              flush=True)
        v.check(f"{tag}: the frozen ground reads a definite L",
                l_def <= 1e-4,
                f"<L2> {l_f:.6f} -- closest L(L+1) within "
                f"{l_def:.2e} (WHICH L is the disclosure)")
        v.note(f"{tag} casci", f"dE {de_meh:+.1f} mEh, E(2) "
               f"{r_f['e2'] * 1e3:+.4f} mEh, <S2> {r_f['s2']:.4f}")
        drop_row(std_row("casci", "free@pur", tag, float(cas_f.e_tot),
                         r_f["e2"] * 1e3,
                         float(cas_f.e_tot) + r_f["e2"],
                         f"{l_f:.4f}", f"{r_f['s2']:.4f}",
                         "definite-L (the class is the disclosure)",
                         conv=conv_of(cas_f.fcisolver.conv_tol, None)))

        # ---- free CASSCF (single root; the escape draw, no gate) ------
        t0 = time.time()
        mc_c, cas_c = free_kernel(mf, ncas, sector, mo_full,
                                  track_iters=True)
        r_c = price(cas_c, expect, roots=[0])[0]
        print(f"  [casscf] e_tot {mc_c.e_tot:.8f}, converged "
              f"{mc_c.converged}, <L2> {r_c['l2']:.6f} "
              f"({time.time() - t0:.1f} s)", flush=True)
        v.note(f"{tag} casscf", f"<L2> {r_c['l2']:.6f}, E(2) "
               f"{r_c['e2'] * 1e3:+.4f} mEh")
        drop_row(std_row("casscf", "free@pur", tag, float(mc_c.e_tot),
                         r_c["e2"] * 1e3,
                         float(mc_c.e_tot) + r_c["e2"],
                         f"{r_c['l2']:.4f}", f"{r_c['s2']:.4f}",
                         "escape draw",
                         conv=conv_of(mc_c.final_de, mc_c.n_macro,
                                      bool(mc_c.converged))))

        # ---- free SA (the term average; PT2 on every root) ------------
        t0 = time.time()
        mc_sa, cas_sa = free_kernel(mf, ncas, sector, mo_full,
                                    weights=weights, nroots=ROOTS,
                                    track_iters=True)
        core_e, e_cas = core_part_of(cas_sa)
        expect_sa = meters(mol, cas_sa.mo_coeff[:, ncore:ncore + ncas],
                           *sector)
        rows_sa = price(cas_sa, expect_sa,
                        roots=range(1) if ground_only else range(ROOTS))
        for r in rows_sa:
            r["e_ref"] = float(e_cas[r["root"]] + core_e)
        l2s_sa = [r["l2"] for r in rows_sa]
        r0 = rows_sa[0]
        dev_m = max(abs(x - L_L2) for x in l2s_sa)
        print(f"  [sa] ensemble {mc_sa.e_tot:.8f}, converged "
              f"{mc_sa.converged}, <L2> "
              f"[{min(l2s_sa):.6f}, {max(l2s_sa):.6f}] "
              f"({time.time() - t0:.1f} s)", flush=True)
        v.note(f"{tag} sa", f"<L2> max dev from 72 {dev_m:.2e}; E(2) "
               f"root 0 {r0['e2'] * 1e3:+.4f} mEh, total "
               f"{r0['e_ref'] + r0['e2']:.4f}")
        drop_row(std_row("sa-casscf", "free@pur", tag,
                         float(e_cas[0] + core_e), r0["e2"] * 1e3,
                         float(r0["e_ref"] + r0["e2"]),
                         fmt_range(l2s_sa),
                         f"{max(r['s2'] for r in rows_sa):.4f}",
                         "5L" if dev_m <= 1e-4 else "draw-dependent",
                         conv=conv_of(mc_sa.final_de, mc_sa.n_macro,
                                      bool(mc_sa.converged))))

    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
