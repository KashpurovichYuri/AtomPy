"""The CAS(4,16)/CAS(4,18)/CAS(4,21) V&V leg -- the correlation-channel
ladder at the certified (4,0) partition (ncore = 5: the 6s2 6p6 7s2
pairs folded; the 4 open electrons = 5f3 6d1).

THE WINDOWS (each in the incumbent- and the @r1-source form where the
sources carry them; the s/p correlation channels at band width):
  (4,16)    5f(7) + 6d(5) + 14p(3) + 7s(1) -- the s channel: the
            folded pair's own correlation direction (the pair MO
            ITSELF is core content here -- the window takes the
            distilled virtual radial); the incumbent form via its
            key, the @r1 form via the tag route;
  (4,18)    5f(7) + 6d(5) + the 6p shell at k=2 -- the incumbent's
            certified 6p depth (the tail certificate REFUSED k=3:
            the occupied-span Gram's third 6p eigenvalue is tail --
            the reference carries 2 significant 6p radials);
  (4,21)@r1 the 6p shell at k=3 from the DENSITY blocks (the 6p
            trio's density is rank-3 -- the composition the
            incumbent cannot host).  The composition is not
            expressible as a key address (the tag's k is uniform per
            window), so these windows are built by DIRECT
            distillation calls with per-shell k -- window_distilled
            and series_window already support per-shell k; src is
            untouched.

THE RUNS per window: frozen CASCI / CASSCF / SA (free@pur) + the
pinned SA (pinned@pur, the wrapper route); PT2 on the ground roots;
the SA meters over all 17 roots (the E(2) cell = root 0's, the table
convention).  The frozen classes DISCLOSED, not gated; the pinned SAs
on the incumbent-style windows may blow up PT2's Sr class (the
diagnosed pinned-gauge pathology) -- the values stay visible as the
witness.

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_4_16_21.py \
        | tee logs/u14e_vv_4_16_21.log
"""

import sys
import time

import numpy as np
from pyscf import mcscf

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))

from _vv_common import (Verdict, conv_of, core_part_of, drop_row,  # noqa: E402
                        fmt_range, free_kernel, load_reference,
                        make_blocks, meters, meters_all, nevpt_rows,
                        pinned_sa_run, price, std_row)

SECTOR = (4, 0)
NCORE = 5
ROOTS = 17
L, PARITY = 8, "u"
BUDGET = {"max_cycle_macro": 200}
SHELL_BANDS = {"5f": ("9f", "10f", "11f"), "6d": ("12d", "13d", "14d"),
               "6p": ("13p", "14p", "15p")}
K_18_R1 = {"5f": 1, "6d": 1, "6p": 2}
K_21_R1 = {"5f": 1, "6d": 1, "6p": 3}


def _print_rep(rep):
    for sh, r in rep.items():
        caps = []
        for c in r["capture"]:
            caps.extend(c) if isinstance(c, (list, tuple)) \
                else caps.append(c)
        caps = [c for c in caps if c]
        cap_s = ("capture n/a" if not caps else
                 f"capture {min(caps):.4f}-{max(caps):.4f}")
        disc = r.get("discarded")
        disc_s = (f", discarded {disc:.4f}, "
                  if disc is not None else ", ")
        print(f"  [window] {sh}: k={r['k']} of {len(r['bands'])} "
              f"bands{disc_s}{cap_s}", flush=True)


def _win_416(key_addr, mol, mf, spec):
    """(4,16): the key-driven windows (the incumbent via the plain
    address, the @r1 via the tag route) -- all four shells are
    declared in the case's shell_radials."""
    from atomterm.cells import CellKey, key_core, window_columns
    key = CellKey.from_address(key_addr)
    nc, _cc, c_pur = key_core(key, mol, mf, spec.shell_bands)
    act, ncas, rep = window_columns(key, mol, mf, spec.shell_bands,
                                    spec.shell_radials, ncol=nc)
    _print_rep(rep)
    return act, ncas, c_pur


def _win_k(mol, mf, spec, source, radials):
    """The 5f/6d/6p windows with per-shell k, by source: the
    incumbent (window_distilled on the occupied-span Gram) or the
    @r1 density blocks (series_window).  The tag grammar's uniform k
    cannot carry these compositions; the signatures support them."""
    from atomterm.windows import window_distilled
    from atomterm.basis import series_window
    from atomterm.core import purified_core_orbitals
    s = mf.get_ovlp()
    c_occ = mf.mo_coeff[:, np.asarray(mf.mo_occ) > 0]
    if source == "inc":
        act, rep = window_distilled(mol, c_occ, s, list(radials),
                                    SHELL_BANDS, radials)
    else:
        act, rep = series_window(mol, [mf], s, ["5f", "6d", "6p"],
                                 SHELL_BANDS, radials, NCORE)
    _print_rep(rep)
    c_pur = purified_core_orbitals(mol, mf, NCORE)[0]
    return act, act.shape[1], c_pur


def main() -> int:
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.operators import window_closure
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    import u14e_vv_unification as uni
    mol_d2h, mf_d2h = uni.wrapper_substrate(mol, mf, spec)
    weights = [1.0 / ROOTS] * ROOTS

    CASES = (
        ("(4,16)", lambda: _win_416(
            "u14e_10f13d14p7s_4-0_casci_ao_-_pur_avg", mol, mf, spec)),
        ("(4,16)@r1", lambda: _win_416(
            "u14e_5f6d14p7s@r1_4-0_casci_ao_-_pur_avg", mol, mf, spec)),
        ("(4,18)@r1", lambda: _win_k(mol, mf, spec, "r1",
                                     K_18_R1)),
        ("(4,21)@r1", lambda: _win_k(mol, mf, spec, "r1",
                                     K_21_R1)),
    )

    for tag, win in CASES:
        print(f"===== {tag} =====", flush=True)
        act, ncas, c_pur = win()
        mo_prod = build_full_mo(mf, occupied_core_idx(mf, NCORE), act,
                                c_core=c_pur)
        mo_win = mo_prod[:, NCORE:NCORE + ncas]
        cl = window_closure(mol, mo_win)
        print(f"  [closure] ||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
              f"{cl['lperp']:.2e}", flush=True)
        expect = meters(mol, mo_win, *SECTOR)
        v.check(f"{tag}: the window is rotation-closed (the "
                "definite-L certificate)",
                cl["lz"] <= 1e-10 and cl["lperp"] <= 1e-10,
                f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
                f"{cl['lperp']:.2e}")

        # ---- frozen CASCI (the pricing class) --------------------------
        t0 = time.time()
        cas_f = mcscf.CASCI(mf, ncas, SECTOR)
        cas_f.canonicalization = False
        cas_f.kernel(mo_prod)
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
        drop_row(std_row("casci", "free@pur", tag, float(cas_f.e_tot),
                         r_f["e2"] * 1e3,
                         float(cas_f.e_tot) + r_f["e2"],
                         f"{l_f:.4f}", f"{r_f['s2']:.4f}",
                         "definite-L (the class is the disclosure)",
                         conv=conv_of(cas_f.fcisolver.conv_tol, None)))

        # ---- free CASSCF (the escape draw) ------------------------------
        t0 = time.time()
        mc_c, cas_c = free_kernel(mf, ncas, SECTOR, mo_prod,
                                  track_iters=True)
        r_c = price(cas_c, expect, roots=[0])[0]
        print(f"  [casscf] e_tot {mc_c.e_tot:.8f}, converged "
              f"{mc_c.converged}, <L2> {r_c['l2']:.6f} "
              f"({time.time() - t0:.1f} s)", flush=True)
        drop_row(std_row("casscf", "free@pur", tag, float(mc_c.e_tot),
                         r_c["e2"] * 1e3,
                         float(mc_c.e_tot) + r_c["e2"],
                         f"{r_c['l2']:.4f}", f"{r_c['s2']:.4f}",
                         "escape draw",
                         conv=conv_of(mc_c.final_de, mc_c.n_macro,
                                      bool(mc_c.converged))))

        # ---- free SA (the term average) ---------------------------------
        t0 = time.time()
        mc_sa, cas_sa = free_kernel(mf, ncas, SECTOR, mo_prod,
                                    weights=weights, nroots=ROOTS,
                                    track_iters=True)
        core_e, e_cas = core_part_of(cas_sa)
        expect_sa = meters(mol, cas_sa.mo_coeff[:, NCORE:NCORE + ncas],
                           *SECTOR)
        l2s_sa = [expect_sa("l2", np.asarray(cas_sa.ci[k]).ravel())
                  for k in range(ROOTS)]
        s2s_sa = [expect_sa("s2", np.asarray(cas_sa.ci[k]).ravel())
                  for k in range(ROOTS)]
        r0 = price(cas_sa, expect_sa, roots=[0])[0]
        print(f"  [sa] ensemble {mc_sa.e_tot:.8f}, converged "
              f"{mc_sa.converged}, <L2> "
              f"[{min(l2s_sa):.6f}, {max(l2s_sa):.6f}] "
              f"({time.time() - t0:.1f} s)", flush=True)
        drop_row(std_row("sa-casscf", "free@pur", tag,
                         float(e_cas[0] + core_e), r0["e2"] * 1e3,
                         float(e_cas[0] + core_e + r0["e2"]),
                         fmt_range(l2s_sa), f"{max(s2s_sa):.4f}",
                         "the term average",
                         conv=conv_of(mc_sa.final_de, mc_sa.n_macro,
                                      bool(mc_sa.converged))))

        # ---- pinned SA (pinned@pur) -------------------------------------
        mo_d2h = build_full_mo(mf_d2h,
                               occupied_core_idx(mf_d2h, NCORE), act,
                               c_core=c_pur)
        nf = NCORE
        c_read = mo_d2h[:, nf:nf + ncas]
        t0 = time.time()
        out, _ = pinned_sa_run(mol_d2h, mf_d2h, mo_d2h, nf, ncas,
                               SECTOR, L, PARITY, c_pur, c_read=c_read,
                               mc_opts=BUDGET)
        blocks, expect_p, _ = make_blocks(mol_d2h, mf_d2h,
                                          np.asarray(out["mo"]), nf,
                                          ncas, SECTOR, L, PARITY,
                                          c_pur, out["roots"],
                                          c_read=c_read)
        l2s_p, s2s_p = meters_all(blocks, expect_p)
        n72 = sum(abs(x - 72.0) <= 0.01 for x in l2s_p)
        n90 = sum(abs(x - 90.0) <= 0.01 for x in l2s_p)
        print(f"  [pinned sa] ensemble {out['e_tot']:.8f}, converged "
              f"{out['converged']} (cycles {out['cycles']}, "
              f"{time.time() - t0:.1f} s); <L2> [{min(l2s_p):.6f}, "
              f"{max(l2s_p):.6f}]; the counts: 72: {n72}, 90: {n90}, "
              f"other: {len(l2s_p) - n72 - n90}", flush=True)
        e2s_p = [row.e_corr for _, mc in blocks
                 for row in nevpt_rows(mc, [0])]
        drop_row(std_row("sa-casscf", "pinned@pur", tag,
                         float(out["e_tot"]),
                         float(np.mean(e2s_p) * 1e3),
                         float(out["e_tot"] + np.mean(e2s_p)),
                         fmt_range(l2s_p), f"{max(s2s_p):.4f}",
                         "draw-dependent partition",
                         conv=conv_of(out["conv_dE"], out["cycles"],
                                      bool(out["converged"]))))

    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
