"""The committed instrument behind the 2026-11-23 V&V verdict
corrections: the measurements the map and the leg disclosures cite,
reproducible from the committed code.

Probe A -- the (6,13) window triangulation (the corrected slice, the
key's own ncore): the assembled span vs the projected span
(``act - P*act``, S-orthonormalized) vs the distilled act -- principal
angles and the closure meters ``||[P,Lz]||``/``||[P,Lperp]||`` on all
three plus ``c_pur``.  Controls: the assembled-vs-projected cosines
all 1 (the assembly is span-preserving), the closure meters ~1e-13
(rotation closure -- the frozen ground's definite-L certificate), the
act-vs-assembled tilt = the 7s projection survival.

Probe B -- the (6,13) frozen ground's true reading (the correct meter
slice): ``<L2>``/``<S2>`` of the CASCI ground and the per-orbital
occupations.  Controls: <L2> = <S2> = 12.000000 (the 7F: L = 3, S = 3,
the d5 high-spin class + one f), the distilled s column unoccupied.

Probe C -- the (4,15) pinned-route mechanism: the distilled p columns'
overlaps with the occupied 6p MOs (the band = the occupied shell's own
radial); the pinned core's per-block content ranks (the C1 reference's
m-mixed occupied 6p spans the promoted slots); the m-pure feed's
``purify_window`` refusal (the vanish is feed-gauge-independent).

Run from the repo root:
    OMP_NUM_THREADS=8 python scripts/probes/probe_u14e_vv_span.py \\
        | tee logs/probe_u14e_vv_span.log
"""

import sys

import numpy as np
from pyscf import mcscf
from pyscf.fci import direct_spin1

sys.path.insert(0, "scripts/u14e")
from _vv_common import (_irrep_projectors, dominant_irreps,  # noqa: E402
                        lane_route, load_reference, meters, purify_window,
                        symmetrize_columns)
from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.cells import CellKey, key_core, window_columns  # noqa: E402
from atomterm.windows import (build_full_mo, lowdin,  # noqa: E402
                              occupied_core_idx)
from atomterm.operators import window_closure  # noqa: E402


def probe_a():
    mol, mf = load_reference()
    spec = campaign_spec("U14e")
    key = CellKey.from_address("u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg")
    nc, _c_core, c_pur = key_core(key, mol, mf, spec.shell_bands)
    assert nc == 4, f"the (6,13) key's core pairs: {nc} != 4"
    act, ncas, _rep = window_columns(key, mol, mf, spec.shell_bands,
                                     spec.shell_radials, ncol=nc)
    S = mf.get_ovlp()
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, nc), act,
                            c_core=c_pur)
    B = mo_prod[:, nc:nc + ncas]
    P = c_pur @ np.linalg.solve(c_pur.T @ S @ c_pur, c_pur.T @ S)
    Vp = lowdin(act - P @ act, S)
    labels = [l.split()[-1] for l in mol.ao_labels()]
    print("[A] act's own column order, dominant AO per column:")
    print("   ", [labels[int(np.argmax(np.abs(act[:, j])))]
                  for j in range(ncas)])

    def pcos(X, Y):
        return np.linalg.svd(X.T @ S @ Y, compute_uv=False)

    print("[A] principal cosines act vs assembled:",
          np.array2string(pcos(act, B), precision=6, max_line_width=200))
    print("[A] principal cosines projected vs assembled:",
          np.array2string(pcos(Vp, B), precision=6, max_line_width=200))
    for name, C in (("act", act), ("projected", Vp), ("assembled", B),
                    ("c_pur", c_pur)):
        cl = window_closure(mol, C)
        print(f"[A] closure[{name:>9s}]: lz={cl['lz']:.2e} "
              f"lperp={cl['lperp']:.2e}")


def probe_b():
    mol, mf = load_reference()
    spec = campaign_spec("U14e")
    key = CellKey.from_address("u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg")
    nc, _c_core, c_pur = key_core(key, mol, mf, spec.shell_bands)
    act, ncas, _rep = window_columns(key, mol, mf, spec.shell_bands,
                                     spec.shell_radials, ncol=nc)
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, nc), act,
                            c_core=c_pur)
    cas = mcscf.CASCI(mf, ncas, (5, 1))
    cas.canonicalization = False
    cas.kernel(mo_prod)
    ok = meters(mol, mo_prod[:, nc:nc + ncas], 5, 1)
    bad = meters(mol, mo_prod[:, 5:5 + ncas], 5, 1)
    ci = np.asarray(cas.ci).ravel()          # nroots=1: the vector itself
    l_ok = ok("l2", ci)
    l_bad = bad("l2", ci)
    l_def = min(abs(l_ok - lx * (lx + 1)) for lx in range(16))
    print(f"[B] frozen CASCI e_tot {cas.e_tot:.8f}")
    print(f"[B] correct slice: <L2> {l_ok:.6f} (closest L(L+1) off "
          f"{l_def:.2e}), <S2> {ok('s2', ci):.6f}")
    print(f"[B] shifted slice (+1 column): <L2> {l_bad:.6f} -- the "
          f"incident's artifact")
    dm1 = np.asarray(direct_spin1.make_rdm1(ci, ncas, (5, 1)))
    occ = (np.real(dm1.diagonal()) if dm1.ndim == 2
           else np.real(dm1[0].diagonal() + dm1[1].diagonal()))
    print("[B] per-orbital occupations (act order; last = the "
          "distilled s): "
          + np.array2string(occ, precision=3, max_line_width=200))


def probe_c():
    spec = campaign_spec("U14e")
    addr = "u14e_10f13d14p@6s6p7s_4-0_sa-casscf_ao_-_pur_avg"
    lmol, lmf, c_core, nf, mo, ncas, act = lane_route(spec, addr, (4, 0))
    S = lmf.get_ovlp()
    labels = [l.split()[-1] for l in lmol.ao_labels()]
    print("[C] act's dominant AO per column:",
          [labels[int(np.argmax(np.abs(act[:, j])))] for j in range(ncas)])
    occ = occupied_core_idx(lmf, nf)
    print(f"[C] the distilled p columns vs the occupied closed MOs "
          f"{occ} (|<v,S,MO>|):")
    for j in range(ncas - 3, ncas):
        v = act[:, j]
        o = [abs(float(v @ S @ lmf.mo_coeff[:, i])) for i in occ]
        print(f"    col {j}: " + " ".join(f"{x:.3f}" for x in o))

    v, _imp = symmetrize_columns(lmol, S, act)
    ids = dominant_irreps(lmol, S, v)
    projs = _irrep_projectors(lmol, S)
    print("[C] the pinned core's per-block content ranks (the C1 "
          "m-mixed occupied 6p, S-projected onto the declared span):")
    for k in sorted(set(ids)):
        cols = [j for j, kid in enumerate(ids) if kid == k]
        d, Ginv_dTS = projs[k]
        Bk = d @ (Ginv_dTS @ c_core)
        gb, Ub = np.linalg.eigh(Bk.T @ S @ Bk)
        keep = gb > 1e-10 * gb.max()
        Bn = Bk @ Ub[:, keep] / np.sqrt(gb[keep])
        sub = v[:, cols] - Bn @ (Bn.T @ S @ v[:, cols])
        post = np.sqrt(np.diag(sub.T @ S @ sub))
        print(f"    {lmol.irrep_name[k]:>4s}: cols {cols}, core rank "
              f"{int(keep.sum())}, post-subtraction norms "
              f"{np.array2string(post, precision=3)}")

    try:
        purify_window(lmol, S, c_core, act)
        print("[C] m-pure feed: PASSES -- the vanish NOT reproduced?")
    except ValueError as e:   # the vanish or the tolerance miss --
        # the refusal's form is draw-dependent, the verdict constant
        print(f"[C] m-pure feed: refused (feed-gauge-independent): {e}")


if __name__ == "__main__":
    probe_a()
    probe_b()
    probe_c()
