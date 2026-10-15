"""SA-deviation decomposition + AVAS comparison + Sc tilt-follow curves
+ purified-core before/after (ticket 16 Phase 2, ticket 13).

Renders, per atom (C / Si / Sc):
1. the SA-family ladder rows with the four-way deviation attribution:
   - many-body commutators (what the term analysis sees),
   - one-particle span defects |[P,Lz]| / |[P,Lp]| (the window itself),
   - the same defects on the rung's INPUT window (the fact-17 baseline),
   - span overlap vs the raw shell span (the identity change gauge),
   - the core-fold witness: refold the rung's converged core with the RAW
     shell window (attributes the deviation to core/fold vs window).
2. the purified-core verdict (ticket 13): the input baselines before
   (canonical core) and after (purified core) -- machine precision -- and
   the *_pur rung rows through the same attribution.
3. the AVAS window rows (independent ecosystem selector).
4. Sc only: frozen tilt-follow curves with <L^2> tracking -- the
   deterministic anatomy of the member-asymmetric descent (fact 18d).

Usage:
    python scripts/sa_decomposition.py [C|Si|Sc|all] [workdir]
"""

import sys
from pathlib import Path

import numpy as np
from pyscf import lib

from atomterm.detbasis import commutator, frobenius_norm
from atomterm.drivers import (
    avas_window,
    build_full_mo,
    purified_core_orbitals,
    radial_basis,
    run_rohf,
    select_valence_ao,
    span_overlap,
)
from atomterm.ladder import run_ladder
from atomterm.operators import build_operators, window_closure
from atomterm.scans import scan_state, scan_window
from atomterm.stages import stage_facts
from atomterm.system import AtomConfig

SPECS = {
    "C": dict(cfg=AtomConfig(), patterns=["2p"], nelec_act=2,
              sa_weights=[1 / 3] * 3, frac=(2, 3), avas_labels=["C 2p"]),
    "Si": dict(cfg=AtomConfig("Si", "ccpvdz", spin=2), patterns=["3p"],
               nelec_act=2, sa_weights=[1 / 3] * 3, frac=(2, 3),
               avas_labels=["Si 3p"]),
    "Sc": dict(cfg=AtomConfig("Sc", "ccpvdz", spin=1), patterns=["3d", "4s"],
               nelec_act=3, sa_weights=[1 / 5] * 5, frac=(3, 6),
               avas_labels=["Sc 3d", "Sc 4s"]),
}
RUNGS = ["casci_avg", "casscf_relax", "casscf_radial",
         "sa", "sa_radial", "sa_pyscf"]
PUR_RUNGS = ["casci_avg_pur", "casscf_relax_pur", "casscf_radial_pur",
             "sa_relax_pur", "sa_pur", "sa_radial_pur"]


def _refold_comm(mol, mf, c_shell, ncore, na, nb, c_core):
    """Core-fold witness: many-body commutators of the rung's converged
    core folded against the RAW shell window (0 => the deviation is in
    the window, not the core/fold)."""
    ops = build_operators(mol, mf, c_shell, ncore, na, nb, c_core=c_core)
    return (float(frobenius_norm(commutator(ops["H"], ops["Lz"]))),
            float(frobenius_norm(commutator(ops["H"], ops["L2"]))))


def decompose(sym: str, workdir: Path) -> list[str]:
    spec = SPECS[sym]
    cfg = spec["cfg"]
    out = workdir
    rungs = run_ladder(cfg, spec["patterns"], spec["nelec_act"],
                       workdir=str(out), sa_weights=spec["sa_weights"],
                       frac_nelec_open=spec["frac"][0],
                       frac_n_open=spec["frac"][1], purified=True,
                       rung_tags=RUNGS + PUR_RUNGS)

    label = cfg.label
    mol, mf = run_rohf(cfg, str(out / f"{label}_rohf.chk"))
    act_idx = select_valence_ao(mol, spec["patterns"])
    ncore = (mol.nelectron - spec["nelec_act"]) // 2
    spin = cfg.spin
    na = (spec["nelec_act"] + spin) // 2
    nb = (spec["nelec_act"] - spin) // 2
    ncas = len(act_idx)
    eye = np.eye(mol.nao_nr())
    c_shell = eye[:, act_idx]
    s = mf.get_ovlp()

    # input windows of the SA family (fact 17 baselines) + the purified
    # substrate before/after (ticket 13: the baseline IS the core span)
    mo_shell = build_full_mo(mf, list(range(ncore)), act_idx)
    mo_rad, _ = radial_basis(mf, list(range(ncore)), act_idx)
    c_pur, contam = purified_core_orbitals(mol, mf, ncore)
    mo_shell_p = build_full_mo(mf, list(range(ncore)), act_idx, c_core=c_pur)
    mo_rad_p, _ = radial_basis(mf, list(range(ncore)), act_idx, c_core=c_pur)
    inputs = {
        "shell-start": window_closure(mol, mo_shell[:, ncore:ncore + ncas]),
        "radial-start": window_closure(mol, mo_rad[:, ncore:ncore + ncas]),
        "canonical-start": window_closure(
            mol, mf.mo_coeff[:, ncore:ncore + ncas]),
    }
    pur_inputs = {
        "pur shell-start": window_closure(
            mol, mo_shell_p[:, ncore:ncore + ncas]),
        "pur radial-start": window_closure(
            mol, mo_rad_p[:, ncore:ncore + ncas]),
    }

    lines = [f"# SA-deviation decomposition: {sym}", "",
             "input-window baselines (one-particle closure defects):"]
    for name, c in inputs.items():
        lines.append(f"  {name:<16} |[P,Lz]| = {c['lz']:9.2e}   "
                     f"|[P,Lp]| = {c['lperp']:9.2e}")
    lines += ["  purified core substrate (ticket 13; the fact-17 fix):"]
    for name, c in pur_inputs.items():
        lines.append(f"  {name:<16} |[P,Lz]| = {c['lz']:9.2e}   "
                     f"|[P,Lp]| = {c['lperp']:9.2e}")
    lines += [f"  core-span closure: "
              f"|[P,Lz]| {window_closure(mol, c_pur)['lz']:.1e}  "
              f"|[P,Lp]| {window_closure(mol, c_pur)['lperp']:.1e}  "
              f"(max foreign-l contamination {contam.max():.1e})", "",
              f"{'rung':<20} {'|[H,Lz]|':>9} {'|[H,L2]|':>9} "
              f"{'|[P,Lz]|':>9} {'|[P,Lp]|':>9} {'sv_min':>7} "
              f"{'refold_lz':>9}  attribution"]
    shell_purity = max(rungs["casci_avg"]["window"].values())
    for tag in RUNGS + PUR_RUNGS:
        if tag not in rungs:
            continue
        r = rungs[tag]
        chk = str(out / f"{label}_{tag}.chk")
        mo = lib.chkfile.load(chk, "mcscf/mo_coeff")
        if mo is None:
            # CASCI stages do not persist mo_coeff: the window actually
            # folded is the rung's INPUT window (raw shell AOs for the
            # canonical rungs, the projected pur window for *_pur); the
            # core block is the stage's Loewdin-rotated one (same span)
            if tag.endswith("_pur"):
                mo = mo_shell_p[:, :ncore + ncas]
            else:
                mo = np.hstack([mo_shell[:, :ncore], c_shell])
        w = mo[:, ncore:ncore + ncas]
        clo = window_closure(mol, w)
        sv = span_overlap(c_shell, w, s)
        rf_lz, _ = _refold_comm(mol, mf, c_shell, ncore, na, nb,
                                mo[:, :ncore])
        purity = max(r["window"].values())
        notes = []
        base_lz = inputs["shell-start"]["lz"]
        base_lp = inputs["shell-start"]["lperp"]
        if clo["lz"] < 1e-10 and clo["lperp"] < 1e-10:
            notes.append("shell span (closed)")
        elif (clo["lperp"] <= 3 * max(base_lp, 1e-15)
              and clo["lz"] <= 3 * max(base_lz, 1e-15)):
            if base_lp > 1e-10 and (clo["lperp"] < base_lp / 10
                                    or clo["lz"] < base_lz / 10):
                notes.append("healed below input baseline")
            else:
                notes.append("at input baseline")
        else:
            if sv < 0.9:
                notes.append("identity change")
            elif purity < shell_purity - 0.05:
                notes.append("l-mixing")
            if clo["lperp"] > 0.1:
                notes.append("member-tilted")
            elif sv < 0.999:
                notes.append("radial relaxation")
        if rf_lz > 1e-8:
            notes.append(f"core/fold {rf_lz:.1e}")
        lines.append(
            f"{tag:<20} {r['comm_norm_h_lz']:>9.2e} {r['comm_norm_h_l2']:>9.2e} "
            f"{clo['lz']:>9.2e} {clo['lperp']:>9.2e} {sv:>7.4f} "
            f"{rf_lz:>9.1e}  " + "; ".join(notes))

    # ---- AVAS row
    from atomterm.drivers import run_casci

    mo_avas, ncas_a, _ = avas_window(mf, spec["avas_labels"], ncore)
    if ncas_a:
        w = mo_avas[:, ncore:ncore + ncas_a]
        cas = run_casci(mf, ncas_a, (na, nb), mo_avas,
                        str(out / f"{label}_avas.chk"), nroots=1, name="avas")
        fa = stage_facts(mol, mf, cas, w, ncore, na, nb, tag="avas")
        clo = window_closure(mol, w)
        sv = span_overlap(c_shell, w, s)
        d_e = fa.e_ground - rungs["casci_avg"]["e_ground"]
        lines += ["",
                  f"AVAS window ({', '.join(spec['avas_labels'])}):",
                  f"  span overlap {sv:.4f}  |[P,Lz]| {clo['lz']:.1e}  "
                  f"|[P,Lp]| {clo['lperp']:.1e}",
                  f"  CASCI E {fa.e_ground:.8f} ({d_e * 1e3:+.1f} mEh vs "
                  f"casci_avg)  |[H,Lz]| {fa.analysis.comm_norm_h_lz:.1e}  "
                  f"|[H,L2]| {fa.analysis.comm_norm_h_l2:.1e}  "
                  f"<L^2> {fa.ground_l2:.4f} ({fa.ground_term})"]
    return lines


def sc_tilt_curves(workdir: Path) -> list[str]:
    """Sc deep-dive: follow the member-asymmetric descent deterministically,
    tracking <L^2> (the term-identity meter)."""
    cfg = SPECS["Sc"]["cfg"]
    out = workdir
    mol, mf = run_rohf(cfg, str(out / f"{cfg.label}_rohf.chk"))
    from atomterm.drivers import select_shell_ao

    d3, d4 = select_shell_ao(mol, "3d"), select_shell_ao(mol, "4d")
    c4s = np.eye(mol.nao_nr())[:, select_shell_ao(mol, "4s")].reshape(-1, 1)
    ncore, na, nb = 9, 2, 1

    def state(thetas, with_l2=True):
        c_act = np.hstack([scan_window(mol, mf, thetas, d3, d4), c4s])
        return scan_state(mol, mf, c_act, ncore, na, nb, with_l2=with_l2)

    e0 = state([0.0] * 5).e_active
    lines = ["", "# Sc tilt-follow curves (frozen, deterministic)", "",
             "delta  subset          dE(uEh)   |[H,Lz]|  |[H,L2]|  <L^2>  L"]
    for subset in [(0, 1, 2), (0, 1, 3, 4)]:
        for delta in np.linspace(0.05, 0.35, 7):
            p = state([delta if i in subset else 0.0 for i in range(5)])
            ell = (-1 + np.sqrt(1 + 4 * p.ground_l2)) / 2
            lines.append(
                f"{delta:5.2f}  {str(subset):<14} "
                f"{(p.e_active - e0) * 1e6:+9.1f}  {p.comm_h_lz:8.2e} "
                f"{p.comm_h_l2:8.2e}  {p.ground_l2:6.3f}  {ell:4.2f}")
    # m-symmetric reference: the real radial direction (rises, fact 18d)
    lines += ["", "m-symmetric reference (all members together):"]
    for delta in (0.05, 0.15, 0.30):
        p = state([delta] * 5)
        ell = (-1 + np.sqrt(1 + 4 * p.ground_l2)) / 2
        lines.append(f"  delta {delta:.2f}: dE {(p.e_active - e0) * 1e3:+8.2f} mEh  "
                     f"<L^2> {p.ground_l2:.3f}  (L {ell:.1f})")
    return lines


def main():
    sym = sys.argv[1] if len(sys.argv) > 1 else "all"
    workdir = Path(sys.argv[2] if len(sys.argv) > 2 else "results")
    workdir.mkdir(parents=True, exist_ok=True)
    syms = list(SPECS) if sym == "all" else [sym]
    for s in syms:
        lines = decompose(s, workdir)
        if s == "Sc":
            lines += sc_tilt_curves(workdir)
        text = "\n".join(lines)
        (workdir / f"{SPECS[s]['cfg'].label}_sa_decomposition.md"
         ).write_text(text + "\n", encoding="utf-8")
        print(text + "\n")


if __name__ == "__main__":
    main()
