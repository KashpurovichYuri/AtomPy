"""The carbon ladder (ticket 08): every orbital-treatment rung through one
analysis pass, compared all-vs-all.

Rungs, in order of increasing orbital freedom:

1. casci_raw    -- CASCI on the raw shell window, ROHF core folded as-is
                   (no spherical averaging): the deliberately symmetry-broken
                   reference. Expect the 3P multiplet split (~4e-5 Eh,
                   {M_L=+-1}/{0} pattern from the quadrupolar core field) and
                   commutators ~ eps ~ 1e-5.
2. casci_avg    -- same but with the spherically averaged core density:
                   exact zeros. The scalar projection of the ROHF core.
3. casci_frac   -- CASCI with the average-of-configuration core
                   (fractional-occupation SCF), folded WITHOUT averaging:
                   the Fock was exactly spherical, so invariance holds by
                   construction; energies must match rung 2 to ~1e-10 --
                   the O(eps^2) test of "averaging = ensemble projection".
4. casscf_relax -- CASSCF with the window span frozen (only core<->virtual
                   rotations): adds core polarization, no escape possible.
5. casscf_shell -- free CASSCF from the raw-shell start: symmetry breaking
                   is physical output (<L^2> ~ 1.67).
6. casscf_pyscf -- free CASSCF from PySCF's default start: start-dependence.
7. sa3          -- SA-CASSCF over the three 3P M_L components: restoration.

Every rung is analyzed with the SAME operator machinery on its OWN window
and core; the comparison matrix below reports energies, ground <L^2>, both
commutator norms, max|dE| vs PySCF, and the window composition.

Sector subtlety: in the maximal-Ms carbon sector (3 determinants) the space
is a COMPLETE L = 1 irrep, so L^2 = 2*Id identically and [H, L^2] = 0 for
ANY H -- it carries no information there. The raw-core fingerprints in this
sector are [H, Lz] (~1e-5) and the multiplet splitting (~4e-5 Eh, with the
{M_L=0} root shifted off the {M_L=+-1} pair). In sectors mixing several
terms (Ms = 0) both commutators are informative.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from atomterm.system import AtomConfig
from atomterm.drivers import (
    build_full_mo,
    run_casci,
    run_casscf,
    run_casscf_core_relax,
    run_casscf_pyscf_start,
    run_frac_occ_scf,
    run_rohf,
    run_sa_casscf,
    select_valence_ao,
    shell_weights,
)
from atomterm.operators import build_operators, core_anisotropy
from atomterm.analysis import analyze

__all__ = ["run_carbon_ladder"]


def _rung(tag: str, kind: str, cas, mol, mf, c_act: np.ndarray, ncore: int,
          na: int, nb: int, spherical_core: bool) -> dict:
    """One analysis pass over a finished stage (shared by all rungs)."""
    ops = build_operators(mol, mf, c_act, ncore, na, nb,
                          spherical_core=spherical_core,
                          c_core=cas.mo_coeff[:, :ncore])
    res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag)
    e_mine = np.linalg.eigvalsh(ops["H"]) + cas.core_energy
    e_pyscf = np.sort(np.atleast_1d(np.asarray(cas.energies)))[: len(e_mine)]
    w_names, w = shell_weights(mol, c_act, mf.get_ovlp())
    ground = res.states[0]
    return {
        "tag": tag,
        "kind": kind,
        "n_det": res.n_det,
        "e_ground": float(e_mine[0]),
        "energies_pyscf": np.asarray(cas.energies).tolist(),
        "max_de_vs_pyscf": float(np.max(np.abs(e_mine[: len(e_pyscf)] - e_pyscf))),
        "ground_term": ground.term,
        "ground_l2": ground.l2,
        "terms": dict(Counter(s.term for s in res.states)),
        "comm_norm_h_lz": res.comm_norm_h_lz,
        "comm_norm_h_l2": res.comm_norm_h_l2,
        "window": {n: round(float(v), 6) for n, v in zip(w_names, w.mean(axis=0))},
    }


def run_carbon_ladder(cfg: AtomConfig | None = None, workdir: str = "results") -> dict:
    """Run all ladder rungs for carbon 2p^2 CAS(2,3), print and persist the
    all-vs-all comparison. Returns ``{rung_tag: payload}``."""
    cfg = cfg or AtomConfig()
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label

    mol, mf = run_rohf(cfg, str(out / f"{label}_rohf.chk"))
    mol_f, mf_f = run_frac_occ_scf(cfg, str(out / f"{label}_fracocc.chk"),
                                   nelec_open=2, n_open=3)

    act_idx = select_valence_ao(mol, ["2p"])
    ncas = len(act_idx)
    na, nb = 2, 0
    ncore = (mol.nelectron - (na + nb)) // 2
    c_shell = np.eye(mol.nao_nr())[:, act_idx]
    chk = lambda name: str(out / f"{label}_{name}.chk")  # noqa: E731

    rungs: dict[str, dict] = {}

    # -- CASCI rungs: identical window, three core treatments ----------------
    casci_specs = [
        ("casci_raw", mf, None, False),
        ("casci_avg", mf, None, True),
        ("casci_frac", mf_f, mf_f.mo_coeff[:, :ncore], False),
    ]
    for tag, mf_ref, c_core, sph in casci_specs:
        mo = build_full_mo(mf_ref, list(range(ncore)), act_idx, c_core=c_core)
        cas = run_casci(mf_ref, ncas, (na, nb), mo, chk(tag), nroots=3)
        rungs[tag] = _rung(tag, "casci", cas, mol, mf_ref, c_shell, ncore,
                           na, nb, spherical_core=sph)
        rungs[tag]["core_anisotropy"] = round(core_anisotropy(
            mol, mf_ref, ncore, c_core=c_core), 12)

    # -- CASSCF rungs ---------------------------------------------------------
    mo_shell = build_full_mo(mf, list(range(ncore)), act_idx)
    casscf_specs = [
        ("casscf_relax", lambda: run_casscf_core_relax(
            mf, ncas, (na, nb), mo_shell, chk("casscf_relax"))),
        ("casscf_shell", lambda: run_casscf(
            mf, ncas, (na, nb), mo_shell, chk("casscf_shell"))),
        ("casscf_pyscf", lambda: run_casscf_pyscf_start(
            mf, ncas, (na, nb), chk("casscf_pyscf"))),
        ("sa3", lambda: run_sa_casscf(
            mf, ncas, (na, nb), mo_shell, [1 / 3] * 3, chk("sa3"))),
    ]
    for tag, run in casscf_specs:
        cas = run()
        c_act = cas.mo_coeff[:, ncore : ncore + ncas]
        rungs[tag] = _rung(tag, "casscf", cas, mol, mf, c_act, ncore,
                           na, nb, spherical_core=True)

    _report(rungs, out, label, mf)
    return rungs


def _report(rungs: dict, out: Path, label: str, mf) -> None:
    """Print and persist the all-vs-all comparison matrix."""
    order = ["casci_raw", "casci_avg", "casci_frac",
             "casscf_relax", "casscf_shell", "casscf_pyscf", "sa3"]
    rungs = {k: rungs[k] for k in order}

    lines = [f"# Carbon ladder: {label}", "",
             f"ROHF reference energy: {mf.e_tot:.8f} Eh", "",
             f"{'rung':<14} {'E_ground':>16} {'<L^2>':>9} {'|[H,Lz]|':>10} "
             f"{'|[H,L2]|':>10} {'max|dE|':>9} {'off-shell':>9}  ground terms"]
    for tag, r in rungs.items():
        off_shell = 1.0 - max(r["window"].values())
        lines.append(
            f"{tag:<14} {r['e_ground']:>16.8f} {r['ground_l2']:>9.5f} "
            f"{r['comm_norm_h_lz']:>10.2e} {r['comm_norm_h_l2']:>10.2e} "
            f"{r['max_de_vs_pyscf']:>9.1e} {off_shell:>9.1e}  {r['terms']}"
        )

    # pairwise ground-state energy differences (all-vs-all)
    lines += ["", "pairwise ground-energy differences (E_i - E_j), Eh:"]
    tags = list(rungs)
    lines.append("            " + " ".join(f"{t[:11]:>12}" for t in tags))
    for ti in tags:
        row = [f"{rungs[ti]['e_ground'] - rungs[tj]['e_ground']:>12.3e}"
               for tj in tags]
        lines.append(f"{ti:<12}" + " ".join(row))

    text = "\n".join(lines)
    print(text)
    (out / f"{label}_ladder.md").write_text(text + "\n", encoding="utf-8")
    (out / f"{label}_ladder.json").write_text(
        json.dumps(rungs, indent=2), encoding="utf-8")
