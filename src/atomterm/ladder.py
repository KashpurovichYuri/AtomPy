"""The orbital-treatment ladder: every rung through one analysis pass,
compared all-vs-all (ticket 08 on carbon, ticket 09 on Si).

Rungs, in order of increasing orbital freedom:

1. casci_raw    -- CASCI on the raw shell window, reference core folded as-is
                   (no spherical averaging): the deliberately symmetry-broken
                   reference. Fingerprint: [H, Lz] ~ eps and the multiplet
                   split {M_L=0} vs {+-1} (quadrupole k=2, m=0 core field;
                   the dipole k=1 is excluded by parity). In a sector that
                   is a complete single irrep L^2 = L(L+1)*Id identically,
                   so [H, L^2] carries no information there.
2. casci_avg    -- spherically averaged core: exact zeros. The scalar
                   projection of the polarized reference core.
3. casci_frac   -- average-of-configuration core (fractional-occupation SCF),
                   folded WITHOUT averaging: spherical at the source.
4. casscf_relax -- window span frozen; only core<->virtual rotations: adds
                   the core's polarization by the correlated active density.
5. casscf_virtmix -- core<->window exchange blocked; window may borrow
                   virtual character: the active pair's radial correlation,
                   isolated from identity-changing rotations.
6. casscf_shell -- free CASSCF from the raw-shell start.
7. casscf_pyscf -- free CASSCF from PySCF's default start (start dependence).
8. sa / sa_pyscf -- SA-CASSCF over the given weights, from the raw-shell and
                   the canonical start (restoration + start dependence).

Divergence guard: every rung reports ``ortho_err`` = max|C^dag S C - I|.
Energies from rungs with ortho_err > 1e-6 are meaningless (the orbital
optimizer lost orthonormality -- observed on Si from raw-shell starts,
where reported energies plunge BELOW the variational bound). Such rungs
are flagged DIVERGED in the report and their numbers must not be read as
physics.

max|dE| dictionary (our reconstruction vs PySCF's own energies):
  0.0          -- identical folds (raw/frac rungs), or a spherical core
                  (SA: the ensemble density is a rotational scalar);
  ~1e-5..5e-5  -- the energy footprint of a POLARIZED core that we average
                  and PySCF uses raw (avg, optimized rungs). A polarization
                  detector, not an error.

Sector subtlety: in a maximal-Ms sector that is a complete single irrep
(carbon/Si 3-det Ms=1), [H, L^2] = 0 for ANY H -- only [H, Lz] and the
splitting are fingerprints there. In sectors mixing several terms both
commutators are informative.
"""

from __future__ import annotations

import json
import math
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
    run_casscf_virtual_mix,
    run_frac_occ_scf,
    run_rohf,
    run_sa_casscf,
    select_valence_ao,
    shell_weights,
)
from atomterm.operators import build_operators, core_anisotropy
from atomterm.analysis import analyze, configuration_summary

__all__ = ["run_carbon_ladder", "run_ladder"]


def _rung(tag: str, kind: str, cas, mol, mf, c_act: np.ndarray, ncore: int,
          na: int, nb: int, spherical_core: bool, s: np.ndarray) -> dict:
    """One analysis pass over a finished stage (shared by all rungs)."""
    ops = build_operators(mol, mf, c_act, ncore, na, nb,
                          spherical_core=spherical_core,
                          c_core=cas.mo_coeff[:, :ncore])
    res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag)
    e_mine = np.linalg.eigvalsh(ops["H"]) + cas.core_energy
    e_pyscf = np.sort(np.atleast_1d(np.asarray(cas.energies)))[: len(e_mine)]
    w_names, w = shell_weights(mol, c_act, s)
    ground = res.states[0]
    mo = cas.mo_coeff
    ortho_err = float(np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max())
    return {
        "tag": tag,
        "kind": kind,
        "n_det": res.n_det,
        "diverged": ortho_err > 1e-6,
        "ortho_err": ortho_err,
        "e_ground": float(e_mine[0]),
        "energies_pyscf": np.asarray(cas.energies).tolist(),
        "max_de_vs_pyscf": float(np.max(np.abs(e_mine[: len(e_pyscf)] - e_pyscf))),
        "ground_term": ground.term,
        "ground_l2": ground.l2,
        "terms": dict(Counter(s_.term for s_ in res.states)),
        "comm_norm_h_lz": res.comm_norm_h_lz,
        "comm_norm_h_l2": res.comm_norm_h_l2,
        "window": {n: round(float(v), 6) for n, v in zip(w_names, w.mean(axis=0))},
        "configurations": configuration_summary(
            mol, s, cas.ci[: len(e_mine)], c_act, na, nb),
    }


def run_ladder(cfg: AtomConfig, shell_patterns: list[str], nelec_act: int,
               workdir: str = "results", sa_weights: list[float] | None = None,
               frac_nelec_open: int | None = None, frac_n_open: int | None = None) -> dict:
    """Run all ladder rungs for one atom's primary Ms sector, print and
    persist the all-vs-all comparison. Returns ``{rung_tag: payload}``."""
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label

    mol, mf = run_rohf(cfg, str(out / f"{label}_rohf.chk"))
    s = mf.get_ovlp()
    mol_f = mf_f = None
    if frac_nelec_open is not None:
        mol_f, mf_f = run_frac_occ_scf(
            cfg, str(out / f"{label}_fracocc.chk"), frac_nelec_open, frac_n_open)

    act_idx = select_valence_ao(mol, shell_patterns)
    ncas = len(act_idx)
    spin = cfg.spin
    na, nb = (nelec_act + spin) // 2, (nelec_act - spin) // 2
    ncore = (mol.nelectron - nelec_act) // 2
    c_shell = np.eye(mol.nao_nr())[:, act_idx]
    chk = lambda name: str(out / f"{label}_{name}.chk")  # noqa: E731

    rungs: dict[str, dict] = {}
    ndet = (math.comb(ncas, na) * math.comb(ncas, nb))

    # -- CASCI rungs: identical window, core treatments ----------------------
    casci_specs = [
        ("casci_raw", mf, None, False),
        ("casci_avg", mf, None, True),
    ]
    if mf_f is not None:
        casci_specs.append(("casci_frac", mf_f, mf_f.mo_coeff[:, :ncore], False))
    for tag, mf_ref, c_core, sph in casci_specs:
        mo = build_full_mo(mf_ref, list(range(ncore)), act_idx, c_core=c_core)
        cas = run_casci(mf_ref, ncas, (na, nb), mo, chk(tag),
                        nroots=min(9, ndet), name=tag)
        rungs[tag] = _rung(tag, "casci", cas, mol, mf_ref, c_shell, ncore,
                           na, nb, spherical_core=sph, s=s)
        rungs[tag]["core_anisotropy"] = round(core_anisotropy(
            mol, mf_ref, ncore, c_core=c_core), 12)

    # -- CASSCF rungs ---------------------------------------------------------
    mo_shell = build_full_mo(mf, list(range(ncore)), act_idx)
    casscf_specs = [
        ("casscf_relax", lambda: run_casscf_core_relax(
            mf, ncas, (na, nb), mo_shell, chk("casscf_relax"))),
        ("casscf_virtmix", lambda: run_casscf_virtual_mix(
            mf, ncas, (na, nb), mo_shell, chk("casscf_virtmix"))),
        ("casscf_shell", lambda: run_casscf(
            mf, ncas, (na, nb), mo_shell, chk("casscf_shell"))),
        ("casscf_pyscf", lambda: run_casscf_pyscf_start(
            mf, ncas, (na, nb), chk("casscf_pyscf"))),
    ]
    if sa_weights:
        casscf_specs.append(("sa", lambda: run_sa_casscf(
            mf, ncas, (na, nb), mo_shell, list(sa_weights), chk("sa"))))
        # canonical-start variant: isolates the optimizer's start dependence
        # (Si's raw-shell start diverges; the canonical one converges)
        casscf_specs.append(("sa_pyscf", lambda: run_sa_casscf(
            mf, ncas, (na, nb), mf.mo_coeff, list(sa_weights), chk("sa_pyscf"))))
    for tag, run in casscf_specs:
        cas = run()
        c_act = cas.mo_coeff[:, ncore : ncore + ncas]
        rungs[tag] = _rung(tag, "casscf", cas, mol, mf, c_act, ncore,
                           na, nb, spherical_core=True, s=s)

    _report(rungs, out, label, mf)
    return rungs


def run_carbon_ladder(cfg: AtomConfig | None = None, workdir: str = "results") -> dict:
    """The carbon ladder (ticket 08): 2p^2 through every rung."""
    return run_ladder(cfg or AtomConfig(), ["2p"], 2, workdir=workdir,
                      sa_weights=[1 / 3] * 3, frac_nelec_open=2, frac_n_open=3)


def _report(rungs: dict, out: Path, label: str, mf) -> None:
    """Print and persist the all-vs-all comparison matrix."""
    order = ["casci_raw", "casci_avg", "casci_frac", "casscf_relax",
             "casscf_virtmix", "casscf_shell", "casscf_pyscf", "sa"]
    tags = [t for t in order if t in rungs]
    rungs = {t: rungs[t] for t in tags}
    ref = rungs.get("casci_avg")

    lines = [f"# Orbital-treatment ladder: {label}", "",
             f"ROHF reference energy: {mf.e_tot:.8f} Eh", "",
             f"{'rung':<14} {'E_ground':>16} {'dE_vs_avg':>10} {'<L^2>':>9} "
             f"{'|[H,Lz]|':>10} {'|[H,L2]|':>10} {'max|dE|':>9} {'ortho':>9}  ground config"]
    for tag in tags:
        r = rungs[tag]
        d_e = r["e_ground"] - ref["e_ground"] if ref else 0.0
        flag = " DIVERGED" if r["diverged"] else ""
        lines.append(
            f"{tag:<14} {r['e_ground']:>16.8f} {d_e:>10.3e} {r['ground_l2']:>9.5f} "
            f"{r['comm_norm_h_lz']:>10.2e} {r['comm_norm_h_l2']:>10.2e} "
            f"{r['max_de_vs_pyscf']:>9.1e} {r['ortho_err']:>9.1e}  "
            f"{r['configurations'][0]}{flag}"
        )

    lines += ["", "pairwise ground-energy differences (E_i - E_j), Eh:",
              "            " + " ".join(f"{t[:11]:>12}" for t in tags)]
    for ti in tags:
        row = [f"{rungs[ti]['e_ground'] - rungs[tj]['e_ground']:>12.3e}"
               for tj in tags]
        lines.append(f"{ti:<12}" + " ".join(row))

    lines += ["", "configurations per root:"]
    for tag in tags:
        for i, cfg_str in enumerate(rungs[tag]["configurations"]):
            lines.append(f"  {tag:<14} root {i}: {cfg_str}")

    text = "\n".join(lines)
    print(text)
    (out / f"{label}_ladder.md").write_text(text + "\n", encoding="utf-8")
    (out / f"{label}_ladder.json").write_text(
        json.dumps(rungs, indent=2), encoding="utf-8")
