"""The declaration layer: one spec per atom case, the campaign registry,
and the launch entry that executes a declaration.

A case is DECLARED, not guessed: a :class:`CaseSpec` fixes the atom, the
reference data, the valence window, and the sectors -- the intent
manifest. :func:`run_case` executes the declaration
(the term tables are the contract, always produced); nothing
auto-generates. ``CAMPAIGN`` is the
per-atom registry the cell layer and the runners derive from --
nothing else re-declares atoms, windows, or the gated references (the
one-declaration rule).

Division of labor with the cell layer: this module declares WHAT is
studied and runs the classic staged pipeline (ROHF -> CASCI per Ms
sector -> CASSCF -> SA-CASSCF, checkpointed under stage tags);
:mod:`atomterm.cells` addresses HOW each single calculation of the
methods grid runs -- one frozen :class:`~atomterm.cells.CellKey` per
axis combination, its grids derived from these specs. The campaign's
two run vocabularies (stage tags here, cell addresses there) are a
recorded phase-5 consideration, not an accident.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from atomterm.detbasis import sector_dimension
from atomterm.system import AtomConfig, ncore_of, sector_na_nb
from atomterm.reference import ReferenceSpec, run_rohf
from atomterm.scans import ScanSpec
from atomterm.windows import build_full_mo, occupied_core_idx, select_valence_ao
from atomterm.casscf import (
    run_casci,
    run_casscf as casscf_stage,
    run_sa_casscf,
)
from atomterm.stages import stage_facts
from atomterm.analysis import format_report

__all__ = ["CaseSpec", "CAMPAIGN", "run_case", "run_atom", "run_carbon",
           "main"]


def _ms_label(na: int, nb: int) -> str:
    """Filename/report label of a sector's M_S = (na - nb) / 2.

    Whole numbers stay bare ("1", "0" -- matching the historic carbon chkfile
    names); halves insert a "p" so paths stay clean ("0p5" for a doublet).
    """
    assert na >= nb, "sectors are built with na >= nb"
    two_ms = na - nb
    if two_ms % 2 == 0:
        return str(two_ms // 2)
    return f"{two_ms / 2:.1f}".replace(".", "p")


@dataclass(frozen=True)
class CaseSpec:
    """The declaration of one atom case -- what to run and what not.

    The term tables (T1/T2) are the contract: every executed case
    produces them. The single-root and SA stages are switched by their
    own flags.

    Attributes
    ----------
    symbol / basis / ecp / spin : the atom (:class:`~atomterm.system.AtomConfig`
        inputs; ``basis`` is a label like ``"ccpvdz"`` or a dict).
    nelec_act / ncore : the active-electron count and the closed-shell
        count (the fold).
    shell_patterns : the valence-shell label prefixes of the standard
        window (e.g. ``["2p"]``; ``["3d", "4s"]``).
    planned_windows : further window labels of the campaign, declared
        but not yet run (Ti: the (4,10)/(4,14) partitions) -- the cell
        grid derives their full factor grids (:func:`atomterm.cells.planned`).
    lattice_windows : the free-ion d2 lattice spans (empty = none) --
        CAS(2,·) cells over growing m-complete shell spans
        (:func:`atomterm.cells.lattice_keys`).
    scan : the declared frozen-scan anatomy
        (:class:`~atomterm.scans.ScanSpec` -- family pair, angle grid,
        fixed extras; None = no scan; the runner is
        ``scripts/radial_scan.py``).
    tag : the workdir/checkpoint label (defaults to the symbol).
    ref : the :class:`~atomterm.reference.ReferenceSpec` (character-forced
        occupations + basin gate) where default SCF filling fails; None
        lets plain aufbau run.
    sa_weights : state-average weights for the primary sector (uniform
        over the ground term's M_L set by convention).
    run_casscf : the single-root stage (off on multi-solution atoms).
    nroots_primary : PySCF roots for the primary CASCI.
    """

    symbol: str
    basis: str | dict
    nelec_act: int
    ncore: int
    shell_patterns: list[str]
    ecp: str | None = None
    spin: int = 0
    tag: str | None = None
    ref: ReferenceSpec | None = None
    sa_weights: list[float] | None = None
    run_casscf: bool = True
    nroots_primary: int = 9
    planned_windows: tuple[str, ...] = ()
    lattice_windows: tuple[str, ...] = ()
    scan: ScanSpec | None = None

    def atom_config(self) -> AtomConfig:
        """The spec as an :class:`AtomConfig` (``tag`` defaults to the
        symbol)."""
        return AtomConfig(symbol=self.symbol, basis=self.basis, ecp=self.ecp,
                          spin=self.spin, tag=self.tag or self.symbol)


# The declared campaign (one spec per atom) -- THE intent manifest: the
# cell grid (cells.window_grid/grid_keys/planned), the regen runner, and
# the table renderers all derive from these; nothing else re-declares
# atoms, windows, or the gated references (the one-declaration rule).
# New atom = new row here.
CAMPAIGN: dict[str, CaseSpec] = {
    "C": CaseSpec("C", "ccpvdz", nelec_act=2, ncore=2,
                  shell_patterns=["2p"], spin=2, sa_weights=[1 / 3] * 3,
                  scan=ScanSpec(fam=("2p", "3p"),
                                grid=(0.0, 0.8, 17))),
    "Si": CaseSpec("Si", "ccpvdz", nelec_act=2, ncore=6,
                   shell_patterns=["3p"], spin=2, sa_weights=[1 / 3] * 3,
                   scan=ScanSpec(fam=("3p", "4p"),
                                 grid=(0.0, 1.1, 23))),
    "Sc": CaseSpec("Sc", "ccpvdz", nelec_act=3, ncore=9,
                   shell_patterns=["3d", "4s"], spin=1,
                   sa_weights=[1 / 5] * 5,
                   scan=ScanSpec(fam=("3d", "4d"), grid=(0.0, 0.4, 9),
                                 # 4s fixed: the 4s->5s direction is
                                 # harmful and would mask the d structure
                                 extra_fixed="4s"),
                   ref=ReferenceSpec(
                       occ_closed=("1s", "2s", "2p", "3s", "3p", "4s"),
                       occ_open=("3d",), occ_e_ref=-759.73625510)),
    "Ti": CaseSpec("Ti", "ccpvdz", nelec_act=4, ncore=9,
                   shell_patterns=["3d", "4s"], spin=2,
                   sa_weights=[1 / 7] * 7,
                   ref=ReferenceSpec(
                       occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
                       occ_open=("3d",), occ_e_ref=-848.40648875),
                   planned_windows=("3d4d", "3d4s4p4d"),
                   lattice_windows=("3d", "3d4p", "3d4p4d")),
}


def run_case(spec: CaseSpec, workdir: str = "elements") -> dict:
    """Execute a declared case: the staged pipeline's launch entry for
    any atom, new or old (an atom is runnable iff declared here).

    Parameters
    ----------
    spec : the :class:`CaseSpec` to execute.
    workdir : where the checkpoints and reports land (generated inserts
        go to ``elements/<Atom>/``).

    Returns
    -------
    the :func:`run_atom` stage dict (``rohf_energy`` + one entry per
    computed stage tag).
    """
    res = run_atom(
        spec.atom_config(),
        list(spec.shell_patterns),
        spec.nelec_act,
        workdir=workdir,
        sa_weights=list(spec.sa_weights) if spec.sa_weights else None,
        run_casscf=spec.run_casscf,
        nroots_primary=spec.nroots_primary,
        ref=spec.ref,
    )
    return res


def run_atom(
    cfg: AtomConfig,
    shell_patterns: list[str],
    nelec_act: int,
    workdir: str = "elements",
    sa_weights: list[float] | None = None,
    run_casscf: bool = True,
    nroots_primary: int = 9,
    ref: ReferenceSpec | None = None,
) -> dict:
    """End-to-end pipeline for any single atom.

    Parameters
    ----------
    cfg : the :class:`AtomConfig` of the atom.
    shell_patterns : valence-shell label prefixes selecting the active AO
        space, concatenated in order (e.g. ``["2p"]`` for C 2p^2,
        ``["3d", "4s"]`` for Sc 3d^1 4s^2). Each pattern must match at least
        one AO; each should contribute a complete m-set (rotation closure).
    nelec_act : electrons treated as active. The Ms sectors follow from
        cfg.spin = 2S: the primary sector (na, nb) = ((n+2S)/2, (n-2S)/2),
        plus an Ms=0 sector (n/2, n/2) when n is even and S > 0 (skipped when
        odd -- no Ms=0 exists, e.g. a doublet).
    sa_weights : run SA-CASSCF over these state weights (primary sector).
    run_casscf : run the single-root CASSCF stage (off for transition metals,
        where orbital optimization is multi-solution).
    nroots_primary : PySCF roots for the primary CASCI (clamped to the sector
        size; the Ms=0 CASCI always diagonalizes the full sector so the
        complete multiplet table is produced).
    ref : passed through to :func:`run_rohf` as a :class:`ReferenceSpec`
        -- the character-forced occupation patterns + reference-basin
        energy gate, mandatory where the default ROHF is an aufbau
        runaway (Ti; see :func:`run_rohf`).
    workdir : where checkpoints and reports land.

    Raises
    ------
    ValueError
        ``spin`` negative or of unequal parity with ``nelec_act`` (no
        such Ms sector exists).

    Every stage is checkpointed under ``workdir`` and analyzed with the exact
    L/S operator matrices; see :func:`run_carbon` for the carbon instance.

    Returns
    -------
    dict
        ``"rohf_energy"`` plus one payload per stage tag
        (``casci_ms*``, ``casscf``, ``sa-casscf``), each carrying the
        analysis object, the run's energies, window weights, max|dE|,
        and the particle-number sanity.
    """
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label
    spin = cfg.spin  # 2S (PySCF convention)
    if spin < 0 or (nelec_act + spin) % 2 != 0:
        raise ValueError(
            f"spin 2S={spin} and nelec_act={nelec_act} must be >= 0 and "
            "of equal parity")

    # ---- stage 1: ROHF (checkpointed) -------------------------------------
    rohf_chk = str(out / f"{label}_rohf.chk")
    mol, mf = run_rohf(cfg, rohf_chk, ref=ref)
    print(f"[ROHF      ] E = {mf.e_tot:.8f} Eh   (chk: {rohf_chk})")

    act_idx = select_valence_ao(mol, shell_patterns)
    ncas = len(act_idx)

    # Ms sectors: (na, nb, nroots) for each CASCI diagonalization
    na_p, nb_p = sector_na_nb(nelec_act, spin)
    sector_defs = [(na_p, nb_p, min(nroots_primary, sector_dimension(ncas, na_p, nb_p)))]
    if nelec_act % 2 == 0 and spin > 0:
        na0 = nelec_act // 2
        sector_defs.append((na0, na0, sector_dimension(ncas, na0, na0)))

    results: dict[str, dict] = {"rohf_energy": float(mf.e_tot)}
    reports: list[str] = []

    def stage_analysis(tag: str, cas, c_act, na: int, nb: int):
        """The shared per-stage render: fold + analyze through
        :func:`stage_facts`, print the report block, accumulate the
        results dict (one closure per run -- it captures the run's
        mol/mf/reports)."""
        ncore = ncore_of(mol.nelectron, na + nb)
        # fold with the stage's OWN core (inside stage_facts: it uses
        # cas.mo_coeff[:, :ncore]): orbital optimization may legally
        # reorganize the inactive span, and a different fold than the one
        # PySCF's energies used breaks E_tot = E_core + e_cas. For CASCI
        # stages the passed core IS the mf core (built by build_full_mo).
        ms = (na - nb) / 2  # the true M_S (chkfile tags label M_S too)
        ms_txt = "0" if ms == 0 else f"{ms:+g}"
        sector = f"Ms={ms_txt} (na={na}, nb={nb})"
        facts = stage_facts(mol, mf, cas, c_act, ncore, na, nb, tag=sector)
        title = (f"[{tag}]  total energies vs PySCF: "
                 f"max|dE| = {facts.max_de_vs_pyscf:.2e} Eh")
        text = format_report(facts.analysis, title)
        off_shell = 1.0 - facts.window_rows.max(axis=1)
        text += (f"\nwindow composition (l-weights): "
                 + ", ".join(f"{n}: {v:.3f}" for n, v in
                             facts.window.items() if v > 1e-4)
                 + f"   (mean off-dominant-l weight {off_shell.mean():.2e};"
                   " l-purity is necessary, not sufficient, for closure)")
        print(text)
        reports.append(text)
        results[tag] = {
            "analysis": facts.analysis,
            "energies_pyscf": np.asarray(facts.energies_pyscf).tolist(),
            "window_weights": facts.window,
            "max_de_vs_pyscf": facts.max_de_vs_pyscf,
            # particle-number sanity of the built determinant space: N is
            # diagonal by construction and must equal na + nb everywhere
            "n_dev": facts.n_dev,
        }

    # ---- stage 2: CASCI in every Ms sector (checkpointed) ------------------
    mo_by_sector: dict[tuple[int, int], np.ndarray] = {}
    for na, nb, nroots in sector_defs:
        ncore = ncore_of(mol.nelectron, na + nb)
        mo = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
        mo_by_sector[(na, nb)] = mo
        tag = f"casci_ms{_ms_label(na, nb)}"
        cas = run_casci(mf, ncas, (na, nb), mo, str(out / f"{label}_{tag}.chk"),
                        nroots=nroots, verbose=0)
        print(f"[CASCI Ms={_ms_label(na, nb)}] roots: {np.round(cas.energies, 6)}")
        stage_analysis(tag, cas, cas.mo_coeff[:, ncore:ncore + ncas], na, nb)

    mo_p = mo_by_sector[(na_p, nb_p)]
    ncore_p = ncore_of(mol.nelectron, na_p + nb_p)

    # ---- stage 3: CASSCF (single root) on the primary sector ---------------
    if run_casscf:
        casscf_chk = str(out / f"{label}_casscf.chk")
        cass = casscf_stage(mf, ncas, (na_p, nb_p), mo_p, casscf_chk, verbose=0)
        print(f"[CASSCF    ] E = {cass.energies[0]:.8f} Eh")
        c_cass = cass.mo_coeff[:, ncore_p : ncore_p + ncas]
        stage_analysis("casscf", cass, c_cass, na_p, nb_p)

    # ---- stage 4: SA-CASSCF over the requested state weights ---------------
    if sa_weights:
        sa_chk = str(out / f"{label}_sa-casscf.chk")
        sa = run_sa_casscf(mf, ncas, (na_p, nb_p), mo_p, list(sa_weights),
                           sa_chk, verbose=0)
        print(f"[SA-CASSCF ] roots: {np.round(sa.energies, 6)}")
        c_sa = sa.mo_coeff[:, ncore_p : ncore_p + ncas]
        stage_analysis("sa-casscf", sa, c_sa, na_p, nb_p)

    # ---- persist -----------------------------------------------------------
    (out / f"{label}_report.md").write_text(
        f"# Atomic term analysis: {label}\n\n"
        f"ROHF/cc-pVDZ energy: {mf.e_tot:.8f} Eh\n\n" + "\n\n".join(reports),
        encoding="utf-8",
    )
    serial = {
        tag: {
            "energies_pyscf": payload["energies_pyscf"],
            "states": [vars(s) for s in payload["analysis"].states],
            "lz_spectrum": payload["analysis"].lz_spectrum,
            "sz_spectrum": payload["analysis"].sz_spectrum,
            "comm_norm_h_lz": payload["analysis"].comm_norm_h_lz,
            "comm_norm_h_l2": payload["analysis"].comm_norm_h_l2,
            "window_weights": payload["window_weights"],
            "max_de_vs_pyscf": payload["max_de_vs_pyscf"],
        }
        for tag, payload in results.items()
        if tag != "rohf_energy"
    }
    (out / f"{label}_report.json").write_text(
        json.dumps({"rohf_energy": results["rohf_energy"], **serial}, indent=2),
        encoding="utf-8",
    )
    print(f"\nreport written to {out}/{label}_report.md and .json")
    return results


def run_carbon(cfg: AtomConfig | None = None, workdir: str = "elements") -> dict:
    """Carbon pilot (2p^2, CAS(2,3)): the generic pipeline with carbon data."""
    return run_atom(
        cfg or AtomConfig(),
        ["2p"],
        2,
        workdir=workdir,
        sa_weights=[1 / 3] * 3,  # SA over the three 3P M_L components (Ms=1)
    )


def main():
    """Run the carbon pilot (the ``python -m`` entry)."""
    run_carbon()


if __name__ == "__main__":
    main()
