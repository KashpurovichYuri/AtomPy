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

Division of labor
-----------------
With the cell layer: this module declares WHAT is studied and runs
the classic staged pipeline (ROHF -> CASCI per Ms
sector -> CASSCF -> SA-CASSCF, checkpointed under stage tags);
:mod:`atomterm.cells` addresses HOW each single calculation of the
methods grid runs -- one frozen :class:`~atomterm.cells.CellKey` per
axis combination, its grids derived from these specs. The campaign's
two run vocabularies (stage tags here, cell addresses there) are
DELIBERATE: the term tables are a stable contract, the grid an
evolving experiment -- coupling them would put the committed tables'
reproducibility at risk.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pyscf import gto

from atomterm.detbasis import sector_dimension
from atomterm.system import AtomConfig, ncore_of, sector_na_nb
from atomterm.reference import ReferenceSpec, run_rohf
from atomterm.scans import ScanSpec
from atomterm.windows import (build_full_mo, closed_electrons, fold_shells,
                              occupied_core_idx, select_valence_ao)
from atomterm.casscf import (
    run_casci,
    run_casscf as casscf_stage,
    run_sa_casscf,
)
from atomterm.stages import stage_facts
from atomterm.analysis import format_report

__all__ = ["CaseSpec", "Scenario", "CAMPAIGN", "scenarios_for",
           "run_case", "run_atom", "run_carbon"]


# the declaration ----------------------------------------------------

@dataclass(frozen=True)
class Scenario:
    """One declared scenario: a named, dated physics claim that a set
    of cells need not run.

    A scenario is DECLARED DATA, interpreted by exactly one consumer
    -- the run gate (:func:`atomterm.cells.void_reason`: a scenario'd
    cell is void-with-reason at its key, still addressable; a full
    regen materializes it).  The vocabulary is generic on purpose:
    what a scenario MEANS physically lives in its ``reason``.  (The
    render half the original spec declared never came -- presentation
    lives with the emitters.)

    Attributes
    ----------
    name : the scenario's provenance key (campaign defaults and
        per-atom declarations merge on it).
    cells : cell ADDRESS PREFIXES (:attr:`atomterm.cells.CellKey.address`)
        -- a cell is scenario'd iff its address starts with one.
    reason : the dated physics claim the gate carries ("... (2026-11-02)");
        ``--full`` re-tests it.
    enabled : the whole scenario's gate.
    """

    name: str
    cells: tuple[str, ...]
    reason: str
    enabled: bool = True

    def __post_init__(self):
        if not self.name:
            raise ValueError("a scenario needs a name")
        if not self.cells:
            raise ValueError(f"scenario {self.name!r} names no cells")
        if not self.reason:
            raise ValueError(f"scenario {self.name!r} carries no reason")

    def gates(self, address: str) -> bool:
        """Whether ``address`` is this scenario's cell (prefix match)."""
        return self.enabled and address.startswith(self.cells)


# campaign-level scenarios applied to EVERY atom unless its row
# overrides them by name (none today -- the merge is the tested seam;
# a campaign-wide scenario arrives with its first campaign-wide claim)
DEFAULT_SCENARIOS: dict[str, Scenario] = {}


def _ms_label(na: int, nb: int) -> str:
    """Filename/report label of a sector's M_S = (na - nb) / 2.

    Whole numbers stay bare ("1", "0" -- matching the committed chkfile
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
    lattice_windows : the lattice spans (empty = none) -- CASCI cells
        over growing m-complete shell spans at ``lattice_sector``
        (:func:`atomterm.cells.lattice_keys`).
    lattice_sector : the lattice series' sector ``(na, nb)`` -- WHAT
        the lattice studies: the free-ion series (the ns2 pair folded
        into the core) or the full valence.  Required when
        ``lattice_windows`` is non-empty.
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
    term : the certified ground term ``(L, parity)`` the pinned lane
        targets (parity = the configuration's total inversion parity);
        None = not certified for lane cells.
    lane_windows : the pinned-lane legs -- window labels carrying
        their own fold (``"3d3p4d@semi"``, :func:`atomterm.windows.fold_shells`
        the tag grammar); each becomes one CI-level lane cell
        (:func:`atomterm.cells.lane_keys`).  The tag marks the lane
        cell: a tagged window NEVER collides with a free-route label.
    scenarios : this atom's :class:`Scenario` declarations (merged
        over :data:`DEFAULT_SCENARIOS` by :func:`scenarios_for`).
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
    lattice_sector: tuple[int, int] | None = None
    scan: ScanSpec | None = None
    term: tuple[int, str] | None = None
    lane_windows: tuple[str, ...] = ()
    scenarios: tuple[Scenario, ...] = ()

    def __post_init__(self):
        """Validate the declaration's own grammar: the term target and
        every lane window (its tag resolves, its fold stays out of its
        window).  Fail here, at the seam that owns the grammar."""
        if self.lane_windows and self.term is None:
            raise ValueError(
                "lane windows need the term target (term=(L, 'g' | 'u'))")
        if self.lattice_windows and self.lattice_sector is None:
            raise ValueError(
                "lattice windows need the lattice sector "
                "(lattice_sector=(na, nb) -- the series' own Ms)")
        na, nb = self.lattice_sector or (1, 0)
        if not (na >= nb >= 0 and na + nb > 0):
            raise ValueError(
                f"lattice_sector {self.lattice_sector!r} must be "
                "(na >= nb >= 0, n > 0)")
        if self.term is not None:
            L, parity = self.term
            if not (isinstance(L, int) and L >= 0) or parity not in ("g", "u"):
                raise ValueError(
                    f"term {self.term!r} must be (L >= 0, 'g' | 'u')")
        open_d = next((p for p in self.shell_patterns
                       if p.endswith("d")), None)
        for label in self.lane_windows:
            span, sep, tag = label.partition("@")
            shells = re.findall(r"\d+[spdfgh]", span)
            if not sep or not tag or "".join(shells) != span:
                raise ValueError(
                    f"lane window {label!r} must be <shells>@<fold> "
                    "(e.g. '3d3p@semi')")
            fold = fold_shells(tag, open_d=open_d)
            overlap = sorted(set(shells) & set(fold))
            if overlap:
                raise ValueError(
                    f"lane window {label!r}: {overlap} sit in both the "
                    "window and its fold")

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

# The ECP atom's raw specs live in data/ (raw basis/ECP text carries no
# path-safe name -- hence the explicit tag). Its window tiers are the
# record's exponent-band lattices expressed as their exact shell
# subsets: this basis is its own exponent grid, so a band IS a set of
# whole shells -- m-complete like any shell window.
_DATA = Path(__file__).resolve().parents[2] / "data"
_U_BASIS = gto.parse((_DATA / "u_dyall_cut_basis.txt").read_text())
_U_ECP = (_DATA / "u_stuttgart_rlc_ecp78.txt").read_text()
_U_TIER = {
    "3f": ("9f", "10f", "11f", "13d"),
    "4f": ("8f", "9f", "10f", "11f", "13d"),
    "2f2d": ("9f", "10f", "12d", "13d"),
    "3f2d": ("9f", "10f", "11f", "12d", "13d"),
    "4f2d": ("8f", "9f", "10f", "11f", "12d", "13d"),
}

# the standing protocol for planned windows, enforced: their OPTIMIZER
# cells stay void until their CASCI anchors have been read (a full
# regen materializes on demand).  The prefixes name every label form
# the grid gives the windows (plain and NO-prefix).
_ANCHOR_FIRST = Scenario(
    name="anchor-first",
    cells=("ti_3d4d_3-1_casscf_", "ti_3d4dno2_3-1_casscf_",
           "ti_3d4d_3-1_sa-casscf_", "ti_3d4dno2_3-1_sa-casscf_",
           "ti_3d4s4p4d_3-1_casscf_", "ti_3d4s4p4dno4_3-1_casscf_",
           "ti_3d4s4p4d_3-1_sa-casscf_", "ti_3d4s4p4dno4_3-1_sa-casscf_"),
    reason="planned windows anchor first: no optimizer cell of a window "
           "whose CASCI anchors are unread (2026-11-02)")

CAMPAIGN: dict[str, CaseSpec] = {
    "C": CaseSpec("C", "ccpvdz", nelec_act=2, ncore=2,
                  shell_patterns=["2p"], spin=2, sa_weights=[1 / 3] * 3,
                  term=(1, "g"), lane_windows=("2p@1s2s",),
                  scan=ScanSpec(fam=("2p", "3p"),
                                grid=(0.0, 0.8, 17))),
    "Si": CaseSpec("Si", "ccpvdz", nelec_act=2, ncore=6,
                   shell_patterns=["3p"], spin=2, sa_weights=[1 / 3] * 3,
                   term=(1, "g"), lane_windows=("3p@1s2s2p3s",),
                   scan=ScanSpec(fam=("3p", "4p"),
                                 grid=(0.0, 1.1, 23))),
    "Sc": CaseSpec("Sc", "ccpvdz", nelec_act=3, ncore=9,
                   shell_patterns=["3d", "4s"], spin=1,
                   sa_weights=[1 / 5] * 5,
                   term=(2, "g"),
                   lane_windows=("3d4s@ar", "3d3p@semi", "3d3p4d@semi",
                                 "3d3p4s@1s2s2p3s",
                                 "3s3p3d@1s2s2p4s"),
                   lattice_windows=("3d", "3d4s", "3d4d", "3d4s4p",
                                    "3d4s4d", "3d4s4p4d",
                                    "3d4s4p4d5d", "3d4s4p4d4f"),
                   lattice_sector=(2, 1),
                   scan=ScanSpec(fam=("3d", "4d"), grid=(0.0, 0.4, 9),
                                 # 4s fixed: the 4s->5s direction is
                                 # harmful and would mask the d structure
                                 extra_fixed="4s"),
                   ref=ReferenceSpec(
                       occ_closed=("1s", "2s", "2p", "3s", "3p", "4s"),
                       occ_open=("3d",), occ_e_ref=-759.73625510,
                       occ_e_ref_uhf=-759.73789123)),
    "Ti": CaseSpec("Ti", "ccpvdz", nelec_act=4, ncore=9,
                   shell_patterns=["3d", "4s"], spin=2,
                   sa_weights=[1 / 7] * 7,
                   term=(3, "g"),
                   lane_windows=("3d4s@ar", "3d4d@ar", "3d4s4p4d@ar",
                                 "3d3p@semi", "3d3p4d@semi",
                                 "3s3p3d@1s2s2p4s",
                                 "3d3p4s@1s2s2p3s"),
                   scenarios=(_ANCHOR_FIRST,),
                   ref=ReferenceSpec(
                       occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
                       occ_open=("3d",), occ_e_ref=-848.40648875,
                       occ_e_ref_uhf=-848.40991227),
                   planned_windows=("3d4d", "3d4s4p4d"),
                   lattice_windows=("3d", "3d4p", "3d4p4d"),
                   lattice_sector=(2, 0)),
    "U": CaseSpec("U", _U_BASIS, nelec_act=4, ncore=5,
                  shell_patterns=list(_U_TIER["3f2d"]), spin=4,
                  ecp=_U_ECP, tag="U_stuttgart78_dyall",
                  sa_weights=[1 / 17] * 17,  # the ground term's M_L set
                  run_casscf=False,          # the multi-solution class
                  nroots_primary=48,
                  planned_windows=tuple("".join(_U_TIER[t])
                                        for t in ("3f", "4f", "2f2d",
                                                  "4f2d")),
                  ref=ReferenceSpec(
                      occ_closed=("l0", "l1"), occ_open=("l2", "l3"),
                      occ_e_ref=-50.38185233,
                      conv_accelerator="newton")),
}


def scenarios_for(spec: CaseSpec) -> tuple[Scenario, ...]:
    """The spec's effective scenarios: the campaign defaults
    (:data:`DEFAULT_SCENARIOS`) merged with the row's own -- the row
    wins by name.  The ONE evaluated policy object the run gate
    reads; disabled scenarios stay in the tuple (their ``enabled``
    flag is the consumer's check).
    """
    merged = {**DEFAULT_SCENARIOS, **{s.name: s for s in spec.scenarios}}
    return tuple(merged.values())


# the executor -------------------------------------------------------

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
    """Carbon pilot (2p^2, CAS(2,3)): the generic pipeline driven by
    the DECLARED row -- the anatomy (window, spin, SA weights) reads
    off :data:`CAMPAIGN`\\ ['C'], never re-declared here."""
    row = CAMPAIGN["C"]
    return run_atom(
        cfg or row.atom_config(),
        row.shell_patterns,
        row.nelec_act,
        workdir=workdir,
        sa_weights=row.sa_weights,
    )


if __name__ == "__main__":  # pragma: no cover
    run_carbon()
