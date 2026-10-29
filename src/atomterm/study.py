"""The declared study runner: the window x treatment product, one
table (ticket 20; ticket 27 named it -- a "study", never "sweep"/"grid").

A :class:`StudySpec` enumerates the product explicitly -- windows (the
:class:`atomterm.windows.Window` objects), treatments (subset of
"casci" / "relax" / "sa-casscf"; the full treatment axis lives in the
:mod:`atomterm.treatments` registry), and per-treatment budgets -- and
:func:`run_study` executes it cell by cell with three guarantees:

- **Quarantine**: a failing cell is recorded as an error and the study
  continues (a crash in one window must not kill the rest);
- **Budgets**: per-treatment optimizer knobs (e.g. ``max_cycle_macro``)
  so long runs are a declared parameter, not a default;
- **Labels**: every cell reports the window's span certificate and --
  in "full" diagnostics mode -- the exact-eigh ground label (PySCF ci
  vectors straggle inside degenerate clusters and are never a label
  source).

Outputs: ``<out>/<label>_study.json`` (all cells) and a markdown
insert (via :func:`render`) meant for the study markers of the
atom's ``elements/<Atom>/TABLES.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path

import numpy as np

from atomterm.casscf import (
    constrained_casscf,
    run_casscf_core_relax,
    run_casci,
)
from atomterm.detbasis import sector_dimension
from atomterm.natural import (
    ci_rdm1,
    natural_orbitals,
    occupation_groups,
    select_no_window,
)
from atomterm.operators import (
    build_operators,
    ci_in_detbasis_order,
    ground_block,
    window_closure,
)
from atomterm.windows import Window, occupied_core_idx

__all__ = ["StudySpec", "CellResult", "plan_study", "run_study",
           "render", "inject_markers"]

TREATMENTS = ("casci", "relax", "sa-casscf")


@dataclass(frozen=True)
class StudySpec:
    """The declaration of a window x treatment study (ticket 20)."""

    label: str                       # chk/store tag (e.g. "u_14e")
    ncore: int
    na: int
    nb: int
    windows: tuple                   # of atomterm.windows.Window
    treatments: tuple = ("casci", "relax", "sa-casscf")
    sa_weights: tuple | None = None  # default: uniform on the measured block
    budgets: dict = field(default_factory=dict)   # treatment -> mc_opts
    diagnostics: str = "light"       # "light" | "full" (exact-eigh labels)
    verbose: int = 0                 # per-cell optimizer chatter (long cells)
    fci_opts: dict | None = None     # e.g. pspace_size=ndet: the exact-CI
                                     # branch per macro (ticket 23's straggle
                                     # fix; clustered multi-root SA needs it)

    def __post_init__(self):
        bad = set(self.treatments) - set(TREATMENTS)
        assert not bad, f"unknown treatments {bad}; choose from {TREATMENTS}"

    def cell_key(self, wtag: str, treatment: str) -> str:
        return f"study_{self.label}_{wtag}_{treatment}"


@dataclass
class CellResult:
    key: str
    status: str                      # "ok" | "error"
    energy: float | None = None
    l2: float | None = None
    span_lz: float | None = None
    span_lperp: float | None = None
    dets: int | None = None
    converged: bool | None = None
    error: str | None = None


def _wtag(w) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in w.label).strip("_")


def plan_study(spec: StudySpec, mol=None) -> str:
    """The intent manifest: exactly what will run, before anything runs."""
    lines = [f"STUDY {spec.label}: sector (na, nb) = ({spec.na}, {spec.nb});"
             f" treatments = {', '.join(spec.treatments)};"
             f" diagnostics = {spec.diagnostics}"]
    for w in spec.windows:
        ncas = (len(w.ao_rows) if w.ao_rows is not None
                else w.c_win.shape[1])
        ndet = sector_dimension(ncas, spec.na, spec.nb)
        gb = ndet ** 2 * 16 / 2 ** 30
        lines.append(f"  [{w.label}] ncas = {ncas}, dets = {ndet}"
                     f" (eigh matrix ~{gb:.1f} GiB), m_complete ="
                     f" {w.m_complete}"
                     f", budgets = {spec.budgets or '{}'}")
    lines.append("  cells: " + "; ".join(
        f"{_wtag(w)}x{t}" for w in spec.windows for t in spec.treatments))
    return "\n".join(lines)


def run_study(spec: StudySpec, mol, mf, out: str | Path,
              c_core: np.ndarray = None) -> dict:
    """Execute the declared study. Quarantined per cell; JSON out."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    chk = str(out / f"{spec.label}_study.chk")   # per-cell atomterm/ keys
    closed = occupied_core_idx(mf, spec.ncore)
    cells: dict[str, dict] = {}
    for w in spec.windows:
        wtag = _wtag(w)
        ncas = (len(w.ao_rows) if w.ao_rows is not None
                else w.c_win.shape[1])
        # --- window substrate (shared by the window's treatments)
        try:
            mo = w.columns(mf, spec.ncore, closed=closed, c_core=c_core)
            ortho = float(np.abs(mo.conj().T @ mf.get_ovlp() @ mo
                                 - np.eye(mo.shape[1])).max())
            assert ortho < 1e-8, f"substrate not orthonormal ({ortho:.2e})"
            c_act = mo[:, spec.ncore:spec.ncore + ncas]
            wc = window_closure(mol, c_act)
        except Exception as exc:                     # noqa: BLE001 (quarantine)
            for t in spec.treatments:
                key = spec.cell_key(wtag, t)
                cells[key] = vars(CellResult(key=key, status="error",
                                             error=f"substrate: {exc}"))
            continue
        for t in spec.treatments:
            key = spec.cell_key(wtag, t)
            print(f"  [cell] {key} starting", flush=True)
            try:
                if t == "casci":
                    cas = run_casci(mf, ncas, (spec.na, spec.nb), mo,
                                    chk, nroots=1, verbose=spec.verbose)
                elif t == "relax":
                    cas = run_casscf_core_relax(mf, ncas, (spec.na, spec.nb),
                                                mo, chk, verbose=spec.verbose,
                                                mc_opts=spec.budgets.get(t))
                elif t == "sa-casscf":
                    if spec.sa_weights is not None:
                        weights = list(spec.sa_weights)
                    else:
                        _, nroots, wts, _ = ground_block(build_operators(
                            mol, mf, c_act, spec.ncore, spec.na, spec.nb,
                            c_core=c_core))
                        weights = list(wts)
                    cas = constrained_casscf(
                        mf, ncas, (spec.na, spec.nb), chk, mo_full=mo,
                        weights=weights, name=f"study_{wtag}_sa_casscf",
                        verbose=spec.verbose, mc_opts=spec.budgets.get(t),
                        fci_opts=spec.fci_opts)
                else:                            # pragma: no cover
                    raise ValueError(t)
                res = {"key": key, "status": "ok",
                       "energy": float(cas.energies[0]),
                       "dets": int(sector_dimension(ncas, spec.na, spec.nb)),
                       "span_lz": wc["lz"], "span_lperp": wc["lperp"],
                       "converged": bool(getattr(cas.mc, "converged",
                                                 True))}
                if spec.diagnostics == "full":
                    ops = build_operators(mol, mf, c_act, spec.ncore,
                                          spec.na, spec.nb, c_core=c_core)
                    # pyscf FCI ordering -> detbasis ordering before
                    # contracting our operator matrices (they differ
                    # for norb >= 4 with >= 2 electrons in a block)
                    ci = ci_in_detbasis_order(cas.ci[0], ncas,
                                              spec.na, spec.nb)
                    res["l2"] = float(ci @ ops["L2"].real @ ci)
                cells[key] = res
            except Exception as exc:             # noqa: BLE001 (quarantine)
                cells[key] = vars(CellResult(key=key, status="error",
                                             error=repr(exc)))
    (out / f"{spec.label}_study.json").write_text(
        json.dumps({"label": spec.label, "cells": cells}, indent=2),
        encoding="utf-8")
    return cells


def render(cells: dict) -> str:
    """The markdown rows for the atom's T3 study markers (between
    ``<!-- study:begin -->`` and ``<!-- study:end -->``)."""
    lines = ["<!-- study:begin --> (generated by atomterm.study)",
             "",
             "| cell | status | E, Eh | ⟨L²⟩ | ‖[P,L̂z]‖ | ‖[P,L̂⊥]‖ |",
             "|---|---|---|---|---|---|"]
    for key, c in cells.items():
        if c["status"] == "ok":
            l2 = f"{c['l2']:.4f}" if c.get("l2") is not None else "--"
            status = "NOT-CONVERGED" if c.get("converged") is False else "ok"
            lines.append(f"| {key} | {status} | {c['energy']:.6f} | {l2} "
                         f"| {c['span_lz']:.1e} | {c['span_lperp']:.1e} |")
        else:
            lines.append(f"| {key} | ERROR | {c['error'][:60]} | | | |")
    lines.append("")
    lines.append("<!-- study:end -->")
    return "\n".join(lines) + "\n"



def inject_markers(path, block: str, tag: str = "study") -> None:
    """Replace the content between the ``<!-- {tag}:begin/end -->``
    markers of an element TABLES.md (the generated-insert contract;
    ticket 28). Separate studies = separate tags (e.g. "study-six13"),
    each its own section of the same file."""
    import re
    p = Path(path)
    src = p.read_text()
    pat = re.compile(
        rf"(<!-- {tag}:begin -->[^\n]*\n)(.*?)(<!-- {tag}:end -->)", re.S)
    assert pat.search(src), f"no {tag} markers in {p}"
    # render() brackets its output with the same markers (standalone
    # use); strip them here so repeated injects do not multiply them
    inner = "\n".join(
        line for line in block.splitlines()
        if not re.fullmatch(rf"<!-- {tag}:(begin|end).*", line.strip()))
    p.write_text(pat.sub(lambda m: m.group(1) + inner + "\n" + m.group(3),
                         src), encoding="utf-8")
    print(f"[inject] {p} updated")
