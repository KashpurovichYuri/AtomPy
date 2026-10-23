"""The window x treatment x start runner: one declared sweep, one table.

A :class:`SweepSpec` enumerates the product explicitly -- windows (the
:class:`atomterm.windows.Window` objects), treatments (subset of
"casci" / "relax" / "sa"), and per-treatment budgets -- and
:func:`run_sweep` executes it cell by cell with three guarantees:

- **Quarantine**: a failing cell is recorded as an error and the sweep
  continues (a crash in one window must not kill the rest);
- **Budgets**: per-treatment optimizer knobs (e.g. ``max_cycle_macro``)
  so long runs are a declared parameter, not a default;
- **Labels**: every cell reports the window's span certificate and --
  in "full" diagnostics mode -- the exact-eigh ground label (PySCF ci
  vectors straggle inside degenerate clusters and are never a label
  source).

Outputs: ``results/<label>_sweep.json`` (all cells) and a markdown
insert (via :func:`render_insert`) meant for the sweep markers of the
atom's ``elements/<Atom>/TABLES.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path

import numpy as np

from atomterm.casscf import (
    constrained_casscf,
    run_casscf_core_relax,
    run_casci,
)
from atomterm.operators import build_operators, window_closure
from atomterm.windows import occupied_core_idx

__all__ = ["SweepSpec", "CellResult", "plan_sweep", "run_sweep",
           "render_insert", "ground_block"]

TREATMENTS = ("casci", "relax", "sa")


@dataclass(frozen=True)
class SweepSpec:
    """The declaration of a window x treatment sweep (ticket 20)."""

    label: str                       # chk/store tag (e.g. "u_14e")
    ncore: int
    na: int
    nb: int
    windows: tuple                   # of atomterm.windows.Window
    treatments: tuple = ("casci", "relax", "sa")
    sa_weights: tuple | None = None  # default: uniform on the measured block
    budgets: dict = field(default_factory=dict)   # treatment -> mc_opts
    diagnostics: str = "light"       # "light" | "full" (exact-eigh labels)

    def __post_init__(self):
        bad = set(self.treatments) - set(TREATMENTS)
        assert not bad, f"unknown treatments {bad}; choose from {TREATMENTS}"

    def cell_key(self, wtag: str, treatment: str) -> str:
        return f"sweep_{self.label}_{wtag}_{treatment}"


@dataclass
class CellResult:
    key: str
    status: str                      # "ok" | "error"
    energy: float | None = None
    l2: float | None = None
    span_lz: float | None = None
    span_lperp: float | None = None
    dets: int | None = None
    error: str | None = None


def _wtag(w) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in w.label).strip("_")


def ground_block(mol, mf, ncore: int, c_act, na: int, nb: int):
    """The frozen ground-term block of a window, MEASURED: the exact
    eigenspace at ev[0] with its own <L^2>. Returns (ev0, nroots,
    weights) -- the SA anchoring rule (never hardcode a term)."""
    ops = build_operators(mol, mf, c_act, ncore, na, nb)
    ev, vec = np.linalg.eigh(ops["H"].real)
    nchk = min(60, len(ev))
    l2 = [vec[:, r] @ ops["L2"].real @ vec[:, r] for r in range(nchk)]
    block = [r for r in range(nchk)
             if abs(l2[r] - l2[0]) < 1e-6 and ev[r] - ev[0] < 1e-6]
    w = np.zeros(len(ev))
    w[block] = 1.0 / len(block)
    return float(ev[0]), len(block), w


def plan_sweep(spec: SweepSpec, mol=None) -> str:
    """The intent manifest: exactly what will run, before anything runs."""
    lines = [f"SWEEP {spec.label}: sector (na, nb) = ({spec.na}, {spec.nb});"
             f" treatments = {', '.join(spec.treatments)};"
             f" diagnostics = {spec.diagnostics}"]
    for w in spec.windows:
        ncas = (len(w.ao_rows) if w.ao_rows is not None
                else w.c_win.shape[1])
        ndet = _ndet(ncas, spec.na, spec.nb)
        gb = ndet ** 2 * 16 / 2 ** 30
        lines.append(f"  [{w.label}] ncas = {ncas}, dets = {ndet}"
                     f" (eigh matrix ~{gb:.1f} GiB), m_complete ="
                     f" {w.m_complete}"
                     f", budgets = {spec.budgets or '{}'}")
    lines.append("  cells: " + "; ".join(
        f"{_wtag(w)}x{t}" for w in spec.windows for t in spec.treatments))
    return "\n".join(lines)


def _ndet(ncas, na, nb):
    from math import comb
    return comb(ncas, na) * comb(ncas, nb)


def run_sweep(spec: SweepSpec, mol, mf, out: str | Path,
              c_core: np.ndarray = None) -> dict:
    """Execute the declared sweep. Quarantined per cell; JSON out."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    chk = str(out / f"{spec.label}_sweep.chk")   # per-cell atomterm/ keys
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
            try:
                if t == "casci":
                    cas = run_casci(mf, ncas, (spec.na, spec.nb), mo,
                                    chk, nroots=1)
                elif t == "relax":
                    cas = run_casscf_core_relax(mf, ncas, (spec.na, spec.nb),
                                                mo, chk,
                                                **spec.budgets.get(t, {}))
                elif t == "sa":
                    if spec.sa_weights is not None:
                        weights = list(spec.sa_weights)
                    else:
                        _, nroots, wts = ground_block(
                            mol, mf, spec.ncore, c_act, spec.na, spec.nb)
                        weights = list(wts)
                    cas = constrained_casscf(
                        mf, ncas, (spec.na, spec.nb), chk, mo_full=mo,
                        weights=weights, name=f"sweep_{wtag}_sa",
                        **spec.budgets.get(t, {}))
                else:                            # pragma: no cover
                    raise ValueError(t)
                res = {"key": key, "status": "ok",
                       "energy": float(cas.energies[0]),
                       "dets": int(_ndet(ncas, spec.na, spec.nb)),
                       "span_lz": wc["lz"], "span_lperp": wc["lperp"]}
                if spec.diagnostics == "full":
                    ops = build_operators(mol, mf, c_act, spec.ncore,
                                          spec.na, spec.nb, c_core=c_core)
                    ci = np.asarray(cas.ci[0]).ravel()
                    res["l2"] = float(ci @ ops["L2"].real @ ci)
                cells[key] = res
            except Exception as exc:             # noqa: BLE001 (quarantine)
                cells[key] = vars(CellResult(key=key, status="error",
                                             error=repr(exc)))
    (out / f"{spec.label}_sweep.json").write_text(
        json.dumps({"label": spec.label, "cells": cells}, indent=2),
        encoding="utf-8")
    return cells


def render_insert(cells: dict) -> str:
    """The markdown rows for the atom's T3 sweep markers (between
    ``<!-- sweep:begin -->`` and ``<!-- sweep:end -->``)."""
    lines = ["<!-- sweep:begin --> (generated by atomterm.sweep)",
             "| cell | status | E, Eh | <L2> | |[P,Lz]| | |[P,Lp]| |",
             "|---|---|---|---|---|---|"]
    for key, c in cells.items():
        if c["status"] == "ok":
            l2 = f"{c['l2']:.4f}" if c.get("l2") is not None else "--"
            lines.append(f"| {key} | ok | {c['energy']:.6f} | {l2} "
                         f"| {c['span_lz']:.1e} | {c['span_lperp']:.1e} |")
        else:
            lines.append(f"| {key} | ERROR | {c['error'][:60]} | | | |")
    lines.append("<!-- sweep:end -->")
    return "\n".join(lines) + "\n"
