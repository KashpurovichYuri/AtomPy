"""The row store and the one rendered table."""

import json  # noqa: E402
from pathlib import Path  # noqa: E402

from vv.reference import rohf_energy  # noqa: E402

EH_EV = 27.211386245988

ROWS = Path(__file__).resolve().parents[1] / "vv.json"
TABLE = Path(__file__).resolve().parents[1] / "VV.md"


# the row store + the one table ---------------------------------------

def de_fmt(de_eh: float) -> str:
    """One delta, both units: ``mEh / eV`` (4 decimals)."""
    return f"{de_eh * 1e3:+.4f} / {de_eh * EH_EV:+.4f}"


def row_key(row: dict) -> str:
    """The store key of one row: ``method|route|cas``."""
    return f"{row['method']}|{row['route']}|{row['cas']}"


def _refuse_unknown(keys) -> None:
    """Fail fast on any key outside the declared table (the keep-set:
    the references, the free routes' CAS ladder, the pinned SA's) --
    an undeclared key used to be appended as a block of its own."""
    bad = [k for k in keys if k not in _ORDER]
    if bad:
        raise ValueError(
            f"row key(s) outside the declared table: {bad} -- the "
            f"grammar is method|route|cas with route in {{'-', "
            f"'free@pur', 'pinned@pur'}} and cas on the ladder "
            f"(_ORDER); print-only disclosures never enter the store")


def drop_row(row: dict) -> None:
    """One leg's result row, merged into the family's store
    (``method`` + ``route`` + ``cas`` is the key; a key outside the
    declared table is refused)."""
    _refuse_unknown([row_key(row)])
    rows = json.loads(ROWS.read_text()) if ROWS.exists() else {}
    rows[row_key(row)] = row
    ROWS.write_text(json.dumps(rows, indent=1) + "\n")
    render_table()


# the CAS ladder per method block: the PLAIN windows first (the
# incumbent sources, size-ordered), then the @r1 versions (the
# density-block source comparison)
_LADDER_PLAIN = ("(4,12)", "(4,15)", "(4,16)", "(6,13)")
_LADDER_R1 = ("(4,12)@r1", "(4,15)@r1", "(4,16)@r1", "(4,21)@r1",
              "(6,13)@r1")
_ORDER = tuple(
    ["rohf|-|-", "uhf|-|-"] +
    [f"{m}|free@pur|{c}"
     for m in ("casci", "casscf", "sa-casscf")
     for c in _LADDER_PLAIN + _LADDER_R1] +
    [f"sa-casscf|pinned@pur|{c}" for c in _LADDER_PLAIN + _LADDER_R1])


def fmt_range(vals, nd: int = 4) -> str:
    """One value when the set agrees, ``a..b`` when it scatters."""
    lo, hi = min(vals), max(vals)
    if hi - lo < 10 ** (-nd) / 2:
        return f"{lo:.{nd}f}"
    return f"{lo:.{nd}f}..{hi:.{nd}f}"


def conv_of(precision, n_iters, converged=True) -> str:
    """The conv cell: ``"<precision>/<n_iters>"`` (the user's form).
    ``precision`` -- the achieved |dE| where the optimizer reports it
    (PySCF's callback ``de``; ``run_irrep_sa``'s ``conv_dE``), else
    the requested tolerance; ``n_iters`` -- the optimizer's macro
    cycle count, ``None`` for the CI solves (no iteration
    structure); ``*`` marks a cycle-cap draw (not converged)."""
    star = "" if converged else "*"
    p = "—" if precision is None else f"{precision:.1e}"
    n = "—" if n_iters is None else str(n_iters)
    return f"{p}/{n}{star}"


def render_table() -> None:
    """The small committed table: the references lead (ROHF is the
    zero of every dE), then one row per declared key in table order
    -- the free routes and the pinned SA, both on the pur core.  The
    reading keys live in README · Results; the wave's verdict records
    live in wayfinder/map.md.  Every energy <= 4 decimals; deltas in
    the ``mEh / eV`` slash form.  A RETRACTED row (a cell withdrawn
    pending a re-measurement) renders dashes in the withdrawn cells;
    a key outside the declared table is refused."""
    rows = json.loads(ROWS.read_text()) if ROWS.exists() else {}
    _refuse_unknown(rows)
    order = [k for k in _ORDER if k in rows]
    e0 = rohf_energy()

    def cell(x, fmt=None):
        return "—" if x is None else (x if fmt is None else fmt(x))
    out = [
        "# U14e -- the V&V family",
        "",
        "| method | CAS | &lt;L²&gt; | &lt;S²&gt; | conv "
        "| E_method (Eh) | dE_method (mEh/eV) | E(2) (mEh/eV) "
        "| E_tot (Eh) | dE_tot (mEh/eV) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    # the method-block separators (the campaign tables' === rule
    # rows): the references share one block, each optimizer class
    # reads its own
    _group = {"rohf": "ref", "uhf": "ref"}
    prev = None
    for k in order:
        r = rows[k]
        g = _group.get(r["method"], r["method"])
        if prev is not None and g != prev:
            out.append("|" + " === |" * 10)
        prev = g
        de_var = r["e_method"] - e0
        de_tot = None if r["e_tot"] is None else r["e_tot"] - e0
        e2 = (f"{r['e2_meh']:+.4f} / {r['e2_ev']:+.4f}"
              if r.get("e2_meh") is not None else "—")
        # the render-time consistency check: the store's own bookkeeping
        if r.get("e2_meh") is not None and abs(
                r["e_method"] + r["e2_meh"] * 1e-3 - r["e_tot"]) > 5e-4:
            # 5e-4 -- an order above the table's 4-decimal rendering
            # scale: the store's bookkeeping must vanish well below it
            raise ValueError(
                f"render check: row {k} violates E_tot = E_method + "
                f"E(2) by {r['e_method'] + r['e2_meh'] * 1e-3 - r['e_tot']:.2e} Eh")
        method = (r["method"] if r["route"] == "-"
                  else f"{r['method']} ({r['route']})")
        out.append(
            f"| {method} | {r['cas']} "
            f"| {cell(r['l2'])} | {r['s2']} | {r.get('conv') or '—'} "
            f"| {r['e_method']:.4f} | {de_fmt(de_var)} | {e2} "
            f"| {cell(r['e_tot'], lambda x: f'{x:.4f}')} "
            f"| {cell(de_tot, de_fmt)} |")
    TABLE.write_text("\n".join(out))


def std_row(method: str, route: str, cas: str, e_method: float,
            e2_meh: float | None, e_tot: float, l2: str, s2: str,
            term: str = "—", conv: str | None = None) -> dict:
    """One table row in the family's schema (the deltas are computed
    at render time against the ROHF row; ``l2``/``s2`` pre-formatted;
    ``e_method`` -- the method's own energy, whatever its semantics;
    ``e_tot`` -- the PT2-corrected total, ``e_method + E(2)``;
    ``term`` -- the row's term-identity disclosure, the pick
    column; ``conv`` -- the ``"<precision>/<n_iters>"`` cell,
    ``conv_of``'s form)."""
    return {
        "method": method, "route": route, "cas": cas,
        "l2": l2, "s2": s2,
        "e_method": e_method,
        "e2_meh": e2_meh,
        "e2_ev": None if e2_meh is None else e2_meh * EH_EV / 1e3,
        "e_tot": e_tot,
        "term": term,
        "conv": conv,
    }
