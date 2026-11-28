"""The table: the routes' closed vocabulary, one row type, the row
store, and the one rendered table.

The store (``vv.json``) and the table (``VV.md``) are the family's
records.  A :class:`Row` carries what a route measured; its total is
COMPUTED (``e_method + E(2)``), never supplied -- the bookkeeping
cannot disagree with itself.  The store accepts only the declared
keys (:func:`table_keys`); anything else is refused, never appended.

Layering
--------
:class:`Route` -- the vocabulary; :class:`Row` -- the record;
``load_rows`` / ``drop_row`` / ``render_table`` -- the store and the
table; ``conv_of`` / ``fmt_range`` / ``de_fmt`` -- the cells' forms.
"""

import json
from dataclasses import asdict, dataclass, fields
from enum import Enum
from pathlib import Path

from vv.reference import rohf_energy

EH_EV = 27.211386245988
HOME = Path(__file__).resolve().parents[1]
ROWS = HOME / "vv.json"
TABLE = HOME / "VV.md"


class Route(Enum):
    """The table's routes -- a closed vocabulary.

    ``FREE`` (``free@pur``): PySCF-direct on the purified core.
    ``PINNED`` (``pinned@pur``): the wfnsym-pinned state average on
    the purified core.  ``RADIAL`` (``radial@pur``): the state
    average constrained to the window's composition (the core stays
    frozen on a rotation-closed span) -- the table's last block, for
    the minimal-window family only.  ``REFERENCE`` (``-``): the
    reference rows, which no route produced.
    """

    REFERENCE = "-"
    FREE = "free@pur"
    PINNED = "pinned@pur"
    RADIAL = "radial@pur"


# the CAS ladder: the PLAIN windows first (the incumbent sources,
# size-ordered), then the @r1 versions (the density-block source
# comparison)
LADDER_PLAIN = ("(4,12)", "(4,15)", "(4,16)", "(6,13)")
LADDER_R1 = ("(4,12)@r1", "(4,15)@r1", "(4,16)@r1", "(4,21)@r1",
             "(6,13)@r1")
# the minimal-window family (the window-variant sweep,
# logs/u14e_vv_variants.log): the free state average's converged
# (4,12)@r1 active block on the key window's purified core, plus the
# channels the suffix names -- one window per active-space size
LADDER_MIN = ("(4,12)m", "(4,13)m+s", "(4,14)m+2s", "(4,15)m+p",
              "(4,16)m+p+s", "(4,17)m+d", "(4,18)m+d+s", "(4,19)m+f",
              "(4,20)m+f+s", "(4,21)m+f+2s", "(4,22)m+f+p",
              "(4,23)m+f+p+s", "(4,24)m+f+d", "(4,25)m+f+d+s",
              "(4,26)m+f+d+2s")
GROUPS = (LADDER_PLAIN, LADDER_R1, LADDER_MIN)
# which ladder groups each route's columns carry: the free routes run
# on every window, the pinned state average on the key windows, the
# constrained one on the minimal family -- every declared key fills
COVERS = {Route.FREE: GROUPS,
          Route.PINNED: (LADDER_PLAIN, LADDER_R1),
          Route.RADIAL: (LADDER_MIN,)}
# the method classes in display order, each with its (method, route)
# columns in display order; a rule row stands between two classes
BLOCKS = (
    ("reference", (("rohf", Route.REFERENCE), ("uhf", Route.REFERENCE))),
    ("casci", (("casci", Route.FREE),)),
    ("casscf", (("casscf", Route.FREE),)),
    ("sa-casscf", (("sa-casscf", Route.FREE),
                   ("sa-casscf", Route.PINNED))),
    ("sa-casscf-radial", (("sa-casscf", Route.RADIAL),)),
)


def table_keys() -> tuple[str, ...]:
    """Every key the table may carry, in display order: the
    references, then each method class column by column -- for each,
    the plain ladder, the @r1 rows, the minimal family."""
    keys = []
    for _cls, columns in BLOCKS:
        for method, route in columns:
            if route is Route.REFERENCE:
                keys.append(f"{method}|-|-")
            else:
                keys += [f"{method}|{route.value}|{cas}"
                         for group in COVERS[route] for cas in group]
    return tuple(keys)


def _class_of(row: "Row") -> str:
    return next(cls for cls, columns in BLOCKS
                if (row.method, row.route) in columns)


@dataclass(frozen=True)
class Row:
    """One table row: what one route measured on one active space.

    ``e_method`` is the energy in the row's own convention (a free
    row: its first root; a pinned SA row: the component ensemble);
    ``e2`` is E(2) in Eh; ``l2`` / ``s2`` are the pre-formatted cells
    over the row's priced states; ``conv`` is the conv cell
    (:func:`conv_of`'s form).

    A reference row carries no E(2) and totals to ``e_method``.  A
    route row with ``e2`` and ``l2`` both ``None`` is WITHDRAWN: its
    priced cells await a re-measurement (``e_total`` is then None).
    """

    method: str
    route: Route
    cas: str
    e_method: float
    e2: float | None
    l2: str | None
    s2: str
    conv: str | None = None

    def __post_init__(self):
        # the withdrawn cells travel together; a reference never prices
        if self.reference:
            if self.e2 is not None:
                raise ValueError(f"reference row {self.key} carries an E(2)")
        elif (self.e2 is None) != (self.l2 is None):
            raise ValueError(
                f"row {self.key}: E(2) and <L2> are withdrawn together "
                f"(got e2={self.e2}, l2={self.l2})")

    @property
    def key(self) -> str:
        """The store key: ``method|route|cas``."""
        return f"{self.method}|{self.route.value}|{self.cas}"

    @property
    def reference(self) -> bool:
        return self.route is Route.REFERENCE

    @property
    def withdrawn(self) -> bool:
        """A route row whose priced cells were withdrawn."""
        return self.e2 is None and not self.reference

    @property
    def e_total(self) -> float | None:
        """E_method + E(2); a reference row totals to E_method; a
        withdrawn row has no total."""
        if self.reference:
            return self.e_method
        return None if self.e2 is None else self.e_method + self.e2

    def to_json(self) -> dict:
        d = asdict(self)
        d["route"] = self.route.value
        return d

    @classmethod
    def from_json(cls, d: dict) -> "Row":
        """The row of one store record; refuses a field outside the
        record schema and a route outside the vocabulary."""
        known = {f.name for f in fields(cls)}
        extra = sorted(set(d) ^ known)
        if extra:
            raise ValueError(
                f"store record {d.get('method')}|{d.get('route')}|"
                f"{d.get('cas')} carries fields {extra} -- the record "
                f"schema is {sorted(known)}")
        return cls(d["method"], Route(d["route"]), d["cas"],
                   d["e_method"], d["e2"], d["l2"], d["s2"], d["conv"])


# --------------------------------------------------------------- the cells

def de_fmt(de_eh: float) -> str:
    """One delta, both units: ``mEh / eV`` (4 decimals)."""
    return f"{de_eh * 1e3:+.4f} / {de_eh * EH_EV:+.4f}"


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


# ------------------------------------------------------- the store + table

def _refuse_unknown(keys) -> None:
    """Fail fast on any key outside the declared table (the keep-set:
    the references, the free routes' CAS ladder, the pinned SA's) --
    an undeclared key used to be appended as a block of its own."""
    bad = [k for k in keys if k not in table_keys()]
    if bad:
        raise ValueError(
            f"row key(s) outside the declared table: {bad} -- the "
            f"grammar is method|route|cas with route in "
            f"{sorted(r.value for r in Route)} and cas on the ladder "
            f"(vv.table.table_keys()); print-only disclosures never "
            f"enter the store")


def load_rows() -> dict[str, Row]:
    """The store's rows by key (an absent store is empty).  Refuses a
    record whose key is undeclared or disagrees with its fields."""
    if not ROWS.exists():
        return {}
    raw = json.loads(ROWS.read_text())
    rows = {k: Row.from_json(r) for k, r in raw.items()}
    _refuse_unknown(rows)
    bad = [k for k, r in rows.items() if r.key != k]
    if bad:
        raise ValueError(f"store record(s) whose key disagrees with "
                         f"its fields: {bad}")
    return rows


def drop_row(row: Row) -> None:
    """Merge one measured row into the store and re-render the table
    (a key outside the declared table is refused -- a route's column
    carries only the ladder groups it covers; nothing is written)."""
    _refuse_unknown([row.key])
    records = {k: r.to_json() for k, r in load_rows().items()}
    records[row.key] = row.to_json()
    ordered = {k: records[k] for k in table_keys() if k in records}
    ROWS.write_text(json.dumps(ordered, indent=1) + "\n")
    render_table()


def _cell(x, fmt=None):
    return "—" if x is None else (x if fmt is None else fmt(x))


def render_table() -> None:
    """Render the store into the table: the references lead (ROHF is
    the zero of every dE), then every declared key present, in
    display order.  Every energy <= 4 decimals; deltas in the
    ``mEh / eV`` slash form.  A withdrawn row renders dashes in its
    withdrawn cells."""
    rows = load_rows()
    e0 = rohf_energy()
    out = [
        "# U14e -- the V&V family",
        "",
        "| method | CAS | &lt;L²&gt; | &lt;S²&gt; | conv "
        "| E_method (Eh) | dE_method (mEh/eV) | E(2) (mEh/eV) "
        "| E_tot (Eh) | dE_tot (mEh/eV) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    prev = None
    for key in table_keys():
        if key not in rows:
            continue
        r = rows[key]
        cls = _class_of(r)
        if prev is not None and cls != prev:
            out.append("|" + " === |" * 10)
        prev = cls
        method = r.method if r.reference else f"{r.method} ({r.route.value})"
        e_tot = r.e_total
        out.append(
            f"| {method} | {r.cas} "
            f"| {_cell(r.l2)} | {r.s2} | {r.conv or '—'} "
            f"| {r.e_method:.4f} | {de_fmt(r.e_method - e0)} "
            f"| {_cell(r.e2, lambda x: f'{x * 1e3:+.4f} / {x * EH_EV:+.4f}')} "
            f"| {_cell(e_tot, lambda x: f'{x:.4f}')} "
            f"| {_cell(None if e_tot is None else e_tot - e0, de_fmt)} |")
    TABLE.write_text("\n".join(out))
