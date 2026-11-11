"""The tables' shared vocabulary: the table constants, the label
helpers, and the write protocols (the verified-move atomic write,
the methods-table edit protocol's marker form).

ONE home for what the scripts' writers share (the methods-table
builder, the regen store writer); the per-table policy -- headings,
column sets, row grouping -- lives with each renderer.  Presentation
only: nothing here touches physics or the store's content.
"""

import math
import re
import sys
from pathlib import Path

from atomterm.windows import NOBLE_CORES, core_shells

# the m-set sizes behind the CAS labels (one l, one column count)
SHELL_COLS = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11}
EV = 27.2114
DASH = "—"


def minus(text: str) -> str:
    """Sign-minus normalization (the tables' typographic convention)."""
    return re.sub(r"(?<!e)-", "−", text)


def core_col(k, spec) -> str:
    """The core column (every table's one grammar): a tagged window's
    own core; an untagged cell's prefix -- the declared chemical shell
    set of the key's EFFECTIVE core (``core_chemistry``: the count
    follows the key's sector -- a promoted-pair window freezes fewer
    pairs), else the noble whose closed count fills ``spec.ncore``
    orbitals; the variants naming their deviation."""
    if k.core_tag is not None:
        return core_label(k.core_tag)
    if spec.core_chemistry:
        n_eff = (2 * spec.ncore + spec.nelec_act
                 - k.sector[0] - k.sector[1]) // 2
        base = spec.core_chemistry.get(n_eff, f"n={n_eff}")
    else:
        base = next((n.capitalize() for n, sh in NOBLE_CORES.items()
                     if sum(SHELL_COLS[s[-1]] for s in sh) == spec.ncore),
                    f"n={spec.ncore}")
    if k.core == "raw":
        return f"{base}·raw"
    if k.core == "frac":
        return f"{base}·frac"
    if k.substrate == "pur":
        return f"{base}·pur"
    return base


def core_label(tag: str) -> str:
    """The core column's README notation: noble names capitalized,
    ``semi`` as itself, explicit lists braced."""
    if tag in NOBLE_CORES:
        return tag.capitalize()
    if tag == "semi":
        return tag
    return "{" + ",".join(core_shells(tag)) + "}"


def cas_label(key, spec=None) -> str:
    """The CAS(nelec, ncas) label: electrons from the key's sector,
    columns from the window's width -- THROUGH THE SHELL LAYER when
    the case declares one, and through the WINDOW COMPRESSION when
    the case declares that (a distilled shell contributes its k
    radials' m-sets -- the CAS label reads the delivered records'
    shapes; the raw band count would read CAS(3,43) for a CAS(3,6)).
    A source-tagged window's k rides the TAG (the budget column must
    state the doubled width -- a k=2 row mislabeled at the k=1 shape
    would hide the double-budget reading)."""
    from atomterm.windows import distilled_ncas

    shells = re.findall(r"\d+[spdfgh]", key.span_label)
    if key.source_tag is not None:
        k = key.source_tag[1]
        cols = distilled_ncas(shells, {sh: k for sh in shells})
    elif spec is not None:
        if getattr(spec, "shell_radials", None) is not None:
            cols = distilled_ncas(shells, spec.shell_radials)
        else:
            from atomterm.windows import realized_patterns
            bands = realized_patterns(shells, spec.shell_bands)
            cols = sum(SHELL_COLS[b[-1]] for b in bands)
    else:
        cols = sum(SHELL_COLS[sh[-1]] for sh in shells)
    return f"({sum(key.sector)},{cols})"


def l2_text(payload: dict) -> str:
    """⟨L²⟩ as one cell: the pinned legs' components (one value when
    they agree, the min–max range when they scramble) or the free
    cells' ground root; the em-dash when the sector sat past the
    dense-meter budget (kernel-only payload)."""
    if "l2" in payload:
        l2 = payload["l2"]
        if max(l2) - min(l2) < 1e-6:
            return f"{sum(l2) / len(l2):.4f}"
        return f"{min(l2):.1f}–{max(l2):.1f}"
    g = payload["ground_l2"]
    return f"{g:.4f}" if g == g else DASH        # NaN check


def comm_cell(payload: dict, field: str, res_field: str) -> str:
    """One symmetry column of a generated table: the dense operator
    norm (whole matrices, plain) where the dense budget allowed it;
    past the budget the stored STATE RESIDUAL (the size-independent
    certificate, when the payload carries one) -- marked ``psi``: the
    commutator applied to that row's converged solution, exact zero
    iff it is a symmetry eigenstate; neither -- the em-dash."""
    v = payload.get(field, float("nan"))
    if v == v:
        return f"{v:.1e}"
    r = payload.get(res_field)
    return f"{r:.1e}ψ" if r is not None else DASH


def conv_cell(payload: dict) -> str:
    """The optimizer rows' convergence column: the final macro-cycle
    |delta e_tot| and the cycle count (``1e-9/7`` -- precision per
    cycle count).  Rows without an optimizer (CI) and stores
    predating the field render the em-dash."""
    de, n = payload.get("conv_dE"), payload.get("cycles")
    if de is None or not n:
        return DASH
    return f"1e{math.floor(math.log10(de))}/{n}" if de > 0 else f"0/{n}"


def payload_of(cells: dict, key) -> dict | None:
    """The key's RUN payload (None for pending and void cells)."""
    p = cells.get(key.address)
    if p is None or "void_reason" in p:
        return None
    return p


def atomic_write_text(path: Path, text: str) -> None:
    """Write through a verified ``.tmp`` move: the file appears whole
    or not at all (a killed writer cannot leave a truncated artifact)."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    assert tmp.read_text(encoding="utf-8") == text
    tmp.replace(path)


def replace_block(src: Path, begin: str, end: str, new_block: str) -> str:
    """The methods-table edit protocol's marker form: splice ``new_block``
    between the section's marker pair, structural-check the splice,
    write through :func:`atomic_write_text`.  Returns the new file
    text; the caller owns the block's own structural assertions.
    """
    text = src.read_text()
    if text.count(begin) != 1 or text.count(end) != 1:
        sys.exit(f"{src}: expected exactly one {begin}/{end} marker pair"
                 " (place them around the generated table on the first"
                 " pass)")
    head, rest = text.split(begin)
    _, tail = rest.split(end)
    new = head + new_block + tail
    assert new.count(begin) == 1 and new.count(end) == 1
    atomic_write_text(src, new)
    return new
