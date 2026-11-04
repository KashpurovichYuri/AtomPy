"""The build_*_tables family's shared vocabulary: the table constants,
the label helpers, and the TABLES edit protocol's marker form.

ONE home for what every emitter formats the same way (the survey/
ladder/lattice/methods tables); the per-table policy -- headings,
column sets, row grouping -- lives with each emitter.  Presentation
only: nothing here touches physics or the store.
"""

import math
import re
import sys
from pathlib import Path

from atomterm.windows import NOBLE_CORES, fold_shells

# the m-set sizes behind the CAS labels (one l, one column count)
SHELL_COLS = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11}
EV = 27.2114
DASH = "—"


def minus(text: str) -> str:
    """Sign-minus normalization (the tables' typographic convention)."""
    return re.sub(r"(?<!e)-", "−", text)


def fold_col(k, spec) -> str:
    """The fold column (every table's one grammar): a tagged window's
    own fold; an untagged cell's prefix -- the noble whose closed
    count fills ``spec.ncore`` orbitals, the variants naming their
    deviation."""
    if k.fold_tag is not None:
        return fold_label(k.fold_tag)
    base = next((n.capitalize() for n, sh in NOBLE_CORES.items()
                 if sum(SHELL_COLS[s[-1]] for s in sh) == spec.ncore),
                f"n={spec.ncore}")
    if k.fold == "raw":
        return f"{base}·raw"
    if k.fold == "frac":
        return f"{base}·frac"
    if k.substrate == "pur":
        return f"{base}·pur"
    return base


def fold_label(tag: str) -> str:
    """The fold column's README notation: noble names capitalized,
    ``semi`` as itself, explicit lists braced."""
    if tag in NOBLE_CORES:
        return tag.capitalize()
    if tag == "semi":
        return tag
    return "{" + ",".join(fold_shells(tag)) + "}"


def cas_label(key) -> str:
    """The CAS(nelec, ncas) label: electrons from the key's sector,
    columns from the window span's summed m-set sizes (an occupation-
    GROUP window names its group count -- ``3d4sno2`` is 2 groups)."""
    if m := re.search(r"no(\d+)$", key.span_label):
        cols = int(m.group(1))
    else:
        cols = sum(SHELL_COLS[l] for l in
                   re.findall(r"\d+([spdfgh])", key.span_label))
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


def replace_block(src: Path, begin: str, end: str, new_block: str) -> str:
    """The TABLES edit protocol's marker form: splice ``new_block``
    between the section's marker pair, structural-check the splice,
    write through a verified ``.md.tmp`` move.  Returns the new file
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
    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    return new
