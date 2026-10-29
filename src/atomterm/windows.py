"""Windows: the active-space selection layer.

PURPOSE: the window is the approximation dial of a correlated run --
it fixes which orbitals are correlated and which are frozen into the
core. GOAL: construct windows with a KNOWN rotational-closure
guarantee and materialize them onto the full MO substrate, so exact
term classification (the [H, L] certificates) is possible wherever
the construction promises it.

The shared contract: a window is exactly rotation-closed iff it is
m-complete -- every (n, l) shell present with its whole 2l+1 m-set on
one shared radial. The certificate (operators.window_closure), not
the construction story, is the criterion. The kinds:

* AO windows -- the generic name for the two AO-row vocabularies:
  - shell: contracted valence shells; labels are true (n,l) names.
    Use on contracted bases (cc-pVnZ family).
  - primitive: the exponent-ranked radial lattice beneath the
    shells; labels are shell COUNTERS there. Use on uncontracted
    bases, where shell labels lie about n.
  Both are m-complete by construction.
* scf_mo -- converged MO columns: open MOs + k virtuals ranked by
  the CALLER (campaign metric: descending valence-row S-weight,
  never energy order). NEVER rotation-closed.
* no -- occupation groups of a term-averaged density; closed iff
  the source density was m-scalar.
* custom -- any S-orthonormal columns (your MO selection, a
  hand-built window): build_full_mo takes them directly.

Layering inside: `shell_slots` is the ONE shell walker (lm_sectors,
primitive_chain, the l-vocabulary of select_span, and angular's
analytic path sit on it); `lowdin` / `span_projector` /
`s_metric_roots` / `complete_orthonormal` are the S-metric algebra;
`build_full_mo` and `radial_basis` are the two substrate
builders (same core treatment, different window gauge); Window +
SELECTORS are the registry on top.

SELECTORS maps kind names to constructors; :class:`Window` carries
the construction result and materializes it through
:meth:`Window.columns`; the bounded-slice contract (active columns at
``[ncore : ncore + ncas]``) is owned here. AVAS (:func:`window_avas`)
is the ecosystem cross-validator -- a candidate future kind.
"""

from __future__ import annotations

import warnings
from typing import NamedTuple

import numpy as np

from atomterm.natural import (ci_rdm1, natural_orbitals,
                              select_no_window)

__all__ = [
    "Window",
    "occupied_core_idx",
    "select_shell_ao",
    "select_valence_ao",
    "select_span",
    "primitive_chain",
    "select_primitives",
    "shell_weights",
    "window_avas",
    "span_overlap",
    "span_projector",
    "s_metric_roots",
    "shell_slots",
    "complete_orthonormal",
    "build_full_mo",
    "lm_sectors",
    "l_per_row",
    "radial_basis",
    "radial_mask",
    "window_shell",
    "window_primitive",
    "window_scf_mo",
    "window_no",
    "SELECTORS",
    "lowdin",
]


# ------------------------------------------------------- the S-metric algebra

def span_projector(s_mat: np.ndarray, rows: list[int]) -> np.ndarray:
    """S-metric projector onto the span of selected basis rows.

    Parameters
    ----------
    s_mat : the AO overlap matrix (the metric).
    rows : indices of the rows spanning the target subspace.

    Returns
    -------
    ``P = S[:, o] S_oo^{-1} S[o, :]``. The S-weight of an
    S-orthonormal column ``c`` on the span is ``c^dag P c`` --
    bounded by 1, unlike the bare quadratic form ``c_o^dag S_oo c_o``
    which is unbounded for diffuse bases.
    """
    s_oo = s_mat[np.ix_(rows, rows)]
    return s_mat[:, rows] @ np.linalg.inv(s_oo) @ s_mat[rows, :]


def lowdin(c: np.ndarray, s: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """S-metric symmetric orthonormalization of the columns of ``c``.

    Parameters
    ----------
    c : (n, k) columns (spans, not necessarily independent-free).
    s : (n, n) the metric.
    tol : singularity threshold (see Raises).

    Returns
    -------
    (n, k) columns with ``c^dag s c = I``: with
    ``ov = c^dag s c = U diag(ev) U^T`` the result is
    ``c U diag(ev)^{-1/2} U^T``.

    Raises
    ------
    ValueError
        If ``min(ev) <= tol``: near-linearly-dependent columns. Pad or
        drop those directions (or raise ``tol``) instead of NaN-ing
        everything downstream.
    """
    if c.shape[1] == 0:
        return np.asarray(c)  # degenerate: nothing to orthonormalize
    ov = c.T @ s @ c
    ev, u = np.linalg.eigh(ov)
    if ev.min() <= tol:
        raise ValueError(
            f"Lowdin on a near-linearly-dependent block: min overlap "
            f"eigenvalue {ev.min():.2e} <= tol {tol:.0e}; pad/drop the "
            f"needle directions or raise tol")
    return c @ (u / np.sqrt(ev)) @ u.T


def s_metric_roots(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The root pair ``(S^{1/2}, S^{-1/2})`` of an overlap matrix --
    the isometry between the S metric and the Euclidean one
    (``<a|b>_S = (S^{1/2} a)^dag (S^{1/2} b)``).

    Returns
    -------
    ``(s_sqrt, s_isqrt)``: ``s_isqrt @ s @ s_isqrt = I``. Matrix
    products of AO operators are exact only in an orthonormal basis;
    this pair is the standard transport between the two.
    """
    ev, u = np.linalg.eigh(s)
    return (u * np.sqrt(ev)) @ u.T, (u / np.sqrt(ev)) @ u.T


def complete_orthonormal(v: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Append an S-orthonormal complement to the columns of ``v``.

    Parameters
    ----------
    v : (nao, k) S-orthonormal columns (the substrate so far).
    s : (nao, nao) the metric.

    Returns
    -------
    (nao, nao): ``v`` plus a basis of its S-orthogonal complement.

    Notes
    -----
    ``S^{1/2}`` is an isometry from the S metric to the Euclidean one
    (``<a|b>_S = (S^{1/2} a)^dag (S^{1/2} b)``): project the
    complement in the Euclidean picture (``I - vx vx^dag`` with
    ``vx = S^{1/2} v``, itself orthonormal there) and map it back
    through ``S^{-1/2}``. The singular-value filter (``sv > 1e-8``)
    counts the complement exactly only for a CLEAN ``v`` -- input
    columns that are near-dependent bleed the spectrum into (0, 1)
    and the count goes wrong (keep ``v`` clean at the source;
    :func:`radial_basis` does).
    """
    s_sqrt, s_isqrt = s_metric_roots(s)
    vx = s_sqrt @ v
    resid = np.eye(v.shape[0]) - vx @ vx.conj().T
    p, sv, _ = np.linalg.svd(resid)
    c = s_isqrt @ p[:, sv > 1e-8]
    return np.hstack([v, c])


def _block_projector(cols: np.ndarray, s: np.ndarray) -> np.ndarray:
    """S-metric projector onto an S-orthonormal block:
    ``P = cols cols^dag s`` (``P x`` removes the block's span)."""
    return cols @ cols.conj().T @ s


# -------------------------------------------------------- the shell substrate

def shell_slots(mol):
    """The ONE shell walker: yields ``(shell, l, slot, rows)`` per
    contraction slot of every shell, ``rows`` = the slot's 2l+1
    member rows. The fused-shell addressing (contraction-outer,
    m-inner) lives here and only here; consumers: lm_sectors,
    primitive_chain, select_span's l-vocabulary, and angular's
    analytic construction."""
    loc = mol.ao_loc_nr()
    for ib in range(mol.nbas):
        l = abs(mol.bas_angular(ib))
        width = 2 * l + 1
        for c in range(mol.bas_nctr(ib)):
            base = loc[ib] + c * width
            yield ib, l, c, [base + m for m in range(width)]


def lm_sectors(mol) -> dict[tuple[int, int], list[int]]:
    """AO rows grouped by (l, member) over :func:`shell_slots`.

    Returns
    -------
    ``{(l, m): rows}`` -- the exactly S-orthogonal sectors of a
    one-center spherical basis: same-sector functions of different n
    OVERLAP, different-sector functions never do (angular
    orthogonality). The decomposition that makes per-sector
    treatments (:func:`radial_basis`, the l-vocabulary of
    :func:`select_span`) globally consistent by construction.
    """
    sectors: dict[tuple[int, int], list[int]] = {}
    for _ib, l, _c, rows in shell_slots(mol):
        for m, r in enumerate(rows):
            sectors.setdefault((l, m), []).append(r)
    return sectors


def l_per_row(mol) -> np.ndarray:
    """Angular momentum ``l`` of every AO row -- the shell walker's
    l-grouping as one ``(nao,)`` int array (no label parsing)."""
    out = np.empty(mol.nao_nr(), dtype=int)
    for _ib, l, _c, rows in shell_slots(mol):
        out[rows] = l
    return out


def shell_weights(mol, c: np.ndarray, s: np.ndarray):
    """S-metric angular-momentum (l) weights of each column.

    Parameters
    ----------
    mol : the molecule (the l decomposition comes from
        :func:`lm_sectors` -- the one l-vocabulary, no label parsing).
    c : (nao, k) columns (a window, a converged MO set).
    s : (nao, nao) the metric.

    Returns
    -------
    ``(names, w)``: ``w[i, k]`` is the S-weight of column ``i`` on
    all AOs of angular momentum ``l_k``; ``names = ['s', 'p', ...]``.

    Notes
    -----
    S-orthonormal columns make every row of ``w`` sum to 1 -- valid
    because different-l subspaces are exactly orthogonal; grouping by
    (n, l) instead would NOT partition (same-l shells overlap and
    double-count). l-purity is necessary but NOT sufficient for
    rotational closure -- only the certificates measure closure.
    """
    l_names = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h", 6: "i"}
    groups: dict[int, list[int]] = {}
    for (l, _m), rows in lm_sectors(mol).items():
        groups.setdefault(l, []).extend(rows)
    order = sorted(groups)
    w = np.zeros((c.shape[1], len(order)))
    for k, l in enumerate(order):
        idx = sorted(groups[l])  # row order fixes the summation order
        ss = s[np.ix_(idx, idx)]
        w[:, k] = np.einsum("ij,ij->j", c[idx, :], ss @ c[idx, :]).real
    return [l_names[l] for l in order], w


# --------------------------------------------------------- the AO vocabulary

def select_shell_ao(mol, pattern: str) -> list[int]:
    """AO indices of a named shell, by plain label prefix.

    Parameters
    ----------
    mol : the molecule (its AO labels are truthful, but PySCF's
        ``search_ao_label`` patterns are regex-flavored and
        over-match: ``"2p*"`` also hits ``"2s"`` -- hence prefix
        matching).
    pattern : label prefix, e.g. ``"2p"``.

    Examples
    --------
    On C/cc-pVDZ (the fused s shell occupies rows 0-2):

    >>> select_shell_ao(mol, "2p")
    [3, 4, 5]
    """
    labels = [line.split()[-1] for line in mol.ao_labels()]
    return [i for i, lbl in enumerate(labels) if lbl.startswith(pattern)]


def select_valence_ao(mol, patterns: list[str]) -> list[int]:
    """AO indices of several shells, concatenated in pattern order --
    the concatenation IS the active space (its order is the window
    slot order).

    Parameters
    ----------
    mol : the molecule.
    patterns : label prefixes, one per shell.

    Returns
    -------
    Rotation-closed when every pattern contributes a complete m-set
    of one (n, l) shell.

    Raises
    ------
    ValueError
        If a pattern matches no AO label (the message lists the
        available labels).

    Examples
    --------
    On Sc/cc-pVDZ -- note 4s is row 3: PySCF lists shells by basis
    block (all s contractions first), NOT by energy; label selection,
    never positions:

    >>> select_valence_ao(mol, ["3d", "4s"])
    [21, 22, 23, 24, 25, 3]
    """
    idx: list[int] = []
    for pattern in patterns:
        matches = select_shell_ao(mol, pattern)
        if not matches:
            raise ValueError(
                f"shell pattern {pattern!r} matched no AO labels; "
                f"available: {[line.split()[-1] for line in mol.ao_labels()]}"
            )
        idx += matches
    return idx


def select_span(mol, pattern: str) -> list[int]:
    """AO rows of one span pattern -- the occupation-span vocabulary.

    Parameters
    ----------
    mol : the molecule.
    pattern : ``"l0".."l9"`` selects a WHOLE angular-momentum
        subspace (via :func:`lm_sectors`): the honest vocabulary on
        uncontracted bases, where AO labels are shell counters and
        lie about ``n`` -- only ``l`` itself is trustworthy there.
        Anything else is a label prefix (:func:`select_shell_ao`).

    Raises
    ------
    ValueError
        If a non-``l`` pattern matches no AO label.
    """
    if len(pattern) == 2 and pattern[0] == "l" and pattern[1].isdigit():
        l = int(pattern[1])
        return sorted(r for (ll, m), rws in lm_sectors(mol).items()
                      if ll == l for r in rws)
    rows = select_shell_ao(mol, pattern)
    if not rows:
        raise ValueError(
            f"span pattern {pattern!r} matched no AO labels; available: "
            f"{[line.split()[-1] for line in mol.ao_labels()]}")
    return rows


def primitive_chain(mol, l: int) -> list[tuple[float, list[int]]]:
    """The radial primitives of one angular momentum, as a chain.

    Parameters
    ----------
    mol : the molecule (walked via :func:`shell_slots`).
    l : the angular momentum.

    Returns
    -------
    ``[(exponent, ao_rows)]`` sorted by exponent ASCENDING (diffuse
    end first) -- the honest radial coordinate of an uncontracted
    basis, where AO labels number shells from the TIGHT end starting
    at the valence name.

    Raises
    ------
    ValueError
        On a true contraction slot (a ``bas_ctr_coeff`` column
        combining several exponents has no single exponent); for
        contracted bases use :func:`select_valence_ao` (a contracted
        shell is one pre-converged radial) or MO/NO selection.
    """
    chain: list[tuple[float, list[int]]] = []
    for ib, ll, c, rows in shell_slots(mol):
        if ll != l:
            continue
        exps = np.asarray(mol.bas_exp(ib))
        coeff = mol.bas_ctr_coeff(ib)  # (nprim, nctr)
        nz = np.flatnonzero(np.abs(coeff[:, c]) > 1e-12)
        if len(nz) != 1:
            raise ValueError(
                f"shell {ib} (l={l}) contraction slot {c} combines "
                f"{len(nz)} exponents -- exponent-ranked selection is "
                "undefined for contractions; use select_valence_ao "
                "(a contracted shell is one pre-converged radial) or "
                "MO/NO-based window selection"
            )
        chain.append((float(exps[nz[0]]), rows))
    chain.sort(key=lambda t: t[0])
    return chain


def select_primitives(mol, l: int, count: int | None = None,
                      exponent_band: tuple[float, float] | None = None
                      ) -> list[int]:
    """AO rows of l-primitives selected on the :func:`primitive_chain`
    -- m-complete by construction (each primitive contributes its
    whole 2l+1 m-set): exactly rotation-closed windows.

    Parameters
    ----------
    mol : the molecule.
    l : the angular momentum.
    count : the ``count`` most diffuse primitives (shorthand for the
        diffuse end of the chain). Exactly one of ``count`` /
        ``exponent_band``.
    exponent_band : ``(lo, hi)`` -- every primitive with
        ``lo <= exp <= hi``. The scan protocol: bands are CONTIGUOUS
        sets built from the measured open-MO primitive weights, never
        hand-picked intervals; conclusions require neighboring-band
        stability (admit one more / one fewer primitive on each side
        -> same ordering).

    Raises
    ------
    ValueError
        If not exactly one selection argument is given, the band is
        empty, or ``count`` exceeds the chain length.
    """
    chain = primitive_chain(mol, l)
    if (count is None) == (exponent_band is None):
        raise ValueError("give exactly one of count / exponent_band")
    if exponent_band is not None:
        lo, hi = exponent_band
        picked = [rows for exp, rows in chain if lo <= exp <= hi]
        if not picked:
            raise ValueError(
                f"no l={l} primitive has an exponent in "
                f"{exponent_band!r}; the chain (ascending): "
                f"{[round(e, 6) for e, _ in chain]}")
    else:
        if count > len(chain):
            raise ValueError(
                f"asked for {count} primitives of l={l}; the chain has "
                f"{len(chain)}")
        picked = [rows for _, rows in chain[:count]]
    return sorted(r for rows in picked for r in rows)


# ------------------------------------------------------------- the meters

def span_overlap(c_a: np.ndarray, c_b: np.ndarray, s: np.ndarray) -> float:
    """Smallest principal cosine between two spans in the S metric --
    the identity-change gauge of the rungs decomposition: 1.0 exactly
    when the spans are identical, dropping first when a window's
    identity changes (independent of closure).
    """
    sv = np.linalg.svd(c_a.conj().T @ s @ c_b, compute_uv=False)
    return float(sv.min())


# ------------------------------------------------------- the substrate builders

def occupied_core_idx(mf, ncore: int) -> list[int]:
    """Indices of the ``ncore`` doubly-occupied MOs -- the honest core
    block for every :func:`build_full_mo` / :func:`radial_basis` call.

    Parameters
    ----------
    mf : converged ROHF (its ``mo_occ``).
    ncore : number of core orbitals wanted.

    Returns
    -------
    The first ``ncore`` entries of the doubly-occupied block.

    Raises
    ------
    ValueError
        If the closed block is smaller than ``ncore``.

    Notes
    -----
    Warns when the closed block is NOT ``range(ncore)``: canonical MO
    ORDER is meaningless on near-degenerate uncontracted bases (U:
    the first canonical MOs are virtual junk while the closed block
    sits at 22/61/88/120/143; ``range(ncore)`` silently folds five
    garbage columns, measured max|dE| = 14.5 Eh). On contracted bases
    the two coincide -- this makes the coincidence a checked
    invariant instead of an assumption.
    """
    closed = np.flatnonzero(np.asarray(mf.mo_occ) == 2.0)
    if len(closed) < ncore:
        raise ValueError(
            f"asked for {ncore} core orbitals; ROHF closed block has "
            f"{len(closed)}")
    if not np.array_equal(closed[:ncore], np.arange(ncore)):
        warnings.warn(
            f"core block is not the first {ncore} MOs (closed at "
            f"{closed[:ncore].tolist()}); folding the occupied block",
            stacklevel=2)
    return closed[:ncore].tolist()


def _core_block(mf, core_idx, c_core, s):
    """Resolve (canonical MOs or ``c_core``) and Loewdin the core."""
    core = mf.mo_coeff[:, core_idx] if c_core is None else c_core
    return lowdin(core, s)


def _to_columns(act, n: int) -> np.ndarray:
    """The representation adapter: AO-row indices -> raw AO columns;
    column matrices pass through unchanged."""
    return act if isinstance(act, np.ndarray) else np.eye(n)[:, act]


def _core_window_block(s: np.ndarray, core: np.ndarray,
                       window: np.ndarray) -> np.ndarray:
    """The project-renormalize-assemble treatment (pure numpy). The
    core block must ALREADY be S-orthonormal (:func:`_core_block`'s
    output): project the window columns against its span, re-Loewdin,
    return the stacked block.
    """
    window = window - _block_projector(core, s) @ window
    window = lowdin(window, s)
    return np.hstack([core, window])


def build_full_mo(mf, core_idx: list[int], act,
                  c_core: np.ndarray = None) -> np.ndarray:
    """Complete S-orthonormal MO matrix whose ``[ncore : ncore+ncas]``
    block is the requested window on the requested core -- the object
    CASCI/CASSCF consume.

    Parameters
    ----------
    mf : converged SCF (the overlap; canonical MOs).
    core_idx : MO column indices of the core block.
    act : the WINDOW -- AO-row indices (the m-complete selectors, the
        common case; adapted to raw AO columns inside) or explicit
        S-orthonormal column vectors (a Window's ``c_win``, canonical
        MO columns, any custom window). ONE treatment for both.
    c_core : alternative core columns -- chiefly the exactly-closed
        purified core (:func:`atomterm.core.purified_core_orbitals`).

    Returns
    -------
    (nao, nao) S-orthonormal; the window sits at
    ``[len(core_idx) : len(core_idx) + ncas]``.

    Raises
    ------
    ValueError
        Propagated from :func:`lowdin` on near-linearly-dependent
        core or window blocks.

    Notes
    -----
    The treatment is block-safe in the AO metric
    (:func:`_core_window_block`): the core block is Loewdin'd within
    its own span (span- and density-preserving); the window columns
    are projected against the core span and renormalized -- for raw
    AO columns the same-l, m-diagonal correction; for MO columns on
    an alternative core the honest re-orthogonalization (canonical
    columns are S-orthonormal against the core they were SOLVED in,
    not against a different one); for core-orthogonal columns a
    machine-precision no-op. The complement is appended exactly
    (:func:`complete_orthonormal`).
    """
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    v = _core_window_block(s, _core_block(mf, core_idx, c_core, s),
                           _to_columns(act, n))
    if v.shape[1] == n:
        return v
    return complete_orthonormal(v, s)


def _sgs(cols, base, s):
    """S-metric sequential Gram-Schmidt (two refinement passes) of the
    columns in ``cols`` on top of the S-orthonormal block(s) in
    ``base`` -- member-preserving: column identity survives (a
    symmetric Loewdin would rotate within degenerate member sets).
    """
    out = []
    cur = np.hstack(base) if base else None
    for c in cols:
        c = np.asarray(c).reshape(-1, 1)  # uniform (n, 1) columns
        for _ in range(2):
            if cur is not None:
                c = c - cur @ (cur.conj().T @ s @ c)
            c = c / float(np.real(np.sqrt((c.conj().T @ s @ c)[0, 0])))
        out.append(c)
        cur = c if cur is None else np.hstack([cur, c])
    return out


def radial_basis(mf, core_idx: list[int], act_idx: list[int],
                 c_core: np.ndarray = None):
    """:func:`build_full_mo`'s sibling with the member-preserving
    gauge: the full MO matrix WITH per-column (l, member) labels --
    the labeled substrate :func:`radial_mask` reads.

    Parameters
    ----------
    mf : converged SCF (mol + the overlap).
    core_idx : MO column indices of the core block (same resolution
        and Loewdin as :func:`build_full_mo`).
    act_idx : AO-row indices of the window slots, one per
        (l, member); column j of the window IS the (l, m) of
        ``act_idx[j]`` by construction.
    c_core : alternative core columns; pass the purified core for
        the exactly-closed substrate.

    Returns
    -------
    ``(mo_full, labels)`` -- the full (nao, nao) S-orthonormal MO
    matrix and ``labels[i] = (l, member) | None`` per column.

    Notes
    -----
    Window slots are built with :func:`_sgs` (sequential
    Gram-Schmidt), NOT a symmetric Loewdin: within a degenerate
    member set the Loewdin gauge is an arbitrary rotation --
    span-preserving but identity-reassigning, and here column
    identity IS the product. Sector residuals join by S-weight:
    strong ones with their label (genuine radial residuals ~8e-2;
    core tails ~3e-6), near-linearly-dependent needles (the 1e-4
    admission) unlabeled through :func:`complete_orthonormal` --
    never rotation partners (a strong-only 1e-3 cutoff was tried and
    reverted: the excluded band destabilizes PySCF's SA kernel). A
    labeled column whose sector weight collapsed is demoted to
    unlabeled padding (rare; printed when it fires).

    Examples
    --------
    >>> mo, lab = radial_basis(mf, core_idx, act_idx)
    >>> mask = radial_mask(lab, ncore, ncas)
    """
    mol = mf.mol
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    sec = lm_sectors(mol)
    row_sector = {row: key for key, rows in sec.items() for row in rows}

    core = _core_block(mf, core_idx, c_core, s)

    # window in member order: raw AO minus the core span, sequentially
    # orthogonalized -- column j IS the (l, m) of act_idx[j]
    p_core = _block_projector(core, s)
    window = _sgs([np.eye(n)[:, ao] - p_core @ np.eye(n)[:, ao]
                   for ao in act_idx], [core], s)
    v = np.hstack([core] + window)
    labels: list = [None] * len(core_idx) + [row_sector[ao] for ao in act_idx]

    # strong sector residuals, same member-preserving treatment;
    # the admission threshold separates genuine radials from needles
    p_v = _block_projector(v, s)
    strong, strong_labels = [], []
    for key, rows in sorted(sec.items()):
        block = np.eye(n)[:, rows]
        r = block - p_v @ block
        gr = r.conj().T @ s @ r
        ev, u = np.linalg.eigh(gr)
        keep = ev > 1e-4
        if keep.any():
            strong += [r @ (u[:, j] / np.sqrt(ev[j]))
                       for j in np.where(keep)[0]]
            strong_labels += [key] * int(keep.sum())
    if strong:
        v = np.hstack([v] + _sgs(strong, [v], s))
        labels += strong_labels

    mo = v if v.shape[1] == n else complete_orthonormal(v, s)
    labels += [None] * (mo.shape[1] - len(labels))

    # gauge check + needle demotion: a labeled column whose sector
    # weight collapsed is a needle the threshold could not separate
    # -- demote to unlabeled padding (never a rotation partner)
    demoted = []
    for i, lab in enumerate(labels):
        if lab is None:
            continue
        rows = sec[lab]
        c = mo[rows, i]
        w = float(np.real(c.conj().T @ s[np.ix_(rows, rows)] @ c))
        if w <= 0.9:
            labels[i] = None
            demoted.append((i, lab, round(w, 3)))
    if demoted:
        print(f"radial_basis: needle directions demoted to unlabeled "
              f"padding: {demoted}")
    assert mo.shape[1] == n
    return mo, labels


def radial_mask(labels: list, ncore: int, ncas: int):
    """``uniq_var_indices`` mask allowing ONLY same-(l, member)
    rotations into a window slot: pure RADIAL mixing.

    Parameters
    ----------
    labels : the per-column labels of a :func:`radial_basis`
        substrate (``(l, member) | None``).
    ncore, ncas : the active-space geometry the mask is built for.

    Returns
    -------
    A mask callable for ``constrained_casscf`` (the closure captures
    ``labels``): 3d_m may borrow 4d_m/5d_m and the 4s any s; the
    l-content, the m-membership of every slot, and the core are
    untouchable (core slots carry no label).

    Notes
    -----
    Caveat (measured): exact rotational closure additionally requires
    the radial rotation applied IDENTICALLY to all m partners --
    independent per-member radial mixing leaves [H, Lz] nonzero at
    the member-difference level; SA over the term's M_L set is the
    symmetry-restoring companion rung.
    """
    lab = list(labels)

    def mask(self, nmo, ncore_, ncas_, frozen):
        m = np.zeros((nmo, nmo), dtype=bool)
        for w in range(ncore, ncore + ncas):
            same = np.array([lab[u] is not None and lab[u] == lab[w]
                             for u in range(nmo)])
            same[w] = False
            m[same, w] = True
        return m

    return mask


# ------------------------------------------------------------- the registry

class Window(NamedTuple):
    """A candidate active space in one of two representations.

    ``ao_rows``: AO indices -- m-complete by construction (shell /
    primitive selectors); materialize through :meth:`columns`.
    ``c_win``: orthonormal orbital columns (SCF-MO / NO selectors,
    custom windows); used directly as the active columns. Exactly one
    is not None. ``m_complete``: the construction's closure GUARANTEE
    (shell / primitive: exact; scf_mo: never; no: only if the source
    density was m-scalar -- verify with window_closure, this flag
    records intent).
    """

    label: str
    ao_rows: list[int] | None
    c_win: np.ndarray | None
    m_complete: bool

    def columns(self, mf, ncore: int, closed: list[int] | None = None,
                c_core: np.ndarray = None) -> np.ndarray:
        """Full MO substrate with this window's active columns at
        ``[ncore : ncore + ncas]`` -- the object CASCI/CASSCF consume.

        Parameters
        ----------
        mf : converged SCF.
        ncore : number of core orbitals.
        closed : the OCCUPIED core indices (canonical order is
            meaningless on near-degenerate bases -- defaults to
            :func:`occupied_core_idx` for BOTH representations;
            never ``range(ncore)``).
        c_core : alternative core columns (e.g. the pur substrate).

        Notes
        -----
        Both representations assemble through :func:`build_full_mo`'s
        ONE treatment (its Notes own the details); ``c_win`` columns
        constructed orthogonal to the core (the NO windows) pass
        through the projection unchanged.
        """
        closed = closed if closed is not None else occupied_core_idx(mf, ncore)
        act = self.ao_rows if self.ao_rows is not None else self.c_win
        return build_full_mo(mf, closed, act, c_core=c_core)


def window_shell(mol, mf, ncore: int, *, patterns) -> Window:
    """Contracted-shell window (AO rows) -- :func:`select_valence_ao`
    wrapped as a :class:`Window`; m-complete."""
    return Window(f"shell {patterns}", select_valence_ao(mol, patterns),
                  None, True)


def window_primitive(mol, mf, ncore: int, *, bands) -> Window:
    """Primitive-band window (AO rows), ``bands = {l: count | (lo, hi)}``
    -- :func:`select_primitives` wrapped as a :class:`Window`;
    m-complete per primitive. Bands follow the scan protocol
    (contiguous sets of the measured open-MO primitive weights)."""
    rows: list[int] = []
    for l in sorted(bands):
        sel = bands[l]
        if isinstance(sel, tuple):
            rows += select_primitives(mol, l, exponent_band=sel)
        else:
            rows += select_primitives(mol, l, count=sel)
    return Window(f"primitive {bands}", rows, None, True)


def window_scf_mo(mol, mf, ncore: int, *, open_cols, rank, k: int) -> Window:
    """Raw-SCF-MO window: the open reference MOs plus the ``k``
    top-``rank`` virtuals.

    Parameters
    ----------
    open_cols : MO column indices of the open reference MOs.
    rank : the caller's virtual ordering -- the campaign metric is
        DESCENDING valence-row S-weight (character, never energy
        order: energy-slice picking runs away on TMs).
    k : how many virtuals to admit.

    Notes
    -----
    NEVER rotation-closed: the span certificate measures 7-9 at
    every practical size (the SCF's orbitals are m-delocalized).
    """
    cols = list(open_cols) + list(rank[:k])
    return Window(f"scf_mo +{k}", None, mf.mo_coeff[:, cols], False)


def window_no(mol, mf, ncore: int, *, c_act, dets, vec, n_max: int,
              weights=None, deg_tol: float = 1e-7) -> Window:
    """NO window: term-average roots -> 1-RDM -> top-occupation NOs
    (the NO mathematics lives in :mod:`atomterm.natural`).

    Parameters
    ----------
    c_act, dets, vec : the parent CASCI's active columns,
        determinant basis, and CI vectors.
    n_max : how many occupation groups to keep.
    weights : root weights for the average -- pass the whole term's
        multiplet for an m-scalar density (a single root yields an
        m-split window: legal but NOT rotation-closed).
    deg_tol : occupation-grouping tolerance (scale-aware inside
        :func:`atomterm.natural.select_no_window`).

    Notes
    -----
    Verify every output with operators.window_closure.
    """
    norb = c_act.shape[1]
    gamma = ci_rdm1(dets, vec, norb, weights)
    c_no, occ = natural_orbitals(c_act, gamma)
    c_win, kept = select_no_window(c_no, occ, n_max, deg_tol)
    n_cols = sum(len(g) for g in kept)
    return Window(f"no {n_cols} cols", None, c_win, False)


def window_avas(mf, labels: list[str], ncore: int, threshold: float = None):
    """Window chosen by AVAS (pyscf.mcscf.avas, arXiv:1701.07862) --
    the independent ecosystem selector, used to cross-validate the
    shell windows (a candidate future registry kind).

    Parameters
    ----------
    mf : converged SCF.
    labels : target AO span, pyscf label list, e.g.
        ``["Sc 3d", "Sc 4s"]``.
    ncore : number of core orbitals.
    threshold : projected-overlap cutoff (AVAS default 0.2).

    Returns
    -------
    ``(mo_full, ncas, nelecas)`` with the active span at
    ``[ncore : ncore+ncas]`` -- an interface stranger to the Window
    constructors (raw pyscf output; electron bookkeeping stays ours
    -- pass the sector you want).

    Notes
    -----
    A full-shell-dimension AVAS window IS the shell span exactly --
    closed by construction; a threshold that drops members of an
    m-set gives a rotation-DEFECTIVE span (the negative-control case
    for the closure meters).
    """
    from pyscf.mcscf import avas

    kw = {} if threshold is None else {"threshold": threshold}
    ncas, nelecas, mo = avas.avas(mf, labels, ncore=ncore, **kw)
    nelecas = np.atleast_1d(np.asarray(nelecas)).ravel()
    return mo, int(ncas), nelecas


SELECTORS = {
    "shell": window_shell,
    "primitive": window_primitive,
    "scf_mo": window_scf_mo,
    "no": window_no,
}
