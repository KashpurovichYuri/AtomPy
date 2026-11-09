"""Windows: the active-space selection layer.

PURPOSE
-------
The window is the approximation dial of a correlated run -- it fixes
which orbitals are correlated and which are frozen into the core.
GOAL: construct windows with a KNOWN rotational-closure guarantee
and materialize them onto the full MO substrate, so exact term
classification (the [H, L] certificates) is possible wherever the
construction promises it.

The shared contract
-------------------
A window is exactly rotation-closed iff it is
m-complete -- every (n, l) shell present with its whole 2l+1 m-set on
one shared radial. The certificate (operators.window_closure), not
the construction story, is the criterion.

The kinds
---------

* AO windows -- the generic name for the two AO-row vocabularies:
  - shell: contracted valence shells; labels are true (n,l) names.
    Use on contracted bases (cc-pVnZ family).
  - primitive: the exponent-ranked radial lattice beneath the
    shells; labels are shell COUNTERS there. Use on uncontracted
    bases, where shell labels lie about n.
  Both are m-complete by construction.
* scf_mo -- converged MO columns, two constructions:
  - the character-selected valence window
    (:func:`select_valence_mo`, start='moc'): the SCF's own valence
    + the virtual complement ranked by the named bands'
    S-projection weight. The physics window.
  - the canonical slice (:func:`canonical_order`, start='mo'):
    closed | open | virtual, energy-ordered -- the m-delocalized
    control (the lattice arm; the reference echo).
  NEVER rotation-closed, either of them.
* no -- occupation groups of a term-averaged density; closed iff
  the source density was m-scalar.
* custom -- any S-orthonormal columns (your MO selection, a
  hand-built window): build_full_mo takes them directly.

Layering inside
---------------
`shell_slots` is the ONE shell walker (lm_sectors,
primitive_chain, the l-vocabulary of select_span, and angular's
analytic path sit on it); `lowdin` / `span_projector` /
`s_metric_roots` / `complete_orthonormal` are the S-metric algebra;
`build_full_mo` and `radial_basis` are the two substrate
builders (same core treatment, different window gauge); :class:`Window`
carries a construction result and materializes it through
:meth:`Window.columns`; the bounded-slice contract (active columns at
``[ncore : ncore + ncas]``) is owned here. AVAS (:func:`window_avas`)
is the ecosystem cross-validator.
"""

from __future__ import annotations

import re
import warnings
from typing import NamedTuple

import numpy as np

__all__ = [
    "Window",
    "occupied_core_idx",
    "canonical_order",
    "select_valence_mo",
    "select_shell_ao",
    "select_valence_ao",
    "select_span",
    "shell_realization",
    "realized_patterns",
    "compressed_ncas",
    "window_compressed",
    "distilled_ncas",
    "window_distilled",
    "primitive_chain",
    "select_primitives",
    "fold_shells",
    "closed_electrons",
    "NOBLE_CORES",
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
    "mask_radial_rotations",
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


# ------------------------------------------------------- the fold vocabulary

# the noble-gas closed-core prefixes, in fill order up to and including
# each noble gas's outermost np (universal -- Z-independent)
NOBLE_CORES = {
    "he": ("1s",),
    "ne": ("1s", "2s", "2p"),
    "ar": ("1s", "2s", "2p", "3s", "3p"),
    "kr": ("1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p"),
    "xe": ("1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p",
           "5s", "4d", "5p"),
    "rn": ("1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p",
           "5s", "4d", "5p", "6s", "4f", "5d", "6p"),
}

_SHELL_M = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11}


# --------------------------------------------------------- the shell layer

def shell_realization(shell: str, layer: dict | None = None) -> tuple[str, ...]:
    """The AO bands realizing one chemical shell -- the shell layer.

    The chemical shell (``"3d"``, ``"6s"``) is the primary object; its
    AO realization is BASIS-SPECIFIC and DECLARED, never read off the
    labels: a contracted basis often holds one block per shell (the
    1:1 identity), but multi-zeta shells split it (a chemical 3d over
    two d contractions) and uncontracted bases count radial primitives
    (``"10d"`` is a counter, not a shell). The declaration lives on
    the campaign (:attr:`atomterm.cases.CaseSpec.shell_bands`),
    certified by the measured composition
    (:func:`scripts.band_composition`-style profiles).

    Parameters
    ----------
    shell : the chemical shell name.
    layer : the declared realization mapping (shell -> band tuple);
        ``None`` or an undeclared shell = the 1:1 identity (the shell
        IS its own band -- the contracted-basis behavior, bit-exact
        with the label selectors).

    Returns
    -------
    the shell's bands, in declared order.

    Examples
    --------
    >>> shell_realization("3d")
    ('3d',)
    >>> shell_realization("3d", {"3d": ("12d", "13d")})
    ('12d', '13d')
    """
    if layer is not None and shell in layer:
        bands = tuple(layer[shell])
        if not bands:
            raise ValueError(
                f"the declared realization of shell {shell!r} is empty -- "
                "an unrepresented shell is a missing band, not an "
                "empty tuple")
        stray = [b for b in bands if b[-1] != shell[-1]]
        if stray:
            raise ValueError(
                f"the declared realization of shell {shell!r} names "
                f"other-l bands {stray} -- a chemical shell is realized "
                "by bands of its own angular momentum")
        return bands
    return (shell,)


def realized_patterns(patterns, layer: dict | None = None) -> list[str]:
    """A shell-pattern list expanded through the shell layer -- the
    form the AO selectors consume (``select_valence_ao``-style
    row lookups run on bands, never on chemical names directly).

    Examples
    --------
    >>> realized_patterns(["3d", "4s"])
    ['3d', '4s']
    >>> realized_patterns(["3d", "4s"], {"3d": ("12d", "13d")})
    ['12d', '13d', '4s']
    """
    return [b for p in patterns for b in shell_realization(p, layer)]


def compressed_ncas(patterns, radials) -> int:
    """The column count of a compressed window declaration: the
    summed ``k * (2l+1)`` over the pattern shells.

    The counting rule of :func:`window_compressed` without building
    the columns (the analysis fallbacks and the table renderers count;
    the run builds).

    Raises
    ------
    ValueError
        A pattern shell absent from ``radials`` (the compression is
        declared per shell; a mixed compressed/raw window is not a
        construct this module offers).

    Examples
    --------
    >>> compressed_ncas(["3d", "4s"], {"3d": 1, "4s": 1})
    6
    """
    undeclared = [p for p in patterns if p not in radials]
    if undeclared:
        raise ValueError(
            f"the compression is undeclared for {undeclared} -- every "
            "window shell carries its k (a mixed compressed/raw "
            "window is not a construct this module offers)")
    return sum(radials[p] * (2 * "spdfgh".index(p[-1]) + 1)
               for p in patterns)


def window_compressed(mol, c_ref, s, patterns, layer=None, radials=None,
                      dominance: float = 0.5, tail: float = 1e-8):
    """The COMPRESSED-SHELL window: each chemical shell as its best k
    reference radials x (2l+1) members -- the chemically-sized window
    on a multi-band realization.

    The problem it solves: on an uncontracted (or multi-zeta) basis a
    shell's declared realization spans many primitive bands, and a
    raw AO window over the family counts every band -- the CAS label
    stops speaking chemistry.  The compression takes the span down to
    the k radial directions the reference's own occupied MOs carry
    (the radial the SCF actually uses -- the uncontracted analog of
    one contraction), replicated on every m-member: m-complete by
    construction at every k, the same member-resolved radial Gram
    logic the purified core span uses
    (:func:`atomterm.core.purified_core_orbitals`).

    Parameters
    ----------
    mol : the built Mole (the shell walker's substrate).
    c_ref : ``(nao, n)`` the reference's OCCUPIED MO columns -- the
        radial signal (the shell-dominated subset is selected here).
    s : ``(nao, nao)`` the overlap metric.
    patterns : the window's chemical shell names, in slot order.
    layer : the declared realization mapping (shell -> band tuple),
        :func:`shell_realization`; the FOLD keeps reading the same
        layer at full family width -- the compression is a
        window-side choice only.
    radials : the declared compression depth, shell -> k.  Every
        pattern shell must declare one.
    dominance : the minimum SPAN-projection weight a reference
        column must carry to count as the shell's radial signal (the
        shell's own MO, selected among the occupied columns).  The
        default is permissive by measurement: on a multi-band
        realization the ns valence MO projects ~0.85 on its family
        (its orthogonality nodes live partly on the inner same-l
        families' rows) while every foreign MO sits far below -- the
        identifiability is the gap, not the absolute height.
    tail : the discarded Gram spectrum must sit below this fraction
        of the leading eigenvalue -- the certificate that the
        reference carries k radials and the rest is tail.

    Returns
    -------
    ``(columns, report)`` : the ``(nao, sum_k k*(2l+1))``
    S-orthonormal window columns (radial-major: each radial's full
    m-set adjacent), and the per-shell report (bands, dominated
    columns, their span weights and capture, the Gram spectrum) --
    the composition log's compression entry.

    Raises
    ------
    ValueError
        ``radials`` is None; a pattern shell undeclared; a declared
        band absent from the basis; a shell with no dominated
        reference column (a virtual window shell has no occupied
        radial signal to compress); member chains whose overlap
        blocks differ (the construction's symmetry premise); the
        kept spectrum at/below ``tail`` (k exceeds the radial
        signal) or the discarded spectrum above it (k would discard
        signal); an assembled window that is not S-orthonormal
        (same-l shell blocks overlap -- the member-preserving
        cross-shell treatment is a future seam).

    Notes
    -----
    Different-(l, m) sectors never overlap, so each member's k
    columns are S-orthonormal against every other member's and every
    different-l shell's: the assembled window is S-orthonormal
    without further treatment, and its span is rotation-closed
    (member-complete per radial) whatever the capture -- the capture
    certificate is a quality reading, never a symmetry one.  The
    honest reading is capture <= band weight: the discarded radials
    are gone by declaration.

    Examples
    --------
    On a 1:1 realization (one band per shell, the contracted-basis
    behavior) k=1 reproduces the shell's own m-set: the compressed
    window IS the label-selected one.
    """
    if radials is None:
        raise ValueError(
            "radials is None -- the uncompressed window is "
            "select_valence_ao's business (realized_patterns -> rows)")
    undeclared = [p for p in patterns if p not in radials]
    if undeclared:
        raise ValueError(
            f"the compression is undeclared for {undeclared} -- every "
            "window shell carries its k")
    labels = [line.split()[-1] for line in mol.ao_labels()]
    n = mol.nao_nr()
    blocks: list[np.ndarray] = []
    report: dict[str, dict] = {}
    for shell in patterns:
        k = radials[shell]
        if k < 1:
            raise ValueError(f"shell {shell!r}: k must be >= 1")
        bands = shell_realization(shell, layer)
        l = "spdfgh".index(shell[-1])
        M = 2 * l + 1
        # the member chains: per member m, the shell's band rows in
        # DECLARED band order (one fixed radial-basis order; every
        # member must share it) -- through the shell walker, never
        # label lists (fused shells)
        band_member: dict[str, dict[int, list[int]]] = {b: {} for b in bands}
        for _ib, sl, _c, rows in shell_slots(mol):
            if sl != l:
                continue
            name = labels[rows[0]]
            name = name[: name.index(shell[-1]) + 1]
            if name in band_member:
                for m, r in enumerate(rows):
                    band_member[name].setdefault(m, []).append(r)
        absent = [b for b in bands if not band_member[b]]
        if absent:
            raise ValueError(
                f"the basis has no band(s) {absent} -- the declared "
                f"realization of shell {shell!r}")
        chain = [[r for b in bands for r in band_member[b][m]]
                 for m in range(M)]
        s_r = s[np.ix_(chain[0], chain[0])]
        for m in range(1, M):
            if not np.allclose(s[np.ix_(chain[m], chain[m])], s_r,
                               atol=1e-10):
                raise ValueError(
                    f"shell {shell!r}: the member chains' overlap "
                    "blocks differ -- the member replication premise "
                    "(a spherical one-center basis) is broken")
        # the dominated reference columns: the S-orthogonal
        # projection weight on the shell's family SPAN (bounded [0,1]
        # -- a raw row-block weight is not: same-l families share the
        # l-chain, and a truncated quadratic form overshoots on the
        # node tails)
        rows_flat = sorted(r for m in range(M) for r in chain[m])
        p = span_projector(s, rows_flat)
        norm = np.einsum("ai,ab,bi->i", c_ref, s, c_ref).real
        w = (np.einsum("ai,ab,bi->i", c_ref, p, c_ref).real / norm)
        dominated = np.flatnonzero(w >= dominance)
        if dominated.size == 0:
            raise ValueError(
                f"shell {shell!r}: no reference column carries "
                f">={dominance} S-weight on its bands {list(bands)} "
                f"(best {w.max():.4f}) -- no occupied radial signal "
                "to compress (a virtual shell cannot be compressed "
                "from the reference)")
        # the member-resolved radial Gram: R = [C^0 C^1 ... C^{M-1}],
        # the dominated columns' rows on each member chain
        R = np.hstack([c_ref[chain[m], :][:, dominated] for m in range(M)])
        G = R.conj().T @ (s_r @ R)
        lam, Q = np.linalg.eigh(G)
        lam, Q = lam[::-1].real, Q[:, ::-1]
        if lam[0] <= 0:
            raise ValueError(
                f"shell {shell!r}: the dominated Gram is numerically "
                f"zero (leading eigenvalue {lam[0]:.2e})")
        if lam[k - 1] <= tail * lam[0]:
            raise ValueError(
                f"shell {shell!r}: the declared k={k} exceeds the "
                f"radial signal -- the k-th Gram eigenvalue "
                f"{lam[k - 1]:.2e} is tail against the leading "
                f"{lam[0]:.6f} (the reference carries fewer radials "
                "than declared; lower k)")
        if k < len(lam) and lam[k] > tail * lam[0]:
            raise ValueError(
                f"shell {shell!r}: the discarded spectrum is NOT tail "
                f"-- eigenvalue {k} reads {lam[k]:.6f} against the "
                f"leading {lam[0]:.6f} (tail fraction {tail}); the "
                "compression would discard signal -- raise k")
        E = R @ Q[:, :k] @ np.diag(lam[:k] ** -0.5)   # S_r-orthonormal
        # member replication, radial-major: each radial's full m-set
        # adjacent (k=1 reads exactly like a label-selected m-set)
        block = np.zeros((n, k * M))
        col = 0
        for j in range(k):
            for m in range(M):
                block[chain[m], col] = E[:, j]
                col += 1
        # the capture certificate: each dominated column's
        # S-projection on the compressed span
        cross = block.T @ (s @ c_ref[:, dominated])
        capture = (np.einsum("ij,ij->j", cross, cross).real
                   / norm[dominated])
        blocks.append(block)
        report[shell] = {
            "bands": list(bands), "k": k,
            "dominated": dominated.tolist(),
            "span_weight": w[dominated].round(4).tolist(),
            "capture": capture.round(6).tolist(),
            "spectrum": lam.round(6).tolist(),
        }
    columns = np.hstack(blocks)
    # the assembled certificate: sector orthogonality should make it
    # exact; same-l shell overlap is the undeclared case
    err = np.abs(columns.T @ s @ columns
                 - np.eye(columns.shape[1])).max()
    if err > 1e-9:
        raise ValueError(
            f"the assembled compressed window is not S-orthonormal "
            f"(max deviation {err:.2e}) -- same-l shell blocks "
            "overlap; the member-preserving cross-shell "
            "orthogonalization is a future seam")
    return columns, report


def distilled_ncas(patterns, radials) -> int:
    """The distilled window's active-space width without building:
    each declared shell contributes ``k x (2l+1)``."""
    total = 0
    for shell in patterns:
        if shell not in radials:
            raise ValueError(
                f"the distillation is undeclared for {shell!r} -- "
                "every window shell carries its k")
        l = "spdfgh".index(shell[-1])
        total += radials[shell] * (2 * l + 1)
    return total


def window_distilled(mol, c_ref, s, patterns, layer=None, radials=None,
                     significance: float = 0.05, tail: float = 1e-8):
    """The DISTILLED window: each chemical shell as its declared k
    radial directions, distilled from the reference's occupied MOs by
    the member-summed radial Gram, replicated on the full m-set -- the
    m-complete, rotation-closed window on a multi-band realization.

    GOAL: a chemically-sized window (the CAS label reads k x (2l+1)
    per shell) whose columns carry the radial content the reference's
    own occupied MOs actually use, with the closure guarantee the
    shared contract demands.  The radial directions are data -- the
    leading eigenvectors of the significant columns' member-summed
    Gram on the shell's declared bands -- never a basis rewrite: the
    construction selects a subspace of the existing AO space.

    Parameters
    ----------
    mol : the built Mole (the shell walker's substrate).
    c_ref : ``(nao, n)`` the reference's OCCUPIED MO columns -- the
        radial signal.
    s : ``(nao, nao)`` the overlap metric.
    patterns : the window's chemical shell names, in slot order.
    layer : the declared realization mapping (shell -> band tuple),
        :func:`shell_realization`; the FOLD keeps reading the same
        layer at full family width -- the distillation is a
        window-side choice only.
    radials : the declared depth, shell -> k.  Every pattern shell
        must declare one.
    significance : the minimum SPAN-projection weight (the
        S-orthogonal projection on the shell's band span -- bounded
        ``[0, 1]``, unlike a raw row-block quadratic form) a
        reference column must carry to enter the shell's radial
        signal.  Default permissive by design: the identifiability is
        the GAP between the shell's own MOs and the foreign ones, not
        the absolute height -- a column below the cut shifts nothing
        (its Gram mass lands in the tail, and the tail certificate
        refuses to discard signal).  A permissive cut is what admits
        a hybrid-open valence: a column projecting weakly on its
        home shell still carries the shell's radial.
    tail : the discarded Gram spectrum must sit below this fraction
        of the leading eigenvalue -- the certificate that the
        reference carries exactly k radials and the rest is tail.

    Returns
    -------
    ``(columns, report)`` : the ``(nao, sum_k k*(2l+1))``
    S-orthonormal window columns (radial-major: each radial's full
    m-set adjacent), and the per-shell report (bands, significant
    columns, their span weights and capture, the Gram spectrum).

    Raises
    ------
    ValueError
        ``radials`` is None; a pattern shell undeclared; a declared
        band absent from the basis; a shell with no significant
        reference column (a virtual window shell has no occupied
        radial signal to distill); member chains whose overlap
        blocks differ (the construction's symmetry premise); the
        kept spectrum at/below ``tail`` (k exceeds the radial
        signal) or the discarded spectrum above it (k would discard
        signal); an assembled window that is not S-orthonormal
        (same-l shell blocks overlap -- the member-preserving
        cross-shell treatment is a future seam).

    Notes
    -----
    Different-(l, m) sectors never overlap, so each member's k
    columns are S-orthonormal against every other member's and every
    different-l shell's: the assembled window is S-orthonormal
    without further treatment, and its span is rotation-closed
    (member-complete per radial) whatever the capture -- the capture
    certificate is a quality reading, never a symmetry one.  The
    honest reading is capture <= band weight: the discarded radials
    are gone by declaration.  Per-column labels (which MO fed which
    radial) are gauge; the physical certificates are the subspace
    reads -- the capture sums and the closure meters.

    Examples
    --------
    On a 1:1 realization (one band per shell, the contracted-basis
    behavior) k=1 reproduces the shell's own m-set: the distilled
    window IS the label-selected one.
    """
    if radials is None:
        raise ValueError(
            "radials is None -- the uncompressed window is "
            "select_valence_ao's business (realized_patterns -> rows)")
    undeclared = [p for p in patterns if p not in radials]
    if undeclared:
        raise ValueError(
            f"the distillation is undeclared for {undeclared} -- every "
            "window shell carries its k")
    labels = [line.split()[-1] for line in mol.ao_labels()]
    n = mol.nao_nr()
    blocks: list[np.ndarray] = []
    report: dict[str, dict] = {}
    for shell in patterns:
        k = radials[shell]
        if k < 1:
            raise ValueError(f"shell {shell!r}: k must be >= 1")
        bands = shell_realization(shell, layer)
        l = "spdfgh".index(shell[-1])
        M = 2 * l + 1
        # the member chains: per member m, the shell's band rows in
        # DECLARED band order (one fixed radial-basis order; every
        # member must share it) -- through the shell walker, never
        # label lists (fused shells)
        band_member: dict[str, dict[int, list[int]]] = {b: {} for b in bands}
        for _ib, sl, _c, rows in shell_slots(mol):
            if sl != l:
                continue
            name = labels[rows[0]]
            name = name[: name.index(shell[-1]) + 1]
            if name in band_member:
                for m, r in enumerate(rows):
                    band_member[name].setdefault(m, []).append(r)
        absent = [b for b in bands if not band_member[b]]
        if absent:
            raise ValueError(
                f"the basis has no band(s) {absent} -- the declared "
                f"realization of shell {shell!r}")
        chain = [[r for b in bands for r in band_member[b][m]]
                 for m in range(M)]
        s_r = s[np.ix_(chain[0], chain[0])]
        for m in range(1, M):
            if not np.allclose(s[np.ix_(chain[m], chain[m])], s_r,
                               atol=1e-10):
                raise ValueError(
                    f"shell {shell!r}: the member chains' overlap "
                    "blocks differ -- the member replication premise "
                    "(a spherical one-center basis) is broken")
        # the significant reference columns: the S-orthogonal
        # projection weight on the shell's family SPAN (bounded [0,1]
        # -- a raw row-block weight is not: same-l families share the
        # l-chain, and a truncated quadratic form overshoots on the
        # node tails)
        rows_flat = sorted(r for m in range(M) for r in chain[m])
        p = span_projector(s, rows_flat)
        norm = np.einsum("ai,ab,bi->i", c_ref, s, c_ref).real
        w = (np.einsum("ai,ab,bi->i", c_ref, p, c_ref).real / norm)
        significant = np.flatnonzero(w >= significance)
        if significant.size == 0:
            raise ValueError(
                f"shell {shell!r}: no reference column carries "
                f">={significance} S-weight on its bands {list(bands)} "
                f"(best {w.max():.4f}) -- no occupied radial signal "
                "to distill (a virtual shell has none)")
        # the member-summed radial Gram: R = [C^0 C^1 ... C^{M-1}],
        # the significant columns' rows on each member chain
        R = np.hstack([c_ref[chain[m], :][:, significant] for m in range(M)])
        G = R.conj().T @ (s_r @ R)
        lam, Q = np.linalg.eigh(G)
        lam, Q = lam[::-1].real, Q[:, ::-1]
        if lam[0] <= 0:
            raise ValueError(
                f"shell {shell!r}: the signal Gram is numerically "
                f"zero (leading eigenvalue {lam[0]:.2e})")
        if lam[k - 1] <= tail * lam[0]:
            raise ValueError(
                f"shell {shell!r}: the declared k={k} exceeds the "
                f"radial signal -- the k-th Gram eigenvalue "
                f"{lam[k - 1]:.2e} is tail against the leading "
                f"{lam[0]:.6f} (the reference carries fewer radials "
                "than declared; lower k)")
        if k < len(lam) and lam[k] > tail * lam[0]:
            raise ValueError(
                f"shell {shell!r}: the discarded spectrum is NOT tail "
                f"-- eigenvalue {k} reads {lam[k]:.6f} against the "
                f"leading {lam[0]:.6f} (tail fraction {tail}); the "
                "distillation would discard signal -- raise k")
        E = R @ Q[:, :k] @ np.diag(lam[:k] ** -0.5)   # S_r-orthonormal
        # member replication, radial-major: each radial's full m-set
        # adjacent (k=1 reads exactly like a label-selected m-set)
        block = np.zeros((n, k * M))
        col = 0
        for j in range(k):
            for m in range(M):
                block[chain[m], col] = E[:, j]
                col += 1
        # the capture certificate: each significant column's
        # S-projection on the distilled span
        cross = block.T @ (s @ c_ref[:, significant])
        capture = (np.einsum("ij,ij->j", cross, cross).real
                   / norm[significant])
        blocks.append(block)
        report[shell] = {
            "bands": list(bands), "k": k,
            "significant": significant.tolist(),
            "span_weight": w[significant].round(4).tolist(),
            "capture": capture.round(6).tolist(),
            "spectrum": lam.round(6).tolist(),
        }
    columns = np.hstack(blocks)
    # the assembled certificate: sector orthogonality should make it
    # exact; same-l shell overlap is the undeclared case
    err = np.abs(columns.T @ s @ columns
                 - np.eye(columns.shape[1])).max()
    if err > 1e-9:
        raise ValueError(
            f"the assembled distilled window is not S-orthonormal "
            f"(max deviation {err:.2e}) -- same-l shell blocks "
            "overlap; the member-preserving cross-shell "
            "orthogonalization is a future seam")
    return columns, report


# which noble core ends at n p -- the semi rule's starting prefix per
# open-d period (3d -> ar, 4d -> kr, ...)
_NOBLE_AT_NP = {n: name for n, name in ((3, "ar"), (4, "kr"),
                                         (5, "xe"), (6, "rn"))}


def closed_electrons(shells: tuple[str, ...]) -> int:
    """Electron count of a closed shell set: two per m-slot
    (``2(2l+1)`` per shell).  The count a closed fold removes from the
    atom's active problem.

    Examples
    --------
    >>> closed_electrons(("1s", "2s", "2p"))
    10
    """
    return sum(2 * _SHELL_M[sh[-1]] for sh in shells)


def fold_shells(name: str, open_d: str | None = None) -> tuple[str, ...]:
    """The shells of a named fold -- the closed-core partition an atom
    is split by (the window is what remains).

    Three rule classes (README · Fold notation is the user-facing
    statement of the SAME rule):

    * a noble-gas name (``"ar"``): that noble gas's closed-shell
      prefix (:data:`NOBLE_CORES`, universal);
    * ``"semi"``: the open-(n)d semi-core -- the noble core ending at
      np, minus np, plus (n+1)s.  The 3d metals: Ar - 3p + 4s =
      ``("1s", "2s", "2p", "3s", "4s")``.  Needs ``open_d`` (the
      atom's open d shell, e.g. ``"3d"``);
    * anything else: an EXPLICIT list, concatenated (``"1s2s2p4s"``)
      or comma-separated (``"1s,2s,2p,4s"``) -- it must consume
      entirely.

    Parameters
    ----------
    name : the fold spec (case-insensitive).
    open_d : the atom's open d shell (``"3d"``); only ``"semi"``
      reads it.

    Returns
    -------
    the fold's shells, in fill order.

    Raises
    ------
    ValueError
        unknown name; ``"semi"`` without ``open_d`` or on a non-d
        shell; a list that does not consume entirely.  Whether the
        shells EXIST in a given basis is the consumer's check (a
        missing shell selects no columns).

    Examples
    --------
    >>> fold_shells("ar")
    ('1s', '2s', '2p', '3s', '3p')
    >>> fold_shells("semi", open_d="3d")
    ('1s', '2s', '2p', '3s', '4s')
    """
    key = name.strip().lower()
    if key in NOBLE_CORES:
        return NOBLE_CORES[key]
    if key == "semi":
        if open_d is None or not re.fullmatch(r"[3-6]d", open_d or ""):
            raise ValueError(
                f"fold 'semi' needs the atom's open d shell ('3d'..'6d'), "
                f"got open_d={open_d!r}")
        core = NOBLE_CORES[_NOBLE_AT_NP[int(open_d[0])]]
        n = int(open_d[0])
        return tuple(s for s in core if s != f"{n}p") + (f"{n + 1}s",)
    stripped = key.replace(",", "").replace(" ", "")
    shells = tuple(re.findall(r"\d+[spdfgh]", stripped))
    if not shells or "".join(shells) != stripped:
        raise ValueError(
            f"fold {name!r} is neither a noble-gas name {sorted(NOBLE_CORES)}"
            f" nor an explicit shell list like '1s2s2p4s'")
    return shells


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
    ORDER is meaningless on near-degenerate uncontracted bases (the
    first canonical MOs can be virtual junk while the closed block
    sits elsewhere; ``range(ncore)`` would silently fold garbage
    columns). On contracted bases
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


def canonical_order(mf, ncore: int) -> list[int]:
    """The canonical matrix reordered into the run convention: the
    closed block first (the core prefix), then the remaining occupied
    columns (the open shell, energy-ordered), then the virtuals
    (energy-ordered).

    The start='mo' substrate. PySCF's raw canonical ORDER is the
    m-block listing, not the occupation/energy order -- on uncontracted
    bases the first columns are diffuse virtual junk while the closed
    block sits elsewhere (the warning :func:`occupied_core_idx` fires),
    and ``columns[:ncore]`` would silently fold garbage. On contracted
    bases the permutation is the identity: the coincidence stays a
    checked invariant, not an assumption.

    Examples
    --------
    >>> perm = canonical_order(mf, 5)      # doctest: +SKIP
    >>> mf.mo_coeff[:, perm][:, :5].shape  # the closed block
    (nao, 5)
    """
    occ = np.asarray(mf.mo_occ)
    n_ao = occ.shape[0]
    closed = list(np.flatnonzero(occ == 2.0))
    if len(closed) < ncore:
        raise ValueError(
            f"asked for {ncore} core orbitals; ROHF closed block has "
            f"{len(closed)}")
    open_sh = sorted((i for i in np.flatnonzero((occ > 0) & (occ < 2))),
                     key=lambda i: mf.mo_energy[i])
    virt = sorted((i for i in np.flatnonzero(occ == 0)),
                  key=lambda i: mf.mo_energy[i])
    perm = closed[:ncore] + closed[ncore:] + open_sh + virt
    if perm != list(range(n_ao)):
        warnings.warn(
            f"canonical order scrambled (closed at {closed[:ncore]}"
            f"); reordering into closed|open|virtual",
            stacklevel=2)
    return perm


def select_valence_mo(mf, patterns: list[str], ncore: int) -> list[int]:
    """Canonical columns of the window's character, in the run
    convention -- the start='moc' substrate (the character-selected
    valence window; ``canonical_order``'s slice is the
    m-delocalized control).

    The column permutation placing the window at ``[ncore :
    ncore + ncas]``: the core prefix (:func:`occupied_core_idx`), the
    SCF's own valence -- the closed tail (a promoted pair's columns)
    then the open shells, energy-ordered -- then the virtual
    complement ranked by the named bands' S-projection weight (ties
    broken by energy); the unselected virtuals trail, energy-ordered.

    The count contract is the AO window's: ``ncas`` =
    ``len(select_valence_ao(mol, patterns))`` -- the same label
    declares the same window size on both starts; the columns are
    the orbitals the raw bands are the material of. The bands, not
    the energy order, drive the ranking: the energy-ordered slice is
    diffuse-dominated on uncontracted bases -- the reference echo
    (the span contains the ROHF's own broken determinant; singles
    are Brillouin-dead; the CASCI reads the mean field back), not a
    valence window.

    Returns
    -------
    The full column permutation -- ``mf.mo_coeff[:, perm]`` carries
    the window at ``[ncore : ncore + ncas]``.

    Raises
    ------
    ValueError
        If the SCF's own valence already exceeds the declared window
        size (the label underdeclares its own electrons), or a
        pattern matches no AO (via :func:`select_valence_ao`).

    Examples
    --------
    >>> perm = select_valence_mo(mf, ["10f", "13d"], 5)  # doctest: +SKIP
    >>> perm[5:5 + 12]  # the window: 4 open shells + 8 f/d virtuals
    [...]              # doctest: +SKIP
    """
    mol = mf.mol
    ncas = len(select_valence_ao(mol, patterns))
    core = occupied_core_idx(mf, ncore)
    core_set = set(core)
    occ = np.asarray(mf.mo_occ)
    closed_tail = sorted((i for i in np.flatnonzero(occ == 2.0)
                          if i not in core_set),
                         key=lambda i: mf.mo_energy[i])
    opens = sorted((i for i in np.flatnonzero((occ > 0) & (occ < 2))),
                   key=lambda i: mf.mo_energy[i])
    valence = closed_tail + opens
    n_virt = ncas - len(valence)
    if n_virt < 0:
        raise ValueError(
            f"the SCF's own valence ({len(valence)} columns: closed "
            f"tail {closed_tail} + opens {opens}) already exceeds the "
            f"declared window size ncas = {ncas}")
    band_rows = [mu for pattern in patterns
                 for mu in select_shell_ao(mol, pattern)]
    # the named bands' S-projection weight (the shell_weights metric
    # at band resolution): the row-block quadratic form -- canonical
    # columns spread over their same-l radials, so the valence
    # virtuals outscore the diffuse junk
    s = mf.get_ovlp()
    ss = s[np.ix_(band_rows, band_rows)]
    cb = mf.mo_coeff[band_rows, :]
    score = np.einsum("ij,ij->j", cb, ss @ cb).real
    virt = sorted((i for i in np.flatnonzero(occ == 0)),
                  key=lambda i: mf.mo_energy[i])
    picked = sorted(virt, key=lambda i: (-score[i], mf.mo_energy[i]))
    rest = picked[n_virt:]
    picked = picked[:n_virt]
    perm = (core + valence
            + sorted(picked, key=lambda i: mf.mo_energy[i])
            + sorted(rest, key=lambda i: mf.mo_energy[i]))
    return perm


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
    the labeled substrate :func:`mask_radial_rotations` reads.

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
    >>> mask = mask_radial_rotations(lab, ncore, ncas)
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


def mask_radial_rotations(labels: list, ncore: int, ncas: int):
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
    symmetry-restoring companion treatment.
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


# ------------------------------------------------------- the Window carrier

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


def window_avas(mf, labels: list[str], ncore: int, threshold: float = None):
    """Window chosen by AVAS (pyscf.mcscf.avas, arXiv:1701.07862) --
    the independent ecosystem selector, used to cross-validate the
    shell windows.

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
