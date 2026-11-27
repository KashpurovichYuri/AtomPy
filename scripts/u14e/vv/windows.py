"""The active spaces: what the routes run on.

An :class:`ActiveSpace` is one active space ASSEMBLED on its core:
``build_full_mo``'s output -- the full S-orthonormal MO matrix in the
partition layout [core | active | virtual], the object CASCI/CASSCF
consume -- with the core's width, the active width and the (na, nb)
sector.  The routes read the core block and the active block off it;
nothing else about the window's construction survives into a route.

A constructor returns a :class:`Built`: the space and its
construction's DISCLOSURE (the lines a leg's log carries -- the
distillation's composition, the core-survival norms).  Two
constructions.  The VISIBLE window is the family's own: the shells'
AO columns, written out and certified elementwise against the
production builder.  The KEY window is the production recipe on a
cell address; its core width and sector are the key's.

Layering
--------
:class:`ActiveSpace` -- the assembled space; :class:`Built` -- the
space and its disclosure; ``window_from_shells`` / ``window_from_key``
(``_key_parts`` / ``_assemble``) -- the two constructions;
``window_6_13_r1`` -- the mixed-source composition; ``lane_window`` --
the lane-core route's substrate and window (the probes' instrument).
"""

import re
from dataclasses import dataclass

import numpy as np

from atomterm.cells import CellKey, key_core, lane_core, window_columns
from atomterm.cases import lane_substrate
from atomterm.operators import window_closure
from atomterm.windows import build_full_mo, lowdin, occupied_core_idx

from vv.reference import Substrate
from vv.symmetry import orthonormal_per_irrep

ORTHO_TOL = 1e-8     # max|C^T S C - I| of the matrix CASCI/CASSCF consume
CLEAN = 1e-10        # a block orthonormal (and orthogonal to the earlier
                    # ones) within this is left untouched by orthonormalize


@dataclass(frozen=True, eq=False)
class ActiveSpace:
    """One active space assembled on its core.

    ``mo`` -- (nao, nao) S-orthonormal, the partition layout
    [core | active | virtual]; ``ncore`` -- the core block's width
    (``ncore`` doubly-occupied orbitals); ``ncas`` -- the active
    width; ``sector`` -- (na, nb) the active electrons.  The core
    block IS the frozen core the routes use: the span the matrix was
    assembled on.
    """

    mo: np.ndarray
    ncore: int
    ncas: int
    sector: tuple[int, int]

    def __post_init__(self):
        # the matrix is shared by every route that runs on the space:
        # an in-place write would silently move a later route's start
        self.mo.setflags(write=False)

    @property
    def core(self) -> np.ndarray:
        return self.mo[:, :self.ncore]

    @property
    def active(self) -> np.ndarray:
        return self.mo[:, self.ncore:self.ncore + self.ncas]

    @property
    def virtual(self) -> np.ndarray:
        return self.mo[:, self.ncore + self.ncas:]

    def on(self, mo, mf) -> "ActiveSpace":
        """The same partition on another orbital set, made exactly
        S-orthonormal (:func:`orthonormalize`) -- a route restarted
        from a converged solution's orbitals."""
        return ActiveSpace(orthonormalize(mo, mf.get_ovlp(), self.ncore,
                                          self.ncas),
                           self.ncore, self.ncas, self.sector)

    def certify(self, mf, tol: float = ORTHO_TOL) -> float:
        """Refuse unless the matrix is S-orthonormal and the partition
        accounts for the molecule's electrons (a non-orthonormal
        matrix fed to an optimizer diverges below the variational
        bound silently).  Returns ``max|C^T S C - I|``.

        Raises
        ------
        ValueError : either clause fails.
        """
        dev = float(np.abs(self.mo.T @ mf.get_ovlp() @ self.mo
                           - np.eye(self.mo.shape[1])).max())
        if dev > tol:
            raise ValueError(
                f"the active space's MO matrix is not S-orthonormal "
                f"(max|C^T S C - I| = {dev:.2e} > {tol:g})")
        n = 2 * self.ncore + sum(self.sector)
        if n != mf.mol.nelectron:
            raise ValueError(
                f"the partition ({self.ncore} core pairs, sector "
                f"{self.sector}) accounts for {n} electrons; the "
                f"molecule has {mf.mol.nelectron}")
        return dev

    def closure(self, mol) -> dict:
        """The one-particle closure meters of the active span:
        ``{"lz": ||[P, Lz]||, "lperp": ||[P, Lperp]||}`` -- both
        vanishing certifies rotation closure completely
        (``operators.window_closure``)."""
        return window_closure(mol, self.active)


@dataclass(frozen=True)
class Built:
    """An active space and its construction's disclosure."""

    space: ActiveSpace
    report: tuple[str, ...]


def orthonormalize(mo, s, ncore: int, ncas: int) -> np.ndarray:
    """An MO matrix in the partition layout made exactly S-orthonormal,
    block by block in the partition's order [core | active |
    virtual]: each block has the earlier blocks projected out and is
    Loewdin-orthonormalized within its own span.

    The optimizer of a symmetry-tagged substrate returns orbitals that
    are orthonormal only to ~4e-05 (observed in the virtual block of
    a pinned state average; the C1 optimizer returns 5e-12).  The
    core's span is kept exactly; a later block's span changes by its
    overlap with the earlier blocks.  A block that is already
    orthonormal within ``CLEAN`` is returned untouched, bitwise.

    Parameters
    ----------
    mo : (nao, nmo) the orbitals, in the partition layout.
    s : (nao, nao) the overlap.
    ncore, ncas : the partition.

    Returns
    -------
    (nao, nmo) with ``C^T S C = I`` to roundoff.

    Raises
    ------
    ValueError : a block is nearly linearly dependent (via ``lowdin``).
    """
    out = np.array(mo, dtype=float)
    for a, b in ((0, ncore), (ncore, ncore + ncas),
                 (ncore + ncas, out.shape[1])):
        if a == b:
            continue
        block = out[:, a:b]
        earlier = out[:, :a]
        cross = earlier.T @ s @ block
        within = block.T @ s @ block - np.eye(b - a)
        if max(np.abs(cross).max(initial=0.0),
               np.abs(within).max(initial=0.0)) <= CLEAN:
            continue          # already orthonormal: returned untouched
        for _ in range(2):    # the projection's roundoff
            block = block - earlier @ (earlier.T @ s @ block)
        out[:, a:b] = lowdin(block, s)
    return out


def orthogonalize(cols, base, s) -> np.ndarray:
    """Columns made S-orthonormal on top of an S-orthonormal
    ``base``: projected against it (two passes), then
    Loewdin-orthonormalized among themselves.  Different (l, member)
    sectors never overlap, so a column only ever mixes with the
    base's columns of its own sector.

    Raises
    ------
    ValueError : nearly linearly dependent (via ``lowdin``).
    """
    out = np.array(cols, dtype=float)
    if base.shape[1]:
        for _ in range(2):
            out = out - base @ (base.T @ s @ out)
    return lowdin(out, s)


def window_with_channels(sub: Substrate, base: Built, blocks=()) -> Built:
    """A base window with correlation channels set beside it.

    Parameters
    ----------
    sub : the reference substrate.
    base : any constructed window (``Built``).
    blocks : ``(columns, report lines)`` per channel, in order --
        :func:`vv.channels.channel_block`'s output.

    Returns
    -------
    The base's active columns followed by the channels', each made
    S-orthonormal on top of what precedes it (:func:`orthogonalize`),
    assembled on the base's core.
    """
    space = base.space
    act, report = space.active, list(base.report)
    s = sub.mf.get_ovlp()
    for cols, lines in blocks:
        act = np.hstack([act, orthogonalize(cols, act, s)])
        report += lines
    return _assemble(sub, space.ncore, act, space.core, space.sector,
                     report)


def dominant_labels(act, mol) -> list[str]:
    """Each column's dominant AO label (its largest raw coefficient)."""
    labels = mol.ao_labels()
    return [labels[int(np.argmax(np.abs(act[:, j])))].split()[-1]
            for j in range(act.shape[1])]


# ------------------------------------------------------ the visible window

def window_from_shells(mol, mf, shells, core, sector) -> Built:
    """The family's own window: the shells' AO columns, the core span
    projected out (the S metric), Loewdin-orthonormalized, the
    complement filled from the canonical virtuals -- written out here
    and certified elementwise against the production builder.

    Parameters
    ----------
    mol, mf : the reference (its canonical MOs fill the complement).
    shells : the window's shells as ``((label prefix, n AOs), ...)`` --
        labels carry member suffixes and are prefix-matched.
    core : (nao, ncore) S-orthonormal core columns (the purified
        span); the core's width is ``ncore``.
    sector : (na, nb) the active electrons.

    Raises
    ------
    ValueError : the label selection is not the declared width, or the
        construction disagrees with ``build_full_mo``.
    """
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    sel = [i for i, lab in enumerate(labels)
           if any(lab.startswith(w) for w, _ in shells)]
    ncas = len(sel)
    if ncas != sum(n for _, n in shells):
        raise ValueError(f"the label selection found {ncas} AOs, "
                         f"expected {sum(n for _, n in shells)} in {shells}")
    nao = mol.nao
    S = mf.get_ovlp()
    act = np.zeros((nao, ncas))
    for k, i in enumerate(sel):
        act[i, k] = 1.0
    P_core = core @ core.T @ S   # the metric projector
    w = act - P_core @ act
    g = w.T @ S @ w
    gv, U = np.linalg.eigh(g)
    c_win = w @ U @ np.diag(gv ** -0.5) @ U.T
    C = np.hstack([core, c_win])
    P_all = C @ C.T @ S   # (C^T S C = I)
    rest = mf.mo_coeff - P_all @ mf.mo_coeff
    g2 = rest.T @ S @ rest
    gv2, U2 = np.linalg.eigh(g2)
    # the complement of an S-orthonormal set is read on the unit
    # scale: its Gram's eigenvalues lie in [0, 1], and an EMPTY
    # complement (the whole space active) leaves only roundoff -- a
    # threshold relative to the maximum would admit that noise
    keep = gv2 > 1e-7
    c_virt = rest @ U2[:, keep] @ np.diag(gv2[keep] ** -0.5)
    if c_virt.shape[1] != nao - C.shape[1]:
        raise ValueError(f"the complement came out {c_virt.shape[1]} "
                         f"wide, expected {nao - C.shape[1]}")
    mo = np.hstack([C, c_virt])
    # the certificate: core + window elementwise against the
    # production builder, the complement by span.  The span band
    # cites the observed diffuse tail's overlap conditioning
    # (1.9e-08); the block equality is the certificate.
    ncore = core.shape[1]
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, ncore), sel,
                            c_core=core)
    nc = C.shape[1]
    dev_blk = np.max(np.abs(mo[:, :nc] - mo_prod[:, :nc]))
    proj = mo_prod[:, nc:]
    resid = mo[:, nc:] - proj @ (proj.T @ S @ mo[:, nc:])
    # no complement (the whole space active): the spans agree by construction
    dev_span = float(np.abs(resid).max()) if resid.size else 0.0
    if dev_blk > 1e-9 or dev_span > 1e-6:
        raise ValueError(
            f"the visible window construction disagrees with the "
            f"production builder (block dev {dev_blk:.2e}, complement "
            f"span dev {dev_span:.2e})")
    line = (f"window: {ncas} AO columns {tuple(w for w, _ in shells)}, "
            f"visible construction == production build_full_mo (block "
            f"dev {dev_blk:.1e}, complement span dev {dev_span:.1e})")
    return Built(ActiveSpace(mo, ncore, ncas, sector), (line,))


# ---------------------------------------------------------- the key window

def _key_parts(sub: Substrate, address: str, radials_ext=None) -> tuple:
    """The production recipe's first two calls on a cell address: the
    key's core resolution (``key_core`` -- the purified span is the
    frozen core) and its window (``window_columns`` -- the declared
    kind).  ``radials_ext`` -- k declarations for shells the case's
    ``shell_radials`` does not carry (a composition's additions live
    in the caller, never in src)."""
    mol, mf, spec = sub.mol, sub.mf, sub.spec
    key = CellKey.from_address(address)
    nc, _c_core, c_pur = key_core(key, mol, mf, spec.shell_bands)
    radials = ({**spec.shell_radials, **radials_ext}
               if radials_ext else spec.shell_radials)
    act, ncas, rep = window_columns(key, mol, mf, spec.shell_bands,
                                    radials, ncol=nc)
    return key, nc, c_pur, act, ncas, rep


def _assemble(sub: Substrate, nc, act, c_pur, sector, report) -> Built:
    """The production recipe's last call: ``build_full_mo`` on the
    key's purified core."""
    mo = build_full_mo(sub.mf, occupied_core_idx(sub.mf, nc), act,
                       c_core=c_pur)
    return Built(ActiveSpace(mo, nc, act.shape[1], sector), tuple(report))


def composition_lines(rep) -> list[str]:
    """The distillation's composition per shell: the distilled k of
    the band chain, the discarded spectrum mass, the capture of the
    head's significant columns."""
    lines = []
    for shell, r in (rep or {}).items():
        caps = [c for group in r["capture"]
                for c in (group if isinstance(group, (list, tuple))
                          else [group]) if c]
        cap = ("capture n/a (no significant source columns)" if not caps
               else f"capture {min(caps):.4f}-{max(caps):.4f}")
        disc = r.get("discarded")
        lines.append(
            f"window {shell}: k={r['k']} of {len(r['bands'])} bands"
            + (f", discarded {disc:.4f}" if disc is not None else "")
            + f", {cap}")
    return lines


def survival_norms(act, c_pur, S) -> np.ndarray:
    """Each column's S-norm after the core span is removed -- 1 for a
    core-orthogonal column, 0 for one inside the core's span (the
    exact per-column remainder norm)."""
    g = c_pur.T @ S @ c_pur
    v = c_pur.T @ S @ act
    proj = np.einsum("ij,ij->j", v, np.linalg.solve(g, v))
    return np.sqrt(np.maximum(np.diag(act.T @ S @ act) - proj, 0.0))


def _survival_line(act, c_pur, S, mol) -> str:
    norms = survival_norms(act, c_pur, S)
    lean = [f"{lab} {x:.4f}"
            for lab, x in zip(dominant_labels(act, mol), norms) if x < 0.99]
    return ("window: the distilled columns' S-norm after removing the "
            "pur core span: " + ("; ".join(lean) if lean else "all at 1.0"))


def window_from_key(sub: Substrate, address: str, radials_ext=None) -> Built:
    """The production recipe on a cell address (``run_cell``'s own
    three calls): ``key_core`` + ``window_columns`` +
    ``build_full_mo``.  The core width and the sector are the KEY's;
    the core is the key's purified span.

    Raises
    ------
    ValueError : the window's certificate refuses the composition
        (propagated from ``window_columns`` -- the declared class).
    """
    key, nc, c_pur, act, _ncas, rep = _key_parts(sub, address, radials_ext)
    report = composition_lines(rep) + [
        _survival_line(act, c_pur, sub.mf.get_ovlp(), sub.mol)]
    return _assemble(sub, nc, act, c_pur, key.sector, report)


# ---------------------------------------------- the (6,13) mixed-source window

def pair_mo(mol, mf, ncore: int) -> tuple[int, np.ndarray]:
    """The reference's own closed s pair outside the ``ncore``-pair
    core: among the occupied columns outside that core set, the one
    dominated by an s band.

    Returns
    -------
    (column index, the column).

    Raises
    ------
    ValueError : the s-dominated column is not unique.
    """
    core = set(occupied_core_idx(mf, ncore))
    occ = [j for j in range(mf.mo_coeff.shape[1])
           if mf.mo_occ[j] > 0 and j not in core]
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    hits = [j for j in occ
            if labels[int(np.argmax(np.abs(np.asarray(mf.mo_coeff)[:, j])))].endswith("s")]
    if len(hits) != 1:
        raise ValueError(
            f"the valence s column is not unique: the occupied "
            f"non-core columns {hits} are s-dominated -- the "
            f"identification needs the character spans, not the labels")
    return hits[0], np.asarray(mf.mo_coeff)[:, hits[0]]


def s_column(act, mol) -> int:
    """The active block's s column: the one whose dominant AO label is
    an s band.

    Raises
    ------
    ValueError : the block's s sector is not exactly one column.
    """
    hits = [j for j, lab in enumerate(dominant_labels(act, mol))
            if re.fullmatch(r"\d+s", lab)]
    if len(hits) != 1:
        raise ValueError(f"the window's s sector is {len(hits)} columns "
                         f"({hits}) -- this composition needs one s radial")
    return hits[0]


def window_6_13_r1(sub: Substrate, address: str) -> Built:
    """The (6,13) MIXED-SOURCE window: the key's ``@r1`` f/d radials
    (the density-block distillation) with the s channel REPLACED by
    the reference's own closed s pair, PROJECTED ONTO ITS PURE-s
    SECTOR.

    No band-limited construction hosts that pair faithfully (it
    spreads over the whole s family); the honest s column is the
    pair's own orbital.  The C1 gauge's pair carries an m-mixed d
    tail, which would break the window's rotation closure, so the
    window takes only the s sector (an l = 0 direction is
    L-invariant).  The replaced column is then S-orthonormalized
    within each dominant D2h irrep on the wrapper: the pinned
    solvers' labeling consumes irrep-pure columns, and a global
    Loewdin would mix across irrep blocks.

    Raises
    ------
    ValueError : :func:`pair_mo`, :func:`s_column` or
        :func:`vv.symmetry.orthonormal_per_irrep` refuses.
    """
    key, nc, c_pur, act, _ncas, rep = _key_parts(sub, address)
    mol, mf = sub.mol, sub.mf
    S = mf.get_ovlp()
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    s_rows = [i for i, lab in enumerate(labels) if re.fullmatch(r"\d+s", lab)]
    d_rows = [i for i, lab in enumerate(labels) if re.match(r"\d+d", lab)]
    j_s, v = pair_mo(mol, mf, nc)
    w_s = float(v[s_rows] @ S[np.ix_(s_rows, s_rows)] @ v[s_rows])
    w_d = float(v[d_rows] @ S[np.ix_(d_rows, d_rows)] @ v[d_rows])
    j_raw = s_column(act, mol)
    old = act[:, j_raw]
    cos2_old = float((v @ S @ old) ** 2 / ((v @ S @ v) * (old @ S @ old)))
    pure = np.zeros_like(v)
    pure[s_rows] = v[s_rows]
    pure = pure / np.sqrt(float(pure @ S @ pure))
    cos2_pure = float((pure @ S @ v) ** 2)
    act = np.array(act, dtype=float)
    act[:, j_raw] = pure
    act, loss = orthonormal_per_irrep(sub.d2h[0], S, act)
    report = composition_lines(rep) + [
        f"s channel: the reference's closed s pair (occupied column "
        f"{j_s}, occ {mf.mo_occ[j_s]:.1f}) replaces the distilled radial "
        f"at {j_raw}; the pair's l-sector weights s {w_s:.6f} / d "
        f"{w_d:.6f} (the m-mixed C1 tail, projected out); the s-pure "
        f"column keeps cos^2 {cos2_pure:.6f} of the pair, whose cos^2 "
        f"against the distilled radial was {cos2_old:.4f}; the "
        f"per-irrep orthonormalization's worst dominant-block weight "
        f"loss {loss:.2e}"]
    return _assemble(sub, nc, act, c_pur, key.sector, report)


# ------------------------------------------------------------- the lane window

def lane_window(spec, address: str) -> tuple[Substrate, Built, np.ndarray]:
    """The lane-core route's substrate and window (the probes'
    instrument; outside the table's keep-set).

    The substrate is the pinned lane's D2h landscape, basin-steered
    (``lane_substrate`` -- an SCF RUNS here); the frozen core is the
    reference's own core columns S-projected (``lane_core``); the
    window is the key's, assembled with ``build_full_mo``.

    Returns
    -------
    (substrate, built, act) : the lane substrate (already D2h-tagged:
        it is its own wrapper); its window carrying the key's sector;
        and the distilled columns BEFORE the assembly's projection
        (the survival diagnostics read them).
    """
    lmol, lmf = lane_substrate(spec)
    key = CellKey.from_address(address)
    _core, nf, c_core, _rep = lane_core(key, lmol, lmf, spec.shell_bands)
    act, ncas, _w = window_columns(key, lmol, lmf, spec.shell_bands,
                                   spec.shell_radials)
    mo = build_full_mo(lmf, list(range(nf)), act, c_core=c_core)
    sub = Substrate(lmol, lmf, spec, d2h=(lmol, lmf))
    return (sub, Built(ActiveSpace(mo, nf, ncas, key.sector),
                       (f"lane window: {address} on a {nf}-column core",)),
            act)
