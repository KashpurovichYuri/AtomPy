"""The basis-preparation layer: the contraction recipe (stage 0b).

The universal window construction's preparation half -- everything
that turns (atom, primitive basis) into RADIAL COLUMNS before any
treatment runs over them.  The treatment comparison (ROHF -> CASCI ->
CASSCF -> SA) consumes the prepared basis as one frozen input; within
a table the basis is identical by construction, the source axis
riding the cell address (``@r<k>`` / ``@s<k>`` window tags,
:attr:`atomterm.cells.CellKey.source_tag`), never a silent per-cell
re-construction.

The recipe: a chemical shell's AO realization is DECLARED (the
certificate-measured band lists, :func:`atomterm.windows.shell_realization`);
the radial prior is the CONFIGURATION LADDER -- the forced-occupation
reference series over the atom's valence configurations -- averaged as
per-configuration member-summed density blocks; k radial directions
come from the s_r-metric eigenproblem on the averaged block, top-k,
replicated on the full m-set (m-complete, rotation-closed whatever the
capture).  k is the zeta dial: k=1 is a single contraction per shell,
k=2 adds the configuration-spread radial (single-configuration blocks
are rank-1; the ladder average carries the spread -- the spectrum
certifies whether k=2 is reachable).

Two distillation families, stated honestly:

* the INCUMBENT path (:func:`atomterm.windows.window_distilled`) --
  the reference-source, span-signal family the stored cells ran;
  untagged windows, frozen (their identity is the store's);
* the PRIOR axis here -- the density-block family: the
  occupation-weighted member-summed block of ONE configuration
  (``@r<k>``, the head) or the ladder average (``@s<k>``), one
  eigenproblem on the block, k per shell from the tag.

The families share the member chains and the m-replication
(:func:`atomterm.windows.member_chains`); they differ in the radial
signal's definition -- the source axis -- and never mix within one
comparison row.

Conventions
-----------
- The ladder's HEAD configuration is the case's own reference (the
  trusted cache): the head flows in as a converged mf, never a fresh
  draw.  Every further configuration runs character-forced ROHF with a
  MANDATORY basin gate (the degenerate-open chaos is measured; an
  ungated draw is a wrong-basin risk the gate refuses).
- Ladder chks live under ``<workdir>/<Key>/rungs/ladder/`` and are
  TRUSTED once present (the reproducibility unit, as every chk in
  this project).
- The report DISCLOSES the budget: per shell the spectrum, the kept
  k, and the discarded fraction (k=2 on a 2-band chain discards
  nothing -- the cap binds only on wider realizations).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from atomterm.reference import ReferenceSpec, run_rohf
from atomterm.windows import member_chains

if TYPE_CHECKING:  # the mol/mf annotations (no cycle; system is leaf)
    from pyscf import scf as _scf

    from atomterm.system import AtomConfig

__all__ = [
    "LadderSpec",
    "run_ladder",
    "valence_block",
    "metric_radials",
    "series_window",
]

# the kept-signal tail: the k-th eigenvalue must carry this fraction of
# the leading one (a below-tail k-th eigenvalue is numerical noise --
# a rank-1 source asked for k=2; window_distilled carries the same
# certificate at the same default)
_TAIL = 1e-4


# the declaration ------------------------------------------------------

class LadderSpec(NamedTuple):
    """One configuration of the ladder: the forced-occupation recipe
    and its basin gate.

    Attributes
    ----------
    tag : the configuration's label (the chk name and the report key,
        e.g. ``"d1s2"``).
    closed / open : the doubly-/singly-occupied AO span patterns
        (:class:`atomterm.reference.ReferenceSpec` semantics -- by
        character, never energy order).
    spin : the configuration's 2S (the :class:`atomterm.system.AtomConfig`
        spin of its sector).
    e_ref : the basin gate energy (Eh).  MANDATORY: the gate matches a
        fresh draw to the declared basin within 1e-4 (guess cycling),
        refusing the wrong-basin draws the degenerate open shell
        admits.  Measured by the census, then trusted.

    Notes
    -----
    The ladder's FIRST entry is the case's own reference configuration
    (the head): its converged mf enters :func:`run_ladder` directly
    from the trusted cache, so the head's ``e_ref`` is the case's
    declared ``occ_e_ref`` -- recorded for completeness, never
    re-gated.
    """

    tag: str
    closed: tuple[str, ...]
    open: tuple[str, ...]
    spin: int
    e_ref: float


# the ladder -----------------------------------------------------------

def run_ladder(cfg: "AtomConfig", workdir: str | Path, key: str, ladder,
               head_mf) -> dict[str, "_scf.RHF"]:
    """Run (or load, trusted) the configuration ladder.

    Parameters
    ----------
    cfg : the atom (:class:`atomterm.system.AtomConfig`) -- the case's
        own (the ladder varies only the occupations and spin).
    workdir : the elements root.
    key : the campaign key (the store dir name; chks land under
        ``<workdir>/<key>/rungs/ladder/<tag>.chk``).
    ladder : the declared :class:`LadderSpec` tuple; its HEAD (first
        entry) must be the case's reference configuration.
    head_mf : the case's trusted converged reference mf (the head's
        converged state -- never re-run, never re-gated).

    Returns
    -------
    tag -> the converged mf per configuration, the head included.

    Raises
    ------
    ValueError
        the head's spin mismatches the case's (the head IS the case
        reference; a mismatched declaration is a wrong ladder).
    RuntimeError
        a fresh configuration did not reach its declared basin (the
        gate's own error, propagated).
    """
    if ladder[0].spin != cfg.spin:
        raise ValueError(
            f"the ladder head {ladder[0].tag!r} declares spin "
            f"{ladder[0].spin} but the case's reference (which the head "
            f"must be) carries spin {cfg.spin} -- a mismatched ladder "
            "declaration")
    out = {ladder[0].tag: head_mf}
    chk_dir = Path(workdir) / key / "rungs" / "ladder"
    chk_dir.mkdir(parents=True, exist_ok=True)
    for spec in ladder[1:]:
        cfg_x = type(cfg)(cfg.symbol, cfg.basis, spin=spec.spin)
        ref = ReferenceSpec(occ_closed=spec.closed, occ_open=spec.open,
                            occ_e_ref=spec.e_ref)
        _m, mf_x = run_rohf(cfg_x, str(chk_dir / f"{spec.tag}.chk"), ref=ref)
        out[spec.tag] = mf_x
    return out


# the radial signal ----------------------------------------------------

def valence_block(mf, chain, ncol: int) -> np.ndarray:
    """One configuration's member-summed shell block of the VALENCE
    density -- the radial signal of the prior axis.

    Occupation-weighted (the configuration's own physics: a doubly
    occupied shell counts twice) and member-summed (m-scalar by
    construction).  The core columns never enter -- the measured
    core-dominance rejection: a core-inclusive block is the core's
    radial, not the valence's.

    Parameters
    ----------
    mf : the configuration's converged SCF object.
    chain : the shell's member chains
        (:func:`atomterm.windows.member_chains`).
    ncol : the fold boundary (core orbital pairs; the valence columns
        start there).

    Returns
    -------
    the symmetric real ``(len(chain[0]), len(chain[0]))`` block.
    """
    cv = mf.mo_coeff[:, ncol:]
    wv = np.asarray(mf.mo_occ)[ncol:]
    d = (cv * wv) @ cv.conj().T
    b = sum(d[np.ix_(chain[m], chain[m])] for m in range(len(chain)))
    return 0.5 * (b + b.T).real


def metric_radials(block: np.ndarray, s_r: np.ndarray, k: int
                   ) -> tuple[np.ndarray, np.ndarray]:
    """The capture-maximizing top-k radials of a density block: the
    generalized eigenproblem ``(block, s_r^-1)`` -- equivalently
    ``eigh(s_r^{1/2} block s_r^{1/2})`` with vectors transported back
    to s_r-orthonormal chain space.

    The functional being maximized is each significant column's
    S-projection on the radial's span (the capture -- the same
    quantity :func:`atomterm.windows.window_distilled`'s column-space
    Gram maximizes; the two constructions agree on a
    single-configuration source, the consistency certificate).  The
    transport: the generalized problem ``(block, s_r^-1)`` reads
    ``eigh(L^T block L)`` in Cholesky coordinates, its eigenvector
    mapped back through ``s_r^-1`` (``L^-T q``) -- the
    capture-maximizing radial in chain coordinates.  An ``s_r^{-1/2}``
    conjugation instead (``eigh(L^-1 block L^-T)``) mixes the
    Euclidean coefficient pairing into the numerator and returns
    s_r^{-1}-distorted radials -- the planted-radial anchor pins the
    direction.

    Returns
    -------
    ``(v, lam)`` : the ``(width, k)`` s_r-orthonormal radial columns
    (descending eigenvalue order) and the full descending spectrum
    (the k certificate: single-configuration blocks are rank-1, the
    ladder average carries the spread).

    Raises
    ------
    AssertionError
    the returned columns fail s_r-orthonormality (a construction
    error, never a physical one).
    """
    l = np.linalg.cholesky(s_r)
    li = np.linalg.inv(l)
    t = l.T @ block @ l
    t = 0.5 * (t + t.T)
    lam, q = np.linalg.eigh(t)
    lam, q = lam[::-1].real, q[:, ::-1].real
    v = li.T @ q[:, :k]                    # s_r-orthonormal columns
    err = np.abs(v.T @ s_r @ v - np.eye(k)).max()
    assert err < 1e-9, f"s_r-orthonormality {err:.2e}"
    return v, lam


# the window -----------------------------------------------------------

def series_window(mol, mfs: list, s: np.ndarray, patterns: list[str],
                  layer: dict | None, radials: dict[str, int], ncol: int
                  ) -> tuple[np.ndarray, dict]:
    """The prior-axis distilled window: each chemical shell as its
    declared k radial directions from the DENSITY-BLOCK family -- the
    head configuration's block (``@r``) or the ladder average
    (``@s``; the caller passes the mf list) -- replicated on the full
    m-set.

    Parameters
    ----------
    mol : the built Mole (the shell walker's substrate).
    mfs : the source configurations' converged mf objects (ONE = the
        ``@r`` head form; several = the ladder average, the ``@s``
        form -- equal weight, the per-configuration block average).
    s : the ``(nao, nao)`` overlap metric.
    patterns : the window's chemical shell names, in slot order.
    layer : the shell's realized bands mapping -- the PRIOR
        realization (the case's ``prior_bands``: the widened,
        certificate-declared band lists k>1 needs), NOT necessarily
        the case's incumbent ``shell_bands``.
    radials : the declared depth, shell -> k (the tag's k, uniform
        over the window's shells by the tag grammar).
    ncol : the fold boundary (core orbital pairs).

    Returns
    -------
    ``(columns, report)`` : the ``(nao, sum k*(2l+1))`` S-orthonormal
    window columns (radial-major) and the per-shell report -- bands,
    k, the normalized spectrum, the DISCARDED fraction (the budget
    disclosure: the spectrum mass beyond k), and the capture of each
    source configuration's significant columns on the distilled span.

    Raises
    ------
    ValueError
        an undeclared shell, an absent band, k exceeding the chain
        width, or an assembled window that is not S-orthonormal
        (same-l shell blocks overlap -- the incumbent path's same
        undeclared case).
    """
    from atomterm.windows import shell_realization

    undeclared = [p for p in patterns if p not in radials]
    if undeclared:
        raise ValueError(
            f"the prior window is undeclared for {undeclared} -- every "
            "window shell carries the tag's k")
    n = mol.nao_nr()
    blocks: list[np.ndarray] = []
    report: dict[str, dict] = {}
    for shell in patterns:
        k = radials[shell]
        bands = shell_realization(shell, layer)
        chain, s_r = member_chains(mol, shell, bands, s)
        if k > len(chain[0]):
            raise ValueError(
                f"shell {shell!r}: the tag's k={k} exceeds its chain "
                f"width {len(chain[0])} (the realization "
                f"{list(bands)} holds fewer radial directions -- widen "
                "the declared prior bands or lower k)")
        avg = sum(valence_block(mf, chain, ncol) for mf in mfs) / len(mfs)
        v, lam = metric_radials(avg, s_r, k)
        if lam[0] <= 0:
            raise ValueError(
                f"shell {shell!r}: the averaged signal block is "
                f"numerically zero -- no occupied radial (a virtual "
                "shell has none)")
        if lam[k - 1] <= _TAIL * lam[0]:
            raise ValueError(
                f"shell {shell!r}: the declared k={k} exceeds the radial "
                f"signal -- the k-th eigenvalue reads {lam[k - 1]:.2e} "
                f"against the leading {lam[0]:.6f} (a single-configuration "
                "block is rank-1; only the LADDER average carries the "
                "second radial -- lower k or average over the ladder)")
        discarded = (float(lam[k:].sum() / lam.sum())
                     if len(lam) > k else 0.0)
        block = np.zeros((n, k * len(chain)))
        col = 0
        for j in range(k):
            for m in range(len(chain)):
                block[chain[m], col] = v[:, j]
                col += 1
        blocks.append(block)
        report[shell] = {
            "bands": list(bands), "k": k,
            "spectrum": (lam / lam[0]).round(5).tolist(),
            "discarded": round(discarded, 6),
            "capture": _capture(block, mfs, s, chain),
        }
    columns = np.hstack(blocks)
    err = np.abs(columns.T @ s @ columns
                 - np.eye(columns.shape[1])).max()
    if err > 1e-9:
        raise ValueError(
            f"the assembled prior window is not S-orthonormal "
            f"(max deviation {err:.2e}) -- same-l shell blocks overlap; "
            "the member-preserving cross-shell orthogonalization is a "
            "future seam")
    return columns, report


def _capture(block: np.ndarray, mfs: list, s: np.ndarray, chain
             ) -> list[list[float]]:
    """Each source configuration's significant columns, captured on
    the distilled span (the honest quality readout, per
    configuration -- the ladder's own columns, not a foreign
    reference's)."""
    from atomterm.windows import span_projector

    rows = sorted(r for m in range(len(chain)) for r in chain[m])
    p = span_projector(s, rows)
    out = []
    for mf in mfs:
        c = mf.mo_coeff[:, np.asarray(mf.mo_occ) > 0]
        norm = np.einsum("ai,ab,bi->i", c, s, c).real
        w = np.einsum("ai,ab,bi->i", c, p, c).real / norm
        sig = np.flatnonzero(w >= 0.05)
        cr = block.T @ (s @ c[:, sig])
        cap = np.einsum("ij,ij->j", cr, cr).real / norm[sig]
        out.append(cap.round(6).tolist())
    return out
