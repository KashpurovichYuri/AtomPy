"""Calculation drivers: ROHF, CASCI, CASSCF, SA-CASSCF with .chkfile caching.

Every stage writes a PySCF checkpoint file (MOs, densities) so re-runs are
fast: if the ROHF chkfile exists it is loaded instead of recomputed. Analysis
artifacts (energies, CI vectors) are stored under the ``atomterm/`` key
prefix of the same chkfile.

Active-space convention: the active orbitals are the AO functions of the
target valence shell (e.g. 2p), placed at ``[ncore : ncore+ncas]`` of the full
MO matrix. Same-shell AOs are orthonormal and exactly rotation-closed, which
is what makes the L-quantum-number analysis exact.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pyscf import lib, mcscf, scf

__all__ = [
    "ReferenceSpec",
    "run_rohf",
    "run_frac_occ_scf",
    "select_shell_ao",
    "select_valence_ao",
    "primitive_chain",
    "select_primitives",
    "occupied_core_idx",
    "shell_weights",
    "purified_s_core",
    "purified_core_orbitals",
    "complete_orthonormal",
    "build_full_mo",
    "lm_sectors",
    "radial_basis",
    "radial_mask",
    "run_casci",
    "constrained_casscf",
    "core_virtual_only_mask",
    "no_core_window_exchange_mask",
    "run_casscf",
    "run_sa_casscf",
    "run_casscf_core_relax",
    "run_casscf_virtual_mix",
    "run_casscf_pyscf_start",
    "avas_window",
    "span_overlap",
    "CASResult",
]


def _ao_span_projector(s_mat: np.ndarray, rows: list[int]) -> np.ndarray:
    """S-metric projector onto an AO row span: ``P = S[:, o] S_oo^{-1} S[o, :]``.

    The weight of an S-orthonormal MO column ``c`` on the span is
    ``c^dag P c`` -- bounded by 1 (a true projector weight), unlike the
    bare quadratic form ``c_o^dag S_oo c_o`` which is unbounded for
    diffuse bases (fact 15c).
    """
    s_oo = s_mat[np.ix_(rows, rows)]
    return s_mat[:, rows] @ np.linalg.inv(s_oo) @ s_mat[rows, :]


def _lowdin(c: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Symmetric Loewdin orthonormalization of the columns of ``c`` in the
    S metric. Span-preserving; the eigenvector gauge within a degenerate
    set is an arbitrary rotation -- fine unless column identity matters
    (:func:`radial_basis` fixes the gauge by Gram-Schmidt instead)."""
    ev, u = np.linalg.eigh(c.conj().T @ s @ c)
    return c @ (u / np.sqrt(ev)) @ u.T


def occupied_core_idx(mf, ncore: int) -> list[int]:
    """Indices of the ``ncore`` doubly-occupied MOs -- the honest core
    block for every build_full_mo/radial_basis call.

    Canonical MO ORDER is meaningless on near-degenerate uncontracted
    bases (U: the first canonical MOs are virtual junk at eps ~ 0.02,
    while the closed block sits at 22/61/88/120/143); ``range(ncore)``
    then silently folds five garbage columns (measured: max|dE| = 14.5
    Eh). On contracted bases the two coincide -- this makes the
    coincidence a checked invariant instead of an assumption.
    """
    closed = np.flatnonzero(np.asarray(mf.mo_occ) == 2.0)
    if len(closed) < ncore:
        raise ValueError(
            f"asked for {ncore} core orbitals; ROHF closed block has "
            f"{len(closed)}")
    if not np.array_equal(closed[:ncore], np.arange(ncore)):
        import warnings

        warnings.warn(
            f"core block is NOT the first {ncore} canonical MOs "
            f"(closed at {closed[:ncore].tolist()}); folding the "
            f"occupied block instead", stacklevel=2)
    return closed[:ncore].tolist()


def _span_rows(mol, pattern: str) -> list[int]:
    """AO rows of one occupation-span pattern: ``"l0".."l9"`` selects a
    WHOLE angular-momentum subspace (the honest vocabulary for uncontracted
    bases, where labels are shell counters -- fact 1); anything else is a
    label prefix (:func:`select_shell_ao`)."""
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


def _forced_rohf_occ(mol, occ_closed: list[str], occ_open: list[str]):
    """Build the character-forced ROHF ``get_occ`` hook (see :func:`run_rohf`).

    Occupations by S-projector WEIGHT on two AO spans, never energy order
    (the counts and spans: :func:`run_rohf`). The open block is ranked by
    ``2 w_open + w_closed``: open-shell character carries a strict premium
    (with whole-l spans a 1:1 sum lets closed-character leftovers tie the
    open MOs -- measured on U, where DIIS then closed garbage), while the
    closed term still lets early-cycle open-shell character hiding in
    core-like MOs be found (the Sc/Ti leakage case).
    """
    s_mat = mol.intor_symmetric("int1e_ovlp")
    p_closed = _ao_span_projector(
        s_mat, sum((_span_rows(mol, p) for p in occ_closed), []))
    p_open = _ao_span_projector(
        s_mat, sum((_span_rows(mol, p) for p in occ_open), []))
    n_closed = (mol.nelectron - mol.spin) // 2
    n_open = mol.spin

    def _occ(mo_energy, mo_coeff=None):
        if mo_coeff is None:
            raise ValueError("character-forced occupations need the MOs")
        occ = np.zeros(len(mo_energy))
        w_closed = np.einsum("ai,ab,bi->i", mo_coeff.conj(), p_closed,
                             mo_coeff).real
        closed = np.argsort(w_closed)[::-1][:n_closed]
        occ[closed] = 2.0
        w_open = np.einsum("ai,ab,bi->i", mo_coeff.conj(), p_open,
                           mo_coeff).real
        rest = np.flatnonzero(occ == 0.0)
        open_mos = sorted(rest,
                          key=lambda i: -(2 * w_open[i] + w_closed[i]))[:n_open]
        occ[open_mos] = 1.0
        return occ

    return _occ


@dataclass(frozen=True)
class ReferenceSpec:
    """How the ROHF reference is produced when plain aufbau filling fails
    (fact 21a). One object bundling character-forced occupations and the
    basin gate, threaded through run_rohf / run_ladder / run_atom.

    occ_closed, occ_open : AO span patterns (label prefixes as
        :func:`select_shell_ao`, or ``"l0".."l9"`` for a WHOLE angular
        momentum subspace -- the honest vocabulary on uncontracted bases
        where labels are shell counters) forcing the occupations by
        orbital CHARACTER instead of energy order -- the integer analog
        of :func:`run_frac_occ_scf`'s ``open_idx`` selection. MANDATORY
        on atoms whose default ROHF is an aufbau runaway (Ti; U needs
        whole-l spans plus Newton).
    occ_e_ref : reference energy (Eh) of the target SCF basin; gates the
        multi-basin chaos by guess cycling (see :func:`run_rohf`).
    scf_method : "diis" (default) or "newton" (second-order SOSCF). On
        heavy uncontracted bases with near-degenerate f/d/p valence
        shells, plain DIIS with forced occupations can fail to converge
        from ANY start (U: E > +1e3, garbage blocks), while Newton lands
        the target basin deterministically from every guess tested.
    """

    occ_closed: tuple[str, ...] | None = None
    occ_open: tuple[str, ...] | None = None
    occ_e_ref: float | None = None
    scf_method: str = "diis"

    def __post_init__(self):
        if (self.occ_closed is None) != (self.occ_open is None):
            raise ValueError("occ_closed and occ_open must be given together")
        if self.occ_e_ref is not None and self.occ_closed is None:
            raise ValueError(
                "occ_e_ref needs the occ_closed/occ_open patterns")
        if self.scf_method not in ("diis", "newton"):
            raise ValueError(f"unknown scf_method: {self.scf_method!r}")


def run_rohf(cfg, chkfile: str, ref: ReferenceSpec | None = None):
    """ROHF for the configured atom, cached in ``chkfile``.

    If the chkfile exists, the converged MOs are loaded and no SCF is run.
    Returns ``(mol, mf)``.

    ref : a :class:`ReferenceSpec` (character-forced occupations +
        optional basin gate). Occupations are assigned each cycle by
        S-projector weight on the ``occ_closed``/``occ_open`` AO spans,
        never by energy order -- mandatory on aufbau-runaway atoms
        (Ti: +123.07 mEh runaway, fact 21a / README §5.5). Even forced,
        a degenerate open shell is multi-basin (Ti: +11.2 mEh
        metastable), so with ``occ_e_ref`` set a FRESH run is retried
        over cycled init guesses until |E - occ_e_ref| <= 1e-4, else
        RuntimeError; cached chkfiles are trusted without the gate. The
        occupation hook is part of the object's identity: set BEFORE
        the cache branch (ticket-10 lesson).
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)

    def fresh_mf(guess: str | None = None):
        m = scf.ROHF(mol)
        if ref is not None and ref.occ_closed is not None:
            m.get_occ = _forced_rohf_occ(mol, list(ref.occ_closed),
                                         list(ref.occ_open))
        m.chkfile = chkfile
        m.conv_tol = 1e-12  # tight core density: keeps symmetry clean
        m.verbose = cfg.verbose
        if guess is not None:
            m.init_guess = guess
        if ref is not None and ref.scf_method == "newton":
            m = m.newton()
        return m

    try:
        _, rec = scf.chkfile.load_scf(chkfile)
        mf = fresh_mf()
        mf.__dict__.update(rec)
        mf.converged = True
        mf.verbose = cfg.verbose
        return mol, mf
    except (OSError, KeyError, ValueError):
        pass
    mf = fresh_mf()
    mf.run()
    if ref is not None and ref.occ_e_ref is not None:
        for guess in ("atom", "huckel"):
            if mf.converged and abs(mf.e_tot - ref.occ_e_ref) <= 1e-4:
                break
            mf = fresh_mf(guess)
            mf.run()
        else:
            if not (mf.converged and abs(mf.e_tot - ref.occ_e_ref) <= 1e-4):
                raise RuntimeError(
                    f"character-forced ROHF did not reach the reference "
                    f"basin E = {ref.occ_e_ref:.8f} (got {mf.e_tot:.8f}, "
                    f"converged={mf.converged}); the basin census or the "
                    f"guess list needs revisiting")
    return mol, mf


def select_shell_ao(mol, pattern: str) -> list[int]:
    """AO indices of a named valence shell, e.g. pattern="2p" -> [3, 4, 5].

    Plain prefix matching on PySCF's AO labels (``search_ao_label`` patterns
    are regex-flavored and over-match: "2p*" also hits "2s").
    """
    labels = [line.split()[-1] for line in mol.ao_labels()]
    return [i for i, lbl in enumerate(labels) if lbl.startswith(pattern)]


def select_valence_ao(mol, patterns: list[str]) -> list[int]:
    """AO indices of several valence shells, concatenated in pattern order,
    e.g. ``["3d", "4s"]`` -> the five 3d AOs followed by the 4s AO.

    The concatenation defines the active space; it is rotation-closed when
    each pattern contributes a complete m-set of one (n, l) shell.
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


def primitive_chain(mol, l: int) -> list[tuple[float, list[int]]]:
    """The radial primitives of angular momentum ``l`` as ``(exponent,
    ao_rows)`` pairs, sorted by exponent ASCENDING (diffuse end first) --
    the honest radial coordinate of an uncontracted basis, where PySCF's
    AO labels number shells from the TIGHT end starting at the valence
    name ("5f" = the exp-36.8 primitive on U; fact-1 amendment).

    Walks the basis exactly as :func:`lm_sectors` does (``bas_angular``
    checked, ``nctr`` unpacked): a primitive's m-set occupies the rows
    ``loc[ib] + c*(2l+1) + m`` of its shell. A contraction slot is a
    RAW primitive iff its ``bas_ctr_coeff`` column carries exactly one
    exponent; true contractions (several exponents combined into one
    function) have no single exponent and RAISE -- for them use
    :func:`select_valence_ao` (a contracted shell is one pre-converged
    radial) or MO/NO-based window selection.
    """
    loc = mol.ao_loc_nr()
    width = 2 * l + 1
    chain: list[tuple[float, list[int]]] = []
    for ib in range(mol.nbas):
        if abs(mol.bas_angular(ib)) != l:
            continue
        exps = np.asarray(mol.bas_exp(ib))
        coeff = mol.bas_ctr_coeff(ib)  # (nprim, nctr)
        for c in range(coeff.shape[1]):
            nz = np.flatnonzero(np.abs(coeff[:, c]) > 1e-12)
            if len(nz) != 1:
                raise ValueError(
                    f"shell {ib} (l={l}) contraction slot {c} combines "
                    f"{len(nz)} exponents -- exponent-ranked selection is "
                    "undefined for contractions; use select_valence_ao "
                    "(a contracted shell is one pre-converged radial) or "
                    "MO/NO-based window selection"
                )
            rows = [loc[ib] + c * width + m for m in range(width)]
            chain.append((float(exps[nz[0]]), rows))
    chain.sort(key=lambda t: t[0])
    return chain


def select_primitives(mol, l: int, count: int = None,
                      exponent_band: tuple[float, float] = None) -> list[int]:
    """AO rows of radial primitives of angular momentum ``l``, selected
    by EXPONENT on the :func:`primitive_chain` -- m-complete by
    construction (each primitive contributes its whole 2l+1 m-set), so
    the window is exactly rotation-closed.

    Exactly one of ``count`` / ``exponent_band``:
    * ``exponent_band=(lo, hi)``: every primitive with ``lo <= exp <=
      hi``. The scan protocol: bands are CONTIGUOUS sets built from the
      measured open-MO primitive weights, never hand-picked intervals,
      and conclusions require neighboring-band stability (admit one
      more / one fewer primitive on each side -> same ordering).
    * ``count``: the ``count`` most diffuse primitives (convenience
      shorthand for the diffuse end of the chain).
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


def avas_window(mf, labels: list[str], ncore: int, threshold: float = None):
    """Window chosen by AVAS (pyscf.mcscf.avas, arXiv:1701.07862): an
    independent ecosystem window selector for cross-validating our raw
    shell AOs (verdicts: fact 19, README §7.3).

    AVAS projects the target AO span (``labels``, e.g. ``["Sc 3d",
    "Sc 4s"]``) out of the converged MOs and keeps eigenvectors of the
    projected overlap above ``threshold`` (its default 0.2). Because
    the construction PROJECTS onto the target span, a full
    shell-dimension AVAS window IS the shell span exactly -- closed by
    construction; a threshold that drops members of an m-set gives a
    rotation-DEFECTIVE span (the negative-control case for the closure
    meters).

    Returns ``(mo_full, ncas, nelecas)`` with the AVAS active span at
    ``[ncore : ncore+ncas]`` (ready for run_casci; electron bookkeeping
    stays ours -- pass the sector you want).
    """
    from pyscf.mcscf import avas

    kw = {} if threshold is None else {"threshold": threshold}
    ncas, nelecas, mo = avas.avas(mf, labels, ncore=ncore, **kw)
    nelecas = np.atleast_1d(np.asarray(nelecas)).ravel()
    return mo, int(ncas), nelecas


def span_overlap(c_a: np.ndarray, c_b: np.ndarray, s: np.ndarray) -> float:
    """Smallest principal cosine between two spans in the S metric
    (1.0 = identical spans, 0.0 = orthogonal direction present).

    The identity-change gauge of the ladder decomposition: sv_min of
    ``C_a^dag S C_b`` is 1 exactly when the converged window kept the
    shell span, and drops first when the window identity changes.
    """
    sv = np.linalg.svd(c_a.conj().T @ s @ c_b, compute_uv=False)
    return float(sv.min())


def run_frac_occ_scf(cfg, chkfile: str, nelec_open: int, n_open: int,
                     open_idx: list[int] | None = None):
    """Average-of-configuration (fractional-occupation) HF, chkfile-cached.

    The open shell is occupied fractionally (``nelec_open / n_open`` per
    orbital, e.g. 2/3 per 2p AO for the carbon 2p^2 configuration), which
    makes the converged Fock exactly spherical -- the scalar-core
    reference that density averaging reproduces, but a DIFFERENT physical
    reference from ROHF (smeared vs localized valence exchange; the
    measured energy relation: README §3.2, fact 15f).

    Control mechanism: ``get_occ`` is PySCF's occupation hook, called
    every SCF cycle; the density and Fock are built from whatever it
    returns, so the SCF converges self-consistently in the ensemble
    density. The returned ``e_tot`` is the ensemble energy of the
    averaged configuration -- use this object as a CORE PROVIDER, never
    as a term energy.

    ``open_idx`` selects the fractional block by CHARACTER (AO indices
    of the open shell): the ``n_open`` non-core orbitals with the largest
    S-projector weight on that shell get the fractional occupation.
    This matters for transition metals: energy-ordered filling is
    unstable where 3d/4p quasi-degenerate (on Sc the aufbau block runs
    away to {4s, 4p, ~2d} -- fact 15c). Without ``open_idx`` the block
    is the HOMO block (energy order; fine for C/Si, where it selects
    the same orbitals as the projection would).

    Closed-shell RHF formalism; spin in ``cfg`` is ignored (only its
    nelectron/spin PARITY is used). Returns ``(mol, mf)``.
    """
    from atomterm.system import build_mol

    # cfg only needs the right nelectron/spin parity: the fractional
    # occupations below define the density; the RHF formalism never sees
    # mol.spin.
    mol = build_mol(cfg)
    s_mat = None
    g_inv = None
    if open_idx is not None:
        assert len(open_idx) == n_open, (
            f"open_idx has {len(open_idx)} AOs but n_open={n_open}")
        s_mat = mol.intor_symmetric("int1e_ovlp")
        # S-projector metric on the open-shell AO span (bounded by 1 for
        # S-orthonormal MOs, unlike a bare quadratic form): _ao_span_projector
        g_inv = np.linalg.inv(s_mat[np.ix_(open_idx, open_idx)])

    def _frac_occ(mo_energy, mo_coeff=None):
        occ = np.zeros(mo_energy.size)
        n_closed = (mol.nelectron - nelec_open) // 2
        occ[:n_closed] = 2.0
        if open_idx is None:
            occ[n_closed : n_closed + n_open] = nelec_open / n_open
        else:
            # valence weight of each MO: projection onto the open shell's
            # AO span (rotation-invariant gauge for the fractional block)
            r = s_mat[open_idx, :] @ np.asarray(mo_coeff)
            w_v = np.einsum("ai,ab,bi->i", r.conj(), g_inv, r).real
            pick = np.argsort(w_v[n_closed:])[::-1][:n_open] + n_closed
            occ[pick] = nelec_open / n_open
        return occ

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.verbose = cfg.verbose
    # the occupation hook is part of the object's identity: set BEFORE
    # the cache branch (ticket-10 lesson)
    mf.get_occ = _frac_occ
    try:
        _, rec = scf.chkfile.load_scf(chkfile)
        mf.__dict__.update(rec)
        mf.converged = True
        return mol, mf
    except (OSError, KeyError, ValueError):
        pass

    mf.chkfile = chkfile
    mf.run()
    return mol, mf


def shell_weights(mol, c: np.ndarray, s: np.ndarray):
    """S-metric projection norms of columns onto angular (l) AO subspaces.

    ``w[i, k] = c[:, i]^dag S_l c[:, i]`` over ALL AOs of angular momentum
    l_k. S-orthonormal columns make every row of ``w`` sum to 1 -- valid
    because different-l subspaces are exactly orthogonal; grouping by
    (n, l) instead would NOT partition (2p and 3p AO subspaces overlap,
    weights double-count). Same AO-metric quantity as the
    ``contamination`` of :func:`purified_s_core`. l-purity is NECESSARY
    but not SUFFICIENT for rotational closure -- only the commutators
    measure closure itself (facts 8, 19).
    Returns ``(names, w)`` with ``names = ['s', 'p', 'd', ...]``.
    """
    import re

    labels = [line.split()[-1] for line in mol.ao_labels()]
    groups: dict[str, list[int]] = {}
    for i, lbl in enumerate(labels):
        m = re.match(r"\d+([spdfgh])", lbl)
        groups.setdefault(m.group(1) if m else lbl, []).append(i)

    names = list(groups)
    w = np.zeros((c.shape[1], len(names)))
    for k, name in enumerate(names):
        idx = groups[name]
        ss = s[np.ix_(idx, idx)]
        w[:, k] = np.einsum("ij,ij->j", c[idx, :], ss @ c[idx, :]).real
    return names, w


def purified_s_core(mol, mf, n_core: int):
    """Canonical core MOs projected onto the s-AO subspace (l = 0 exactly).

    The converged "s" MOs of an open-shell atom carry a tiny p/d tail
    (measured on C/cc-pVDZ: 1s 1.9e-7, 2s 6.3e-6 weight deficit) which
    breaks two preconditions at once -- orthonormality (<MO|2p_AO>
    ~ 1e-3, so second quantization is ill-defined on such columns) and
    rotational closure ([H, Lz] fails at ~1e-4 instead of machine
    precision). Projection + a Loewdin inside the s subspace repairs
    both at second-order cost (~1e-10 Eh) -- the orbital analog of the
    density-space spherical average. Full causal chain: fact 9.

    Returns ``(C, contamination)``: (nao, n_core) orthonormal columns
    and the dropped fraction ``1 - <c_s|S|c>`` per canonical MO (the SCF
    symmetry-breaking diagnostic; MO weights live in the AO metric,
    c^dag S c = 1, not the Euclidean norm).
    """
    import re

    labels = [line.split()[-1] for line in mol.ao_labels()]
    s_idx = [i for i, lbl in enumerate(labels) if re.fullmatch(r"\d+s", lbl)]

    c_can = mf.mo_coeff[:, occupied_core_idx(mf, n_core)]
    c_s = c_can[s_idx, :]                                # drop the p/d tail
    # MO normalization lives in the AO metric (c^dag S c = 1), so the
    # s-subspace weight is c_s^dag S_ss c_s -- not the Euclidean norm.
    s_ov = mf.get_ovlp()[np.ix_(s_idx, s_idx)]
    weight = np.einsum("ij,ij->j", c_s, s_ov @ c_s)
    contamination = np.clip(1.0 - weight, 0.0, None)     # fp noise can dip < 0
    c_s = c_s / np.sqrt(weight)                          # renormalize the rest

    c_s = _lowdin(c_s, s_ov)

    c = np.zeros((mol.nao_nr(), n_core), dtype=c_s.dtype)
    c[s_idx, :] = c_s
    return c, contamination


def purified_core_orbitals(mol, mf, n_core: int):
    """Canonical core MOs replaced by an exactly rotation-closed core span.

    Generalizes :func:`purified_s_core` to p/d cores (Si [Ne]+3s, Sc
    [Ar]). The canonical core MOs of an open-shell atom carry foreign-l
    tails and m-ASYMMETRIC same-l content; both break the L-closure of
    the core SPAN, and every window orthogonalized against that span
    inherits the defect (the fact-17 input-window baseline, ~1e-3 Si /
    ~1e-2 Sc). Closure theorem (E (x) members is L-invariant for any
    radial subspace E): fact 20, README §5.3. Construction per l:

    1. member-resolved radial blocks ``R = [r_{i,m}]`` (each column =
       MO i's coefficients on member m's radial chain);
    2. E_l = top-k_l span of the S-metric column Gram ``R^T S_rad R``,
       S-orthonormalized -- eigenvalues = S-weight x member multiplicity
       (~2l+1 per occupied radial), so the beyond-k spectrum is the tail
       content, asserted as such;
    3. every radial of E_l replicated on ALL member sectors of l
       (k_l x (2l+1) columns -- m-complete by construction).

    One global Loewdin then orthonormalizes the set (S = 0 between
    different l keeps l-purity; within-l mixing is span-preserving).
    TRAPS (fact 20): the S-generalized eigenproblem of ``R R^T`` returns
    S^{-1}-warped radials -- the column Gram is the correct object; and
    canonical core MOs are NOT shell-confined (C 2s: 66% of S-weight off
    its shell) -- per-shell projection is wrong. For l = 0 the
    construction reduces exactly to :func:`purified_s_core`.

    Returns ``(C, contamination)``: (nao, n_core) S-orthonormal columns
    and per-MO ``1 - own-l weight`` (the foreign-l SCF tail diagnostic).
    """
    s = mf.get_ovlp()
    c_can = np.asarray(
        mf.mo_coeff[:, occupied_core_idx(mf, n_core)])
    sectors = lm_sectors(mol)

    # dominant l per core MO (lumped-l weights partition unity exactly;
    # the foreign-l weight is the tiny numerical tail, fact 9)
    l_names, w_l = shell_weights(mol, c_can, s)
    dom_l = [l_names[i] for i in w_l.argmax(axis=1)]
    for i, wmax in enumerate(w_l.max(axis=1)):
        assert wmax > 0.95, (i, wmax)  # a mixed-l core MO = wrong core

    c_pur = np.zeros((mol.nao_nr(), n_core), dtype=c_can.dtype)
    contamination = np.zeros(n_core)
    col = 0
    for k_col, l_name in enumerate(l_names):
        l_val = "spdfgh".index(l_name)
        mos = [i for i, d in enumerate(dom_l) if d == l_name]
        if not mos:
            continue
        members = sorted(m for (ll, m) in sectors if ll == l_val)
        k = len(mos) // len(members)   # radial copies (n_mos / m-multiplicity)
        assert k * len(members) == len(mos), (
            f"core MOs of l={l_name} do not form complete m-sets: {len(mos)}")
        chain = sectors[(l_val, members[0])]
        s_rad = s[np.ix_(chain, chain)]

        # foreign-l contamination diagnostic
        for i in mos:
            contamination[i] = 1.0 - w_l[i, k_col]

        # dominant radial subspace E_l: top-k of the S-metric column Gram
        # of the member-resolved radial blocks
        rr = np.hstack([c_can[sectors[(l_val, m)], :][:, mos]
                        for m in members])
        ev, q = np.linalg.eigh(rr.conj().T @ s_rad @ rr)
        order = np.argsort(ev)[::-1]
        top = order[:k]
        assert ev[top[-1]] > 1e-8, (l_name, ev)
        if len(ev) > k:  # discarded directions are tails, not radials
            assert ev[order[k]] < 1e-2 * ev[top[-1]], (l_name, ev)
        v_all = rr @ q[:, top] @ np.diag(ev[top] ** -0.5)  # S-orthonormal
        for j in range(k):
            for m in members:  # shared radial on EVERY member (m-complete)
                c_pur[sectors[(l_val, m)], col] = v_all[:, j]
                col += 1
    assert col == n_core

    # one global Loewdin (S = 0 between different l -> l-purity survives;
    # within-l mixing is span-preserving)
    return _lowdin(c_pur, s), contamination


def complete_orthonormal(v: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Append an orthonormal-in-metric-``s`` complement to the columns of ``v``."""
    ev, u = np.linalg.eigh(s)
    s_half = (u * np.sqrt(ev)) @ u.T
    s_ihalf = (u / np.sqrt(ev)) @ u.T
    vx = s_half @ v
    resid = np.eye(v.shape[0]) - vx @ vx.conj().T
    p, sv, _ = np.linalg.svd(resid)
    c = s_ihalf @ p[:, sv > 1e-8]
    return np.hstack([v, c])


def build_full_mo(mf, core_idx: list[int], act_idx: list[int], c_core: np.ndarray = None) -> np.ndarray:
    """Full orthonormal MO matrix with chosen core and active columns placed
    so that CASCI/CASSCF's window ``[ncore : ncore+ncas]`` selects the active
    space.

    Core columns are canonical SCF MOs by default; ``c_core`` supplies
    alternative columns (e.g. fractional-occupation core). Active columns
    are the RAW AO functions of the target shell -- exactly rotation-closed,
    essential for the L-quantum-number analysis.

    Construction is block-safe in the AO metric S (fact 14): (1) the core
    block is Loewdin-orthonormalized within its own span (span- and
    density-preserving); (2) the active columns are projected against the
    core span and renormalized (raw valence AOs overlap the canonical
    core's same-l members: <2p_MO|3p_AO> can reach ~0.3) -- a same-l,
    m-diagonal correction, so window l-purity and m-set structure
    survive; (3) the complement is appended through an exact projector --
    the naive `I - P` filter is exact only for S-orthonormal v and
    silently returns extra columns otherwise (observed: (nao, nao+6)
    on Si).
    """
    v = _core_window_block(mf, core_idx, act_idx, c_core)
    if v.shape[1] == mf.mo_coeff.shape[0]:
        return v
    return complete_orthonormal(v, mf.get_ovlp())


def _core_window_block(mf, core_idx: list[int], act_idx: list[int],
                       c_core: np.ndarray = None) -> np.ndarray:
    """Steps (1)+(2) of :func:`build_full_mo`: the S-orthonormal
    core+window block that both it and :func:`radial_basis` start from."""
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    core = mf.mo_coeff[:, core_idx] if c_core is None else c_core
    act = np.eye(n)[:, act_idx]

    core = _lowdin(core, s)

    p_core = core @ core.conj().T @ s
    act = act - p_core @ act
    act = _lowdin(act, s)

    return np.hstack([core, act])


def lm_sectors(mol) -> dict[tuple[int, int], list[int]]:
    """AO rows grouped by (l, member): the EXACTLY S-orthogonal sectors of
    a one-center spherical basis -- the AO space is a direct sum of radial
    chains, one per (l, member); same-sector functions of different n
    OVERLAP, different-sector functions never do (angular orthogonality,
    fact 16e). This decomposition is what makes a per-sector treatment
    (:func:`radial_basis`) globally S-orthonormal by construction.
    """
    sectors: dict[tuple[int, int], list[int]] = {}
    loc = mol.ao_loc_nr()
    for ib in range(mol.nbas):
        l = abs(mol.bas_angular(ib))
        nctr = mol.bas_nctr(ib)
        for c in range(nctr):
            for member in range(2 * l + 1):
                # fused shells (nctr > 1) are contraction-outer, m-inner
                sectors.setdefault((l, member), []).append(
                    loc[ib] + c * (2 * l + 1) + member)
    return sectors


def radial_basis(mf, core_idx: list[int], act_idx: list[int],
                 c_core: np.ndarray = None):
    """Full MO matrix plus per-slot (l, member) labels -- the substrate for
    :func:`radial_mask`. ``c_core`` overrides the core block; pass
    :func:`purified_core_orbitals` output for the exactly-closed
    substrate (the fact-17 fix).

    Everything non-core is built in a MEMBER-PRESERVING gauge (sequential
    S-metric Gram-Schmidt, not a symmetric Loewdin): the Loewdin gauge
    within a degenerate member set is an arbitrary rotation --
    span-preserving, but it reassigns which column is which (l, m), and
    here column identity IS the product (fact 16e, README §3.3). Window
    slots carry their AO's sector label; STRONG sector residuals
    (>= 8e-2 of S-weight; core tails <= 3e-6) join with their label; the
    near-linearly-dependent needles (Sc 6p: 1.7e-4) pad through
    :func:`complete_orthonormal` with label ``None`` -- never rotation
    partners.

    Returns ``(mo_full, labels)`` with ``labels[i] = (l, member) | None``.
    """
    mol = mf.mol
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    sec = lm_sectors(mol)
    row_sector = {row: key for key, rows in sec.items() for row in rows}

    core = mf.mo_coeff[:, core_idx] if c_core is None else c_core
    core = _lowdin(core, s)

    def sgs(cols, base):
        """S-metric Gram-Schmidt (two refinement passes) of the columns in
        ``cols`` on top of the S-orthonormal block(s) in ``base``;
        preserves column identity (unlike a symmetric Loewdin)."""
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

    # window in member order: raw AO minus the core span, sequentially
    # orthogonalized -- column j IS the (l, m) of act_idx[j] by construction
    p_core = core @ core.conj().T @ s
    window = sgs([np.eye(n)[:, ao] - p_core @ np.eye(n)[:, ao]
                  for ao in act_idx], [core])
    v = np.hstack([core] + window)
    labels: list = [None] * len(core_idx) + [row_sector[ao] for ao in act_idx]

    # strong sector residuals, same member-preserving treatment. Strength
    # cutoff 8e-2 / needle class 1.7e-4: docstring + fact 16e (the 3s of
    # cc-pVDZ C survives at 0.085; the Sc 6p needle does not)
    p_v = v @ v.conj().T @ s
    strong, strong_labels = [], []
    for key, rows in sorted(sec.items()):
        block = np.eye(n)[:, rows]
        r = block - p_v @ block
        gr = r.conj().T @ s @ r
        ev, u = np.linalg.eigh(gr)
        keep = ev > 1e-4
        if keep.any():
            strong += [r @ (u[:, j] / np.sqrt(ev[j])) for j in np.where(keep)[0]]
            strong_labels += [key] * int(keep.sum())
    if strong:
        v = np.hstack([v] + sgs(strong, [v]))
        labels += strong_labels

    mo = v if v.shape[1] == n else complete_orthonormal(v, s)
    labels += [None] * (mo.shape[1] - len(labels))

    # gauge assertion: every labeled column really is dominated by its own
    # sector (catches any gauge surprise at the purity level)
    for i, lab in enumerate(labels):
        if lab is None:
            continue
        rows = sec[lab]
        c = mo[rows, i]
        w = float(np.real(c.conj().T @ s[np.ix_(rows, rows)] @ c))
        assert w > 0.9, (i, lab, w)
    assert mo.shape[1] == n
    return mo, labels


def radial_mask(labels: list, ncore: int, ncas: int):
    """``uniq_var_indices`` mask allowing ONLY same-(l, member) rotations
    into a window slot: pure RADIAL mixing (3d_m may borrow 4d_m/5d_m; the
    4s any s), l-content and m-membership of every slot untouchable, and
    no core<->window exchange at all (core slots carry no label).

    Physics caveat (measured, fact 16): exact rotational closure of the
    window additionally requires the radial rotation applied IDENTICALLY
    to all m partners; independent per-member radial mixing leaves
    [H, Lz] nonzero at the member-difference level. SA over the term's
    M_L set is the symmetry-restoring companion rung.
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


@dataclass
class CASResult:
    """Uniform result container for CASCI / CASSCF / SA-CASSCF."""

    kind: str
    energies: np.ndarray  # total energies of the roots
    ci: list = field(default_factory=list)
    mo_coeff: np.ndarray = None
    mc: object = None  # live pyscf object
    core_energy: float = 0.0  # total minus active-space energy


def _persist(chkfile: str, name: str, energies, ci) -> None:
    if chkfile is None:
        return
    lib.chkfile.save(chkfile, f"atomterm/{name}/energies", np.asarray(energies))
    lib.chkfile.save(chkfile, f"atomterm/{name}/ci", np.asarray(ci, dtype=object))


def _core_energy(mc) -> float:
    """total(e_tot) minus active-space energy (weight-aware for SA runs)."""
    e_cas = np.atleast_1d(np.asarray(mc.e_cas)).ravel().astype(float)
    e_tot = float(np.atleast_1d(np.asarray(mc.e_tot)).ravel()[0])
    weights = getattr(mc.fcisolver, "weights", None)
    if weights is not None and len(weights) == len(e_cas):
        avg_active = float(np.dot(weights, e_cas))
    else:
        avg_active = float(e_cas[0])
    return e_tot - avg_active


def _root_energies(mc) -> np.ndarray:
    """Total energies of all roots (active-space e_cas + common core part)."""
    e_cas = np.atleast_1d(np.asarray(mc.e_cas)).ravel().astype(float)
    return e_cas + _core_energy(mc)


def _run(mc, nroots: int, chkfile: str, name: str, verbose: int) -> CASResult:
    mc.fcisolver.nroots = nroots
    mc.verbose = verbose
    mc.chkfile = chkfile
    mc.kernel()
    energies = _root_energies(mc)  # nroots=1 reproduces float(mc.e_tot)
    ci = [mc.ci] if nroots == 1 else list(mc.ci)
    _persist(chkfile, name, energies, ci)
    return CASResult(name, energies, ci, mc.mo_coeff, mc, _core_energy(mc))


def run_casci(mf, ncas, nelecas, mo_full, chkfile, nroots=1, verbose=0,
              name="casci", fci_opts: dict = None) -> CASResult:
    """CASCI in the given window. ``name`` is the chkfile store key -- keep
    it distinct per variation so differently-folded runs never collide.
    ``fci_opts``: attributes set on the fcisolver before the kernel (e.g.
    ``{"pspace_size": ndet}`` forces pyscf's EXACT full-diagonalization
    branch -- Davidson straggles inside degenerate clusters: measured on
    U 2f2d/Ms=2, 48 roots: up to 38 mEh within 13-17-fold multiplets
    while both exact spectra agree to 1e-7)."""
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
    for k, v in (fci_opts or {}).items():
        setattr(mc.fcisolver, k, v)
    return _run(mc, nroots, chkfile, name, verbose)


# ---------------------------------------------------------------------------
# Orbital-rotation masks (seam of :func:`constrained_casscf`).
#
# PySCF's CASSCF optimizer enumerates its rotation variables through the
# ``uniq_var_indices(nmo, ncore, ncas, frozen)`` hook, which returns a
# boolean (nmo, nmo) pair mask: the generator entry (i, j) is optimized
# iff mask[i, j]. Rebinding that method constrains the optimization WITHOUT
# touching the optimizer -- each mask below is one declarative constraint.
# ---------------------------------------------------------------------------


def core_virtual_only_mask(self, nmo, ncore, ncas, frozen):
    """ONLY inactive<->virtual rotations allowed (the relax rung)."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, :ncore] = True
    return mask


def no_core_window_exchange_mask(self, nmo, ncore, ncas, frozen):
    """Everything EXCEPT active<->inactive rotations: no orbital ever
    trades places between the core and the active window (the virtmix
    rung)."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, : ncore + ncas] = True
    return mask


def constrained_casscf(mf, ncas, nelecas, chkfile, *, mo_full=None, mask=None,
                       weights=None, verbose=0, name="casscf",
                       mc_opts: dict = None) -> CASResult:
    """CASSCF with declaratively constrained orbital rotations -- the one
    deep driver every orbital-optimization rung goes through.

    mo_full : starting MOs (e.g. raw shell columns via :func:`build_full_mo`);
        ``None`` = PySCF's default start (canonical MOs, canonical window) --
        the start-variant axis.
    mask : an ``uniq_var_indices``-shaped callable (see the masks above, or
        build one per window from AO labels); ``None`` = free optimization.
    weights : state-average over these weights. PySCF's SA kernel returns
        only the ENSEMBLE energy, so per-root energies come from a CASCI
        pass in the SA-converged window.
    name : chkfile store key under ``atomterm/`` -- keep it distinct per
        variation so differently-constrained runs never collide.
    mc_opts : attributes set on the mc object before the kernel (optimizer
        knobs: ``max_cycle_macro``, ``conv_tol``, ...).
    """
    import types

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    if mo_full is not None:
        mc.mo_coeff = mo_full
    for k, v in (mc_opts or {}).items():
        setattr(mc, k, v)
    if mask is not None:
        mc.uniq_var_indices = types.MethodType(mask, mc)
    if weights is None:
        return _run(mc, 1, chkfile, name, verbose)
    mc.state_average_(list(weights))
    res = _run(mc, len(weights), chkfile, name, verbose)
    roots = mcscf.CASCI(mf, ncas, nelecas)
    roots.mo_coeff = mc.mo_coeff
    roots.canonicalization = False  # report the SA window exactly as converged
    rr = _run(roots, len(weights), chkfile, f"{name}_roots", verbose)
    res.energies = rr.energies
    res.core_energy = rr.core_energy
    return res


def run_casscf(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """Free CASSCF from the given start (raw-shell columns by convention)."""
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              name="casscf", verbose=verbose)


def run_sa_casscf(mf, ncas, nelecas, mo_full, weights, chkfile, verbose=0,
                  mask=None) -> CASResult:
    """SA-CASSCF over the given state weights (primary sector). PySCF's SA
    kernel returns only the ensemble energy; per-root energies come from
    the CASCI pass in the converged window -- see :func:`constrained_casscf`
    (``weights``)."""
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              weights=list(weights), mask=mask,
                              name="sa_casscf", verbose=verbose)


def run_casscf_core_relax(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with the active-window span frozen ("CASCI + core relaxation").

    Mask :func:`core_virtual_only_mask`: only inactive<->virtual
    rotations -- the core polarizes in response to the correlated CI
    density (real physics, typically ~mEh) but no orbital ever enters
    or leaves the window: state identity and window closure preserved
    by construction. Peculiarity: the relaxed core is slightly
    anisotropic; the analysis fold averages it, so the rung's
    ``max|dE|`` vs PySCF measures that anisotropy -- reported, not
    hidden.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=core_virtual_only_mask,
                              name="casscf_relax", verbose=verbose)


def run_casscf_virtual_mix(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with core<->window exchange blocked ("virtmix" rung).

    Mask :func:`no_core_window_exchange_mask`: active<->virtual and
    inactive<->virtual allowed, active<->inactive blocked -- the window
    borrows radial/polarization character from the virtuals (real
    dynamic correlation; appears as radial correlation in the natural
    orbitals) but the state's identity cannot change. relax < virtmix <
    free CASSCF brackets the contributions (README §3.3).
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=no_core_window_exchange_mask,
                              name="casscf_virtmix", verbose=verbose)


def run_casscf_pyscf_start(mf, ncas, nelecas, chkfile, verbose=0) -> CASResult:
    """CASSCF from PySCF's own default guess (canonical ROHF MOs, active =
    the canonical window ``[ncore:ncore+ncas]``), NOT the raw shell columns.
    Isolates the start-dependence of the optimizer basin (compare with
    :func:`run_casscf` on the same system)."""
    return constrained_casscf(mf, ncas, nelecas, chkfile,
                              name="casscf_pyscf", verbose=verbose)
