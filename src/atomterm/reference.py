"""The SCF reference engine: ROHF with character-forced occupations and
the reference-basin gate, its spin-broken UHF sibling, plus the
fractional-occupation ensemble.

PURPOSE
-------
Every downstream stage consumes one (mol, mf) pair produced here, and
the reference's honesty bounds every window statement made on it.
Improper SCF convergence has two measured faces (README · Methods):
the wrong configuration (the aufbau runaway) and a metastable basin at
fixed occupations. The toolbox, escalating: character-forced filling
(S-projector weights, never energy order -- :func:`run_rohf`), the
basin gate (cycled guesses against a reference energy --
:class:`ReferenceSpec`), and the fractional-occupation ensemble (the
spherical reference with smeared valence -- :func:`run_frac_occ`).

The UHF class (:func:`run_uhf`) runs the same discipline per spin:
at the campaign's maximal-Ms sectors the forced occupations have the
same per-spin counts as ROHF, so any breaking lives in the orbital
shapes -- closed pairs polarizing around the open shell -- and
:func:`uhf_contamination` measures it. The UHF object is a
REFERENCE-CLASS product here, deliberately not a cells-layer axis:
every fold path consumes one closed-shell core density (``2 C C^dag``)
and one MO set, so consuming a UHF reference downstream is the
declared broken-core fold step (the ``frac`` precedent does not
carry: the fractional mf is RHF-shaped and drops in unchanged).

Conventions
-----------
- Caching: each run is chkfile-cached and a cached chkfile is TRUSTED
  (loaded, marked converged, never re-gated) -- the unit of
  reproducibility, because from-scratch SCF is not bit-reproducible
  (threaded BLAS draws last-digit-different MOs, and orbital-optimized
  stages amplify the draw chaotically: root splittings wobble between
  identical reruns). Equivalence checks pin the same chkfile,
  never two fresh runs.
- The occupation hook (``get_occ``) is part of the SCF object's
  identity: attached BEFORE the cache branch, so cache hits carry it
  too.
- Occupation metrics are S-projector weights ``c^dag P c`` (bounded by
  1 on S-orthonormal columns); the bare quadratic form is unbounded on
  diffuse bases. One home:
  :func:`atomterm.windows.span_projector`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import scf

from atomterm.windows import select_span, span_projector

__all__ = ["ReferenceSpec", "run_rohf", "run_uhf", "uhf_basin_census",
           "run_frac_occ", "uhf_contamination", "run_folded"]


# the occupation contract --------------------------------------------

def _forced_rohf_occ(mol, occ_closed: list[str], occ_open: list[str]):
    """The character-forced ROHF occupation hook (a ``get_occ``
    replacement built for :func:`run_rohf`).

    Parameters
    ----------
    mol : the built Mole (its nelectron/spin fix the block sizes).
    occ_closed, occ_open : list[str]
        AO span patterns (the vocabulary: :class:`ReferenceSpec`).
        Both blocks are filled by S-projector WEIGHT on those spans,
        never by energy order.

    Returns
    -------
    callable
        ``occ(mo_energy, mo_coeff) -> occupations``: the ``n_closed``
        MOs of highest closed-span weight get 2.0; then, among the
        rest, the ``n_open`` MOs ranked by

            rank = 2 w_open + w_closed

        get 1.0. The premium on open-shell character is strict: with
        whole-l spans a 1:1 sum lets closed-character leftovers TIE
        the open MOs (measured on U: DIIS then closed s/p garbage,
        E > +1e3), while the closed term still finds early-cycle
        open-shell character hiding in core-like MOs (the Sc/Ti
        leakage case).

    Raises
    ------
    ValueError
        If the hook is called without MO coefficients (character
        needs them), or a pattern matches no AO labels.
    """
    s_mat = mol.intor_symmetric("int1e_ovlp")
    p_closed = span_projector(
        s_mat, sum((select_span(mol, p) for p in occ_closed), []))
    p_open = span_projector(
        s_mat, sum((select_span(mol, p) for p in occ_open), []))
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


def _forced_uhf_occ(mol, occ_closed: list[str], occ_open: list[str]):
    """The character-forced UHF occupation hook (the per-spin twin of
    :func:`_forced_rohf_occ`, built for :func:`run_uhf`).

    Parameters
    ----------
    mol : the built Mole (its nelectron/spin fix the per-spin block
        sizes).
    occ_closed, occ_open : list[str]
        AO span patterns -- the shared :class:`ReferenceSpec`
        vocabulary, interpreted per spin.

    Returns
    -------
    callable
        ``occ(mo_energy, mo_coeff) -> (occ_alpha, occ_beta)``: per
        spin, the ``n_pair`` MOs of highest closed-span weight get
        1.0; alpha additionally takes the ``mol.spin`` unpaired
        electrons on the same ``2 w_open + w_closed`` premium as
        :func:`_forced_rohf_occ`. The beta open block is empty at the
        maximal-Ms sectors (``n_pair`` pairs + ``mol.spin`` alpha
        singles exhaust the electrons) -- a beta-open UHF reference
        would need the premium generalized, not this hook.

    Raises
    ------
    ValueError
        If the hook is called without MO coefficients (character
        needs them), or a pattern matches no AO labels.
    """
    s_mat = mol.intor_symmetric("int1e_ovlp")
    p_closed = span_projector(
        s_mat, sum((select_span(mol, p) for p in occ_closed), []))
    p_open = span_projector(
        s_mat, sum((select_span(mol, p) for p in occ_open), []))
    n_pair = (mol.nelectron - mol.spin) // 2
    n_open_a = mol.spin

    def _occ(mo_energy, mo_coeff=None):
        if mo_coeff is None:
            raise ValueError("character-forced occupations need the MOs")
        occ = np.zeros((2, mo_coeff[0].shape[1]))
        w_closed_a = None
        for s in (0, 1):
            c = mo_coeff[s]
            w_closed = np.einsum("ai,ab,bi->i", c.conj(), p_closed,
                                 c).real
            closed = np.argsort(w_closed)[::-1][:n_pair]
            occ[s, closed] = 1.0
            if s == 0:
                w_closed_a = w_closed
        c = mo_coeff[0]
        w_open = np.einsum("ai,ab,bi->i", c.conj(), p_open, c).real
        rest = np.flatnonzero(occ[0] == 0.0)
        open_mos = sorted(
            rest, key=lambda i: -(2 * w_open[i] + w_closed_a[i])
        )[:n_open_a]
        occ[0, open_mos] = 1.0
        return occ

    return _occ


@dataclass(frozen=True)
class ReferenceSpec:
    """The declared recipe for an ROHF reference that plain aufbau
    filling cannot produce: which AO spans force the occupations (by
    character, never energy order), and which energy pins the target
    basin. Threaded through run_rohf / run_case / run_atom.

    Attributes
    ----------
    occ_closed, occ_open : tuple[str, ...] | None
        AO span patterns forcing the occupations by orbital CHARACTER
        instead of energy order -- label prefixes as
        :func:`select_shell_ao`, or ``"l0".."l9"`` for a WHOLE angular
        momentum subspace (the honest vocabulary on uncontracted
        bases, where labels are shell counters). The integer analog of
        :func:`run_frac_occ`'s ``open_idx`` selection. MANDATORY
        on aufbau-runaway atoms (Ti; U needs whole-l spans plus
        Newton). Both or neither.
    occ_e_ref : float | None
        Reference energy (Eh) of the target SCF basin; gates the
        multi-basin chaos by guess cycling (see :func:`run_rohf`).
    occ_e_ref_uhf : float | None
        The UHF gate anchor: the census energy of the intended UHF
        basin (:func:`uhf_basin_census` measures it;
        :func:`run_uhf` gates on it). NOT the ROHF ``occ_e_ref`` --
        a spin-broken solution sits below it and would trip that
        gate by design.
    conv_accelerator : str
        "diis" (default; the plain SCF), "newton" (second-order SOSCF
        from the start), or "diis+newton" (DIIS linearly first, then
        the Newton switch). On heavy uncontracted bases with
        near-degenerate f/d/p valence shells, plain DIIS with forced
        occupations can fail to converge from ANY start (U: E > +1e3,
        garbage blocks), while Newton lands the target basin
        deterministically from every guess tested.

    Raises
    ------
    ValueError
        On a pattern given without its pair, a basin gate without
        patterns, or an unknown accelerator.
    """

    occ_closed: tuple[str, ...] | None = None
    occ_open: tuple[str, ...] | None = None
    occ_e_ref: float | None = None
    occ_e_ref_uhf: float | None = None
    conv_accelerator: str = "diis"

    def __post_init__(self):
        if (self.occ_closed is None) != (self.occ_open is None):
            raise ValueError("occ_closed and occ_open must be given together")
        if self.occ_e_ref is not None and self.occ_closed is None:
            raise ValueError(
                "occ_e_ref needs the occ_closed/occ_open patterns")
        if self.occ_e_ref_uhf is not None and self.occ_closed is None:
            raise ValueError(
                "occ_e_ref_uhf needs the occ_closed/occ_open patterns")
        if self.conv_accelerator not in ("diis", "newton", "diis+newton"):
            raise ValueError(
                f"unknown conv_accelerator: {self.conv_accelerator!r}")


# the runners --------------------------------------------------------

def _accelerated(m, ref):
    """Run one SCF under the declared acceleration policy -- the ONE
    dispatch, shared by :func:`run_rohf` / :func:`run_uhf`:
    "diis" plain, "newton" second-order from the start,
    "diis+newton" DIIS linearly first, then the Newton switch."""
    acc = ref.conv_accelerator if ref is not None else "diis"
    if acc == "newton":
        m = m.newton()
    elif acc == "diis+newton":
        m.run()
        m = m.newton()
    m.run()
    return m


def _converged_cached(mol, fresh_mf, chkfile, ref, e_gate, cfg):
    """The cache-then-gate skeleton shared by the SCF runners: an
    existing chkfile is loaded and TRUSTED (marked converged, never
    re-gated -- the reproducibility unit); a fresh run solves under
    ``ref``'s acceleration and, with a gate energy ``e_gate`` given,
    retries over cycled guesses (default, then "atom", then
    "huckel") until ``|E - e_gate| <= 1e-4``, else RuntimeError.
    Returns ``(mol, mf)``."""
    try:
        _, rec = scf.chkfile.load_scf(chkfile)
        mf = fresh_mf(None)
        mf.__dict__.update(rec)
        mf.converged = True
        mf.verbose = cfg.verbose
        return mol, mf
    except (OSError, KeyError, ValueError):
        pass
    mf = _accelerated(fresh_mf(None), ref)
    if e_gate is not None:
        for guess in ("atom", "huckel"):
            if mf.converged and abs(mf.e_tot - e_gate) <= 1e-4:
                break
            mf = _accelerated(fresh_mf(guess), ref)
        else:
            if not (mf.converged and abs(mf.e_tot - e_gate) <= 1e-4):
                raise RuntimeError(
                    f"character-forced SCF did not reach the reference "
                    f"basin E = {e_gate:.8f} (got {mf.e_tot:.8f}, "
                    f"converged={mf.converged}); the basin census or the "
                    f"guess list needs revisiting")
    return mol, mf


def run_rohf(cfg, chkfile: str, ref: ReferenceSpec | None = None):
    """ROHF for the configured atom, chkfile-cached.

    Parameters
    ----------
    cfg : AtomConfig
        The atom (element, basis, charge, spin).
    chkfile : str
        Cache location. An existing chkfile is loaded and TRUSTED: no
        SCF, no gate (the reproducibility unit).
    ref : ReferenceSpec, optional
        Character-forced occupations + optional basin gate.
        Occupations are assigned each cycle by S-projector weight on
        the span patterns, never by energy order -- mandatory on
        aufbau-runaway atoms. Even forced, a degenerate open shell
        is multi-basin (a metastable basin at the SAME occupations):
        with ``occ_e_ref``
        set, a fresh run retries over cycled init guesses (default,
        then "atom", then "huckel") until |E - occ_e_ref| <= 1e-4,
        else RuntimeError.

    Returns
    -------
    (mol, mf)
        The built molecule and the converged (or cache-loaded) object.

    Notes
    -----
    The gated Ti/cc-pVDZ reference lands E = -848.40648875 with
    pure-3d open MOs; its good basin carries an input-window defect
    of ~1e-8 (vs 1.6e-2 in the metastable one), which is what lets
    the frozen-window relax treatment close at 2.3e-8 there.
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)

    def fresh_mf(guess: str | None = None):
        """One unconverged ROHF object: occupation hook attached,
        chkfile/verbosity set, optional starting guess."""
        m = scf.ROHF(mol)
        if ref is not None and ref.occ_closed is not None:
            m.get_occ = _forced_rohf_occ(mol, list(ref.occ_closed),
                                         list(ref.occ_open))
        m.chkfile = chkfile
        m.conv_tol = 1e-12  # tight core density: keeps symmetry clean
        m.verbose = cfg.verbose
        if guess is not None:
            m.init_guess = guess
        return m

    gate = ref.occ_e_ref if ref is not None else None
    return _converged_cached(mol, fresh_mf, chkfile, ref, gate, cfg)


def run_uhf(cfg, chkfile: str, ref: ReferenceSpec | None = None):
    """UHF for the configured atom, chkfile-cached -- the spin-broken
    sibling of :func:`run_rohf`, same discipline (the per-spin
    character-forced hook, the basin gate, the accelerator, the
    trusted cache). Read the spin breaking with
    :func:`uhf_contamination`.

    Parameters
    ----------
    cfg, chkfile, ref : as :func:`run_rohf` (the same
        :class:`ReferenceSpec`; the occupation patterns apply per
        spin).

    Returns
    -------
    (mol, mf)
        The built molecule and the converged (or cache-loaded) UHF
        object (``mo_coeff`` is ``(2, nao, nao)``, ``mo_occ`` is
        ``(2, nao)``).

    Notes
    -----
    - At the campaign's maximal-Ms sectors the forced occupations
      have the same per-spin counts as ROHF (the open shell entirely
      alpha): any breaking is in the orbital SHAPES -- closed pairs
      polarizing around the open shell -- not a different
      configuration. The energy can then sit BELOW the ROHF basin;
      the gate anchors on :attr:`ReferenceSpec.occ_e_ref_uhf`, the
      census energy of the intended UHF basin
      (:func:`uhf_basin_census` measures it) -- the ROHF
      ``occ_e_ref`` would gate the UNBROKEN basin only and a broken
      solution trips it BY DESIGN.
    - The returned object is a REFERENCE-CLASS product, not a cells
      axis: downstream folds consume one closed-shell core density
      and one MO set (the module docstring's declined-seam note).
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)

    def fresh_mf(guess: str | None = None):
        """One unconverged UHF object: the per-spin occupation hook
        attached, chkfile/verbosity set, optional starting guess."""
        return _uhf_object(mol, chkfile, ref, guess, verbose=cfg.verbose)

    gate = ref.occ_e_ref_uhf if ref is not None else None
    return _converged_cached(mol, fresh_mf, chkfile, ref, gate, cfg)


def _uhf_object(mol, chkfile, ref, guess=None, break_=None,
                verbose: int = 0):
    """One unconverged UHF object -- the ONE construction, shared by
    :func:`run_uhf` (the production path) and
    :func:`uhf_basin_census` (the measurement arms): the per-spin
    occupation hook attached, chkfile/verbosity set, optional
    starting guess and spin-break amplitude (PySCF's
    ``init_guess_break`` -- the deliberate symmetry breaking of the
    start). ``chkfile=None`` disables the checkpoint (the census
    writes nothing)."""
    m = scf.UHF(mol)
    if ref is not None and ref.occ_closed is not None:
        m.get_occ = _forced_uhf_occ(mol, list(ref.occ_closed),
                                    list(ref.occ_open))
    m.chkfile = chkfile
    m.conv_tol = 1e-12  # tight core density: keeps symmetry clean
    m.verbose = verbose
    if guess is not None:
        m.init_guess = guess
    if break_ is not None:
        m.init_guess_break = break_
    return m


# the census default arms: the production start, the named PySCF
# starts, and one deliberately strong spin break
UHF_CENSUS_GUESSES = (("default", None, None),
                      ("minao", "minao", None),
                      ("atom", "atom", None),
                      ("huckel", "huckel", None),
                      ("broken", None, 2.0))


def uhf_basin_census(cfg, ref: ReferenceSpec | None = None,
                     guesses: tuple = UHF_CENSUS_GUESSES) -> list[dict]:
    """The UHF basin census: one converged SCF per starting guess,
    the basin table side by side (energy, spin honesty, breaking).

    The measurement behind the gate: UHF's multi-basin class is the
    widest of the reference formalisms (every closed pair is a
    breaking site), so the intended basin is DECLARED by anchoring
    :attr:`ReferenceSpec.occ_e_ref_uhf` on one census row -- the
    census measures, the declaration picks. The occupations stay
    character-forced throughout (``ref``): the basins differ in
    orbital SHAPE at fixed occupations, which is exactly what the
    ``sigma`` spectrum reads. No checkpoint is written. The arms run
    under ``ref``'s declared acceleration policy (the
    :func:`_accelerated` dispatch -- the production path's own): an
    atom whose SCF class needs the Newton accelerator diverges under
    plain DIIS in EVERY arm, which is the evidence the declaration
    rides (U14e: all five DIIS arms diverge, 2026-11-05).

    Parameters
    ----------
    cfg : AtomConfig
        The atom.
    ref : ReferenceSpec, optional
        The occupation forcing, as :func:`run_uhf`.
    guesses : (label, init_guess, init_guess_break) triples
        The arms; the default sweeps the production start, the named
        PySCF starts, and one strong spin break.

    Returns
    -------
    list[dict]
        Per arm: ``label``, ``e_tot``, ``converged``, ``s2``,
        ``sigma_min`` (the weakest corresponding-orbital pair),
        ``sigma`` (the full pairing spectrum).
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)
    out = []
    for label, guess, break_ in guesses:
        mf = _accelerated(_uhf_object(mol, None, ref, guess, break_), ref)
        s2, sigma = uhf_contamination(mf)
        out.append({"label": label, "e_tot": float(mf.e_tot),
                    "converged": bool(mf.converged), "s2": float(s2),
                    "sigma_min": float(sigma.min()) if sigma.size else 1.0,
                    "sigma": sigma.tolist()})
    return out


def uhf_contamination(mf) -> tuple[float, np.ndarray]:
    """``<S^2>`` of a UHF determinant and its pairing spectrum -- the
    reference's own spin-honesty meter.

    The overlap formula: ``<S^2> = Sz (Sz + 1) + N_beta -
    sum_ij |<phi_a_i | phi_b_j>|^2`` over the occupied blocks (the
    AO-metric overlap). A spin-paired determinant (the ROHF pattern:
    the beta occupied span inside alpha's) makes the sum ``N_beta``
    and ``<S^2> = S (S+1)`` exactly -- the formula's fixed point;
    any excess is spin contamination.

    Parameters
    ----------
    mf : a UHF object (``mo_coeff`` ``(2, nao, nao)``, ``mo_occ``
        ``(2, nao)``, ``get_ovlp``) -- converged or not.

    Returns
    -------
    (s2, sigma)
        ``s2`` the expectation value; ``sigma`` the singular values
        of the occupied overlap matrix (the corresponding-orbital
        pairs), descending: 1 = a paired couple, 0 = fully broken.
        The core pairs sitting at 1 is the spin-paired-core
        certificate any future UHF fold admits on.

    Notes
    -----
    Certified against the exact determinant-space S^2 (the unit-tier
    pins): the formula is a diagonal matrix element of
    :func:`atomterm.operators.s2_operator_parts` on the UHF
    determinant.
    """
    ca, cb = mf.mo_coeff
    occ = np.asarray(mf.mo_occ)
    a_idx = np.flatnonzero(occ[0] > 0)
    b_idx = np.flatnonzero(occ[1] > 0)
    na, nb = len(a_idx), len(b_idx)
    x = ca[:, a_idx].conj().T @ mf.get_ovlp() @ cb[:, b_idx]
    sz = (na - nb) / 2.0
    s2 = sz * (sz + 1.0) + nb - float(np.sum(np.abs(x) ** 2))
    return s2, np.linalg.svd(x, compute_uv=False)


def run_frac_occ(cfg, chkfile: str, nelec_open: int, n_open: int,
                 open_idx: list[int] | None = None,
                 ref: "ReferenceSpec | None" = None):
    """Average-of-configuration (fractional-occupation) HF,
    chkfile-cached.

    The open shell is occupied fractionally (``nelec_open / n_open``
    per orbital, e.g. 2/3 per 2p AO for the carbon 2p^2 configuration),
    which makes the converged Fock exactly spherical -- the scalar-core
    reference that density averaging reproduces, but a DIFFERENT
    physical reference from ROHF (smeared vs localized valence
    exchange; the measured energy relation: README · Motivation).

    Parameters
    ----------
    cfg, chkfile : as :func:`run_rohf` (same caching contract).
    nelec_open, n_open : electrons in, and orbital count of, the open
        shell.
    open_idx : list[int], optional
        AO indices of the open shell: the fractional block is selected
        by S-projector CHARACTER weight on that span (one home:
        :func:`atomterm.windows.span_projector`), never energy order --
        mandatory where the valence is quasi-degenerate (on Sc the
        aufbau block runs away to {4s, 4p, ~2d}, a converged but
        p-containing ensemble). ``None`` = the HOMO block by
        energy (fine for C/Si, where both rules select the same
        orbitals).
    ref : ReferenceSpec, optional
        ONLY its acceleration policy is consumed
        (:func:`_accelerated` -- an atom whose SCF class needs the
        Newton accelerator diverges under plain DIIS here exactly as
        in the gated ROHF; U14e evidence 2026-11-05). The occ patterns
        and gates do not apply: the occupations ARE the fractional
        block.

    Returns
    -------
    (mol, mf)
        Closed-shell RHF formalism; spin in ``cfg`` is ignored beyond
        nelectron parity. ``e_tot`` is the ENSEMBLE energy of the
        averaged configuration -- use the object as a core provider,
        never as a term energy.

    Notes
    -----
    Mechanism: ``get_occ`` is PySCF's per-cycle occupation hook; the
    density and Fock are built from whatever it returns, so the SCF
    converges self-consistently in the ensemble density. The selection
    is gauge-invariant WITHIN a cycle but not anchored ACROSS cycles:
    if a future atom's near-degenerate open block oscillates between
    cycles, the fix is MOM-anchored selection (follow the previous
    cycle's MOs), not a tighter conv_tol.
    """
    from atomterm.system import build_mol

    # cfg only needs the right nelectron/spin parity: the fractional
    # occupations below define the density; the RHF formalism never sees
    # mol.spin.
    mol = build_mol(cfg)
    p_open = None
    if open_idx is not None:
        if len(open_idx) != n_open:
            raise ValueError(
                f"open_idx has {len(open_idx)} AOs but n_open={n_open}")
        p_open = span_projector(
            mol.intor_symmetric("int1e_ovlp"), list(open_idx))

    def _frac_occ(mo_energy, mo_coeff=None):
        occ = np.zeros(mo_energy.size)
        n_closed = (mol.nelectron - nelec_open) // 2
        occ[:n_closed] = 2.0
        if p_open is None:
            occ[n_closed : n_closed + n_open] = nelec_open / n_open
        else:
            # valence weight of each MO: S-projector weight on the open
            # shell's AO span (rotation-invariant gauge for the block)
            c = np.asarray(mo_coeff)
            w_v = np.einsum("ai,ab,bi->i", c.conj(), p_open, c).real
            pick = np.argsort(w_v[n_closed:])[::-1][:n_open] + n_closed
            occ[pick] = nelec_open / n_open
        return occ

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.verbose = cfg.verbose
    mf.get_occ = _frac_occ
    try:
        _, rec = scf.chkfile.load_scf(chkfile)
        mf.__dict__.update(rec)
        mf.converged = True
        return mol, mf
    except (OSError, KeyError, ValueError):
        pass

    mf.chkfile = chkfile
    mf = _accelerated(mf, ref)
    return mol, mf


# the folded reference -------------------------------------------------

def _closed_field(mf, dm: np.ndarray) -> np.ndarray:
    """The spin-summed closed-shell field ``V[dm]`` of the substrate
    (PySCF's ROHF ``get_veff`` returns the per-spin pair -- the two
    entries coincide for a closed-shell density; take one)."""
    v = mf.get_veff(mf.mol, dm)
    return np.asarray(v[0] if np.ndim(v) == 3 else v)


def _spin_summed(dm: np.ndarray) -> np.ndarray:
    """PySCF's ROHF kernel hands its hooks either the spin-summed
    density or the per-spin pair, by call site -- normalize to the
    summed form."""
    d = np.asarray(dm)
    return d.sum(axis=0) if d.ndim == 3 else d


def run_folded(mf, c_span: np.ndarray, n_pairs: int,
               occ_closed: tuple[str, ...], occ_open: tuple[str, ...]):
    """The folded reference: the valence ROHF self-consistent in the
    frozen confined-core field.

    The lane's reference class.  The fold's confined core
    (:func:`atomterm.core.confined_core_orbitals` -- the span's own
    closed-shell Roothaan, the ion-like field) is FROZEN, and the
    valence electrons solve ROHF in its field, in the complement of
    the confined columns.  This is the reference the lane's
    construction assumes: a window cut on it carries the ion-field
    valence, not the neutral-field one.  Deterministic given its
    inputs (the confined solve reads no converged orbitals; the
    valence SCF seeds from the substrate's own valence, projected --
    same trusted cache, same folded reference), so it needs no
    trusted cache of its own.

    Mechanics: PySCF's ROHF kernel runs UNMODIFIED on the complement
    basis -- an orthonormal Euclidean space (``get_ovlp = I``,
    ``get_hcore`` = the projected frozen-field h, ``get_veff`` = the
    shim translating the complement density to the AO field and
    projecting the result back).  The complement localization is what
    makes the solve converge: a pinned core block inside the cycle
    would carry a NEVER-dying error vector (the confined core
    diagonalizes its own Fock, never the total one), so the block
    must exist in no cycle object at all.  The filling is
    CHARACTER-FORCED (the :func:`_forced_rohf_occ` rule, restricted
    to the complement): the valence closed pairs by closed-span
    weight, then the opens by the ``2 w_open + w_closed`` premium,
    three-level picks (class, then energy, then the MOM anchor
    across cycles -- the :func:`run_frac_occ` cure for the
    degenerate-shell gauge chase).  The solve SEEDS from the
    substrate's own valence orbitals, projected against the confined
    core -- the basin-census cure (this atom's d^2 open-shell
    problems diverge from every scratch start; the projected valence
    starts the DIIS inside the declared basin).  Plain DIIS.

    The bookkeeping (the total, computed explicitly each cycle):

        E = E_conf + Tr[D_v h_eff] + 1/2 Tr[D_v V[D_v]]

    with ``E_conf`` = :func:`atomterm.core.closed_shell_core_energy`
    of the confined core and ``h_eff = h + V[D_conf]`` -- the core
    energy enters once, the valence pays its interaction with the
    frozen field once, no double counting.  The variational bound:
    the frozen core constrains the full same-filling problem, so
    the folded total sits at or above the full ROHF.

    Parameters
    ----------
    mf : the reference substrate (``mol``, ``get_hcore``,
        ``get_ovlp``, ``get_veff``; no converged orbitals are read).
    c_span : (nao, M) S-orthonormal fold-span columns -- the SAME
        realization the lane's fold declares.
    n_pairs : the confined pairs (the fold's chemistry).
    occ_closed, occ_open : AO span patterns
        (:class:`ReferenceSpec` vocabulary) -- the character spans
        of the open-shell premium.

    Returns
    -------
    (mol, mf)
        The molecule and the folded ROHF object (converged; the
        explicit ``e_tot``; ``mo_coeff`` = [the confined core | the
        valence complement], S-orthonormal; ``mo_energy``'s frozen
        block = the confined Fock's occupied eigenvalues).  The
        accessors are restored after the solve -- the object reads
        as a standard ROHF on the real mol.  The open pair's gauge
        is the solve's own choice -- the fold and window
        constructions are gauge-blind to it.

    Raises
    ------
    ValueError
        If the complement cannot host the valence filling, or the
        patterns match no AO labels.
    RuntimeError
        If the valence SCF does not converge under plain DIIS.
    """
    from pyscf import scf

    from atomterm.core import (closed_shell_core_energy,
                               confined_core_orbitals)
    from atomterm.windows import (complete_orthonormal, select_span,
                                  span_projector)

    mol = mf.mol
    n_open = mol.spin
    n_closed_val = (mol.nelectron - 2 * n_pairs - mol.spin) // 2
    if n_closed_val < 0:
        raise ValueError(
            f"{n_pairs} confined pairs leave {mol.nelectron - 2 * n_pairs} "
            f"valence electrons but the sector wants {mol.spin} open -- "
            "the fold over-freezes the atom")

    c_conf, _ = confined_core_orbitals(mf, c_span, n_pairs)
    s = mf.get_ovlp()
    dm_conf = 2.0 * (c_conf @ c_conf.conj().T)
    h_eff = mf.get_hcore() + _closed_field(mf, dm_conf)
    e_conf = closed_shell_core_energy(mf, c_conf)

    # the complement: the S-orthonormal columns completing the frozen block
    compl = complete_orthonormal(c_conf, s)[:, n_pairs:]
    if compl.shape[1] < n_closed_val + n_open:
        raise ValueError(
            f"the complement of {n_pairs} confined columns holds "
            f"{compl.shape[1]} columns; the valence filling wants "
            f"{n_closed_val + n_open} occupied")
    p_open = span_projector(
        s, sum((select_span(mol, p) for p in occ_open), []))
    p_closed = span_projector(
        s, sum((select_span(mol, p) for p in occ_closed), []))
    eps_conf = np.einsum("ai,ab,bi->i", c_conf.conj(), h_eff,
                         c_conf).real
    prev_open = [None]   # the previous cycle's open columns (the MOM anchor)
    prev_closed = [None]  # ditto for the valence closed pair

    # the seed: the substrate's own valence orbitals, projected against
    # the confined core (normalized -- the projection removes the
    # confined-core content).  The MOM anchors start on the seed, so
    # the DIIS never leaves the projected basin.
    occ_rows = [i for i, o in enumerate(mf.mo_occ) if o > 0]
    n_seed = n_closed_val + n_open
    if len(occ_rows) < n_pairs + n_seed:
        raise ValueError(
            f"the substrate fills {len(occ_rows)} orbitals; the fold's "
            f"{n_pairs} pairs + the valence {n_seed} need more -- wrong "
            "reference for this fold")
    val_cols = occ_rows[n_pairs:n_pairs + n_seed]
    c_seed = (compl.conj().T @ s @ mf.mo_coeff[:, val_cols])
    c_seed = c_seed / np.linalg.norm(c_seed, axis=0)
    dm_seed = (c_seed * np.concatenate([np.full(n_closed_val, 2.0),
                                        np.ones(n_open)])) @ c_seed.conj().T
    prev_closed[0] = c_seed[:, :n_closed_val].copy()
    prev_open[0] = c_seed[:, n_closed_val:].copy()

    hc_compl = compl.conj().T @ h_eff @ compl
    eye = np.eye(compl.shape[1])

    fm = scf.ROHF(mol)
    fm.conv_tol = 1e-12
    fm.verbose = 0
    fm.get_hcore = lambda *args, **kw: hc_compl
    fm.get_ovlp = lambda *args, **kw: eye

    def _veff(mol_=None, dm=None, dm_last=0, v_last=None, hermi=1):
        d_ao = compl @ _spin_summed(dm) @ compl.conj().T
        v = compl.conj().T @ _closed_field(mf, d_ao) @ compl
        return np.array((v, v))

    fm.get_veff = _veff

    def _init_guess(mol_=None, key=None, **kw):
        # the projected-valence seed (see above): the substrate's own
        # valence, confined-core-orthogonal, in the declared filling
        return dm_seed

    fm.get_init_guess = _init_guess

    def _occ(mo_energy, mo_coeff):
        # the filling is CHARACTER-FORCED (the _forced_rohf_occ rule,
        # restricted to the complement): the closed pairs by
        # closed-span weight, then the opens by the premium.
        occ = np.zeros(len(mo_energy))
        c_ao = compl @ mo_coeff
        w_open = np.einsum("ai,ab,bi->i", c_ao.conj(), p_open,
                           c_ao).real
        w_closed = np.einsum("ai,ab,bi->i", c_ao.conj(), p_closed,
                             c_ao).real
        # the picks are three-level: the CHARACTER class (rounded -- a
        # degenerate shell ties exactly), then the ENERGY ascending
        # (the tight core-shell radials in the complement carry the
        # same pure s/p character as the valence pair -- only the
        # energy separates them; the ion field's d-below-s inversion
        # is why the closed pair must not be picked by energy ALONE),
        # then the MOM anchor (the previous cycle's columns -- the
        # run_frac_occ cure for the degenerate-shell gauge chase).
        mom_c = np.zeros(len(occ))
        mom_o = np.zeros(len(occ))
        if prev_open[0] is not None:
            mom_o = np.abs(prev_open[0].conj().T @ mo_coeff).sum(axis=0)
        if prev_closed[0] is not None:
            mom_c = np.abs(prev_closed[0].conj().T @ mo_coeff).sum(axis=0)
        cls_c = np.round(w_closed, 8)
        closed = sorted(range(len(occ)),
                        key=lambda i: (-cls_c[i],
                                       round(float(mo_energy[i]), 6),
                                       -mom_c[i]))[:n_closed_val]
        occ[closed] = 2.0
        rest = [i for i in range(len(occ)) if occ[i] == 0.0]
        cls_o = np.round(2 * w_open + w_closed, 8)
        opens = sorted(rest,
                       key=lambda i: (-cls_o[i],
                                      round(float(mo_energy[i]), 6),
                                      -mom_o[i]))[:n_open]
        occ[opens] = 1.0
        prev_open[0] = mo_coeff[:, opens].copy()
        prev_closed[0] = mo_coeff[:, closed].copy()
        return occ

    fm.get_occ = _occ

    def _energy_elec(dm=None, h1e=None, vhf=None):
        d_v = compl @ _spin_summed(dm) @ compl.conj().T
        e_val = (float(np.einsum("ij,ji->", h_eff, d_v).real)
                 + 0.5 * float(
                     np.einsum("ij,ji->", _closed_field(mf, d_v),
                               d_v).real))
        fm.scf_summary["e1"] = e_val
        fm.scf_summary["e2"] = 0.5 * e_val
        return e_conf + e_val, e_val

    fm.energy_elec = _energy_elec

    fm.kernel()
    if not fm.converged:
        raise RuntimeError(
            f"the folded valence SCF did not converge "
            f"(E = {fm.e_tot:.8f}) -- the complement-local DIIS "
            "needs revisiting")
    # reassemble on the real substrate: the pinned block prepended,
    # the accessors restored (the consumers see a standard ROHF)
    fm.mo_coeff = np.hstack([c_conf, compl @ fm.mo_coeff])
    fm.mo_energy = np.concatenate([eps_conf, fm.mo_energy])
    fm.mo_occ = np.concatenate([np.full(n_pairs, 2.0), fm.mo_occ])
    for name in ("get_occ", "get_init_guess", "energy_elec",
                 "get_hcore", "get_ovlp", "get_veff"):
        del fm.__dict__[name]
    return mol, fm

