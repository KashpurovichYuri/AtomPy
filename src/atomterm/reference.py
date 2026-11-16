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
:class:`ReferenceSpec`), the SEEDED-MOM mode on ill-conditioned bases
(:class:`ReferenceSpec`'s ``seed_basis``: a projected source carries
the configuration, the occupations FOLLOW it -- the certified repair
for the fresh-SCF disease), and the fractional-occupation ensemble
(the spherical reference with smeared valence -- :func:`run_frac_occ`).

The UHF class (:func:`run_uhf`) runs the same discipline per spin:
at the campaign's maximal-Ms sectors the forced occupations have the
same per-spin counts as ROHF, so any breaking lives in the orbital
shapes -- closed pairs polarizing around the open shell -- and
:func:`uhf_contamination` measures it. The UHF object is a
REFERENCE-CLASS product here, deliberately not a cells-layer axis:
every core path consumes one closed-shell core density (``2 C C^dag``)
and one MO set, so consuming a UHF reference downstream is the
declared broken-core core step (the ``frac`` precedent does not
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
  too.  THREE filling policies share the one insertion point: PySCF's
  aufbau (no patterns), the weight hook (:func:`_forced_rohf_occ` --
  certified on contracted bases), and the MOM-follow hook
  (:func:`_mom_occ` -- the seeded mode; the seed's occupied columns
  win an S-overlap contest each cycle, so the metric-amplified
  virtual tail never competes).
- Occupation metrics are S-projector weights ``c^dag P c`` (bounded by
  1 on S-orthonormal columns); the bare quadratic form is unbounded on
  diffuse bases. One home:
  :func:`atomterm.windows.span_projector`.
- The seeded mode's honesty meter is the CHARACTER-OF-OCCUPIED
  post-check (:func:`occupation_postcheck`): the converged occupied
  columns must carry the declared spans by majority weight -- it
  never ranks all columns, because on the diseased bases the virtual
  tail carries span weights ~1 too and any ranking degrades to
  noise comparison.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import gto, scf

from atomterm.windows import (pattern_projector, project_orbitals,
                             span_projector)

__all__ = ["ReferenceSpec", "run_rohf", "run_uhf", "occupation_postcheck",
           "projected_source", "uhf_seed_density", "rohf_basin_census",
           "uhf_basin_census", "run_frac_occ", "uhf_contamination",
           "run_folded"]


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
    p_closed = pattern_projector(mol, occ_closed)
    p_open = pattern_projector(mol, occ_open)
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
        prem = _open_premium(w_open, w_closed)
        open_mos = sorted(rest, key=lambda i: -prem[i])[:n_open]
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
    p_closed = pattern_projector(mol, occ_closed)
    p_open = pattern_projector(mol, occ_open)
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
        prem = _open_premium(w_open, w_closed_a)
        open_mos = sorted(rest, key=lambda i: -prem[i])[:n_open_a]
        occ[0, open_mos] = 1.0
        return occ

    return _occ


def _mom_occ(mol, c_seed, occ_seed):
    """The MOM-anchored occupation hook (a ``get_occ`` replacement,
    the seeded mode's filler -- built for :func:`run_rohf` under
    ``ReferenceSpec.seed_basis``).

    The occupations FOLLOW the seed's occupied columns by max
    S-overlap each cycle: every current MO competes for every seeded
    slot and the highest-overlap column takes the seed's occupation.
    Near-null metric combinations are S-orthogonal to physical
    columns, so the metric-amplified virtual tail never wins the
    contest -- the seed alone carries the configuration (the weight
    hook's ranking, by contrast, degrades to noise comparison on the
    diseased bases).

    Parameters
    ----------
    mol : Mole
        The atom (the overlap metric).
    c_seed : (nao, n_seed) the seed's occupied columns (S-orthonormal
        in ``mol`` -- :func:`project_orbitals`'s product).
    occ_seed : (n_seed,) the seed's occupations per column (2.0 pairs,
        1.0 opens); assignment runs seed-order, closed pairs first.

    Returns
    -------
    callable
        ``occ(mo_energy, mo_coeff) -> occupations``.

    Notes
    -----
    The contest state RESETS each call: MO indices re-diagonalize
    every cycle, so a cross-cycle 'taken' set would assign cycle 2 to
    whatever is left (the bug the first prototype caught).
    """
    s = mol.intor_symmetric("int1e_ovlp")
    c_seed = np.asarray(c_seed)
    occ_seed = np.asarray(occ_seed)
    order = np.argsort(-occ_seed)          # closed pairs first

    def _occ(mo_energy, mo_coeff=None):
        if mo_coeff is None:
            raise ValueError("MOM occupations need the MOs")
        occ = np.zeros(len(mo_energy))
        ov = np.abs(mo_coeff.conj().T @ s @ c_seed)
        taken = np.zeros(len(mo_energy), dtype=bool)
        for j in order[occ_seed > 0]:
            free = np.flatnonzero(~taken)
            bi = free[np.argmax(ov[free, j])]
            occ[bi] = occ_seed[j]
            taken[bi] = True
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
    seed_basis : str | None
        The SEEDED-MOM mode's source basis (e.g. ``"6-31g"``) -- the
        certified repair for ill-conditioned bases where the weight
        hook's ranking degrades to noise (the fresh-SCF disease).
        Set ⇒ :func:`run_rohf` builds the projected source
        (:func:`projected_source`), seeds the density with it, and
        arms the MOM-follow hook (:func:`_mom_occ`) -- the
        occupations FOLLOW the seed's occupied columns instead of
        ranking against the span patterns; the patterns serve only
        the character post-check. Seedless ⇒ the weight hook stays
        (certified on contracted bases).
    source_closed, source_open : tuple[str, ...] | None
        The SEEDED-MOM source's own forced-hook patterns -- the
        configuration declared in the SOURCE basis's vocabulary.
        Declared whenever the two bases speak different vocabularies:
        on a contracted source the whole-l spans' closed RANK can tie
        occupied valence against virtual polarization shells and
        steal virtuals into the closed block (the rank is a weight
        statement, not a shell identity), while the label patterns
        pin the configuration. ``None`` ⇒ the source uses
        ``occ_closed``/``occ_open`` (the same vocabulary serves both
        bases -- the contracted-campaign cases).  Requires
        ``seed_basis`` (the fields are the seeded mode's own).
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
        patterns, a seed without patterns (the post-check needs the
        declared character), or an unknown accelerator.
    """

    occ_closed: tuple[str, ...] | None = None
    occ_open: tuple[str, ...] | None = None
    occ_e_ref: float | None = None
    occ_e_ref_uhf: float | None = None
    seed_basis: str | None = None
    source_closed: tuple[str, ...] | None = None
    source_open: tuple[str, ...] | None = None
    conv_accelerator: str = "diis"

    def __post_init__(self):
        if (self.occ_closed is None) != (self.occ_open is None):
            raise ValueError("occ_closed and occ_open must be given together")
        if (self.source_closed is None) != (self.source_open is None):
            raise ValueError(
                "source_closed and source_open must be given together")
        if self.occ_e_ref is not None and self.occ_closed is None:
            raise ValueError(
                "occ_e_ref needs the occ_closed/occ_open patterns")
        if self.occ_e_ref_uhf is not None and self.occ_closed is None:
            raise ValueError(
                "occ_e_ref_uhf needs the occ_closed/occ_open patterns")
        if self.seed_basis is not None and self.occ_closed is None:
            raise ValueError(
                "seed_basis needs the occ_closed/occ_open patterns "
                "(the character post-check verifies the converged "
                "solution against them)")
        if ((self.source_closed is not None)
                and self.seed_basis is None):
            raise ValueError(
                "source_closed/source_open are the seeded mode's own "
                "-- they need seed_basis (seedless, the weight hook "
                "never reads them)")
        if self.conv_accelerator not in ("diis", "newton", "diis+newton"):
            raise ValueError(
                f"unknown conv_accelerator: {self.conv_accelerator!r}")


def _open_premium(w_open, w_closed):
    """The open-block selection key: twice the open-span character plus
    the closed-span -- one shared formula across the filling policies
    (the ROHF/UHF hooks' sort keys) and the folded gates (the
    qualification floor and the character class)."""
    return 2 * w_open + w_closed


# the runners --------------------------------------------------------

def _accelerated(m, ref, dm0=None):
    """Run one SCF under the declared acceleration policy -- the ONE
    dispatch, shared by :func:`run_rohf` / :func:`run_uhf`:
    "diis" plain, "newton" second-order from the start,
    "diis+newton" DIIS linearly first, then the Newton switch.
    ``dm0`` threads the seeded mode's initial density (one signature,
    one thread; PySCF consumes it in the first Fock build)."""
    acc = ref.conv_accelerator if ref is not None else "diis"
    if acc == "newton":
        m = m.newton()
    elif acc == "diis+newton":
        m.run(dm0=dm0)
        m = m.newton()
    m.run(dm0=dm0)
    return m


def _converged_cached(mol, fresh_mf, chkfile, ref, e_gate, cfg,
                      dm0=None):
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
    mf = _accelerated(fresh_mf(None), ref, dm0)
    if e_gate is not None:
        for guess in ("atom", "huckel"):
            if mf.converged and abs(mf.e_tot - e_gate) <= 1e-4:
                break
            mf = _accelerated(fresh_mf(guess), ref, dm0)
        else:
            if not (mf.converged and abs(mf.e_tot - e_gate) <= 1e-4):
                raise RuntimeError(
                    f"character-forced SCF did not reach the reference "
                    f"basin E = {e_gate:.8f} (got {mf.e_tot:.8f}, "
                    f"converged={mf.converged}); the basin census or the "
                    f"guess list needs revisiting")
    return mol, mf


def run_rohf(cfg, chkfile: str, ref: ReferenceSpec | None = None,
             dm0=None):
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
        With ``seed_basis`` set, the SEEDED-MOM mode arms instead
        (:class:`ReferenceSpec`) -- the projected source carries the
        configuration, the occupations follow it, and the character
        post-check verifies the converged solution.
    dm0 : (nao, nao) ndarray, optional
        The caller-supplied initial density (the explicit-seed form --
        e.g. a reference chk's density feeding a spin-broken run).
        The two masters refuse: ``dm0`` together with
        ``ref.seed_basis`` is a ValueError -- the seed owns the guess
        it builds.

    Returns
    -------
    (mol, mf)
        The built molecule and the converged (or cache-loaded) object.

    Raises
    ------
    ValueError
        ``dm0`` given together with ``ref.seed_basis``.
    RuntimeError
        The seeded mode's character post-check failed on the
        converged solution (the configuration drifted from the
        seed -- the frozen-assignment disease; the basin census or
        the seed needs revisiting).

    Notes
    -----
    The gated Ti/cc-pVDZ reference lands E = -848.40648875 with
    pure-3d open MOs; its good basin carries an input-window defect
    of ~1e-8 (vs 1.6e-2 in the metastable one), which is what lets
    the frozen-window relax treatment close at 2.3e-8 there.
    """
    from atomterm.system import build_mol

    if dm0 is not None and ref is not None and ref.seed_basis is not None:
        raise ValueError(
            "dm0 and ref.seed_basis are two masters for one guess -- "
            "the seed owns the guess it builds; pass dm0 only for the "
            "caller-supplied form (seedless ref)")

    mol = build_mol(cfg)

    seeded = ref is not None and ref.seed_basis is not None
    if seeded:
        src_chk = f"{chkfile}.{ref.seed_basis}.source.chk"
        c_seed, occ_seed, _ = projected_source(
            cfg, ref, source_basis=ref.seed_basis, chkfile=src_chk)
        dm0 = (c_seed * occ_seed) @ c_seed.conj().T
        dm0 = dm0.real.astype(np.float64)

    def fresh_mf(guess: str | None = None):
        """One unconverged ROHF object: occupation hook attached,
        chkfile/verbosity set, optional starting guess."""
        m = scf.ROHF(mol)
        if seeded:
            m.get_occ = _mom_occ(mol, c_seed, occ_seed)
        elif ref is not None and ref.occ_closed is not None:
            m.get_occ = _forced_rohf_occ(mol, list(ref.occ_closed),
                                         list(ref.occ_open))
        m.chkfile = chkfile
        m.conv_tol = 1e-12  # tight core density: keeps symmetry clean
        m.verbose = cfg.verbose
        if guess is not None:
            m.init_guess = guess
        return m

    gate = ref.occ_e_ref if ref is not None else None
    mol, mf = _converged_cached(mol, fresh_mf, chkfile, ref, gate, cfg,
                                dm0=dm0)
    if seeded and not occupation_postcheck(mol, mf, ref.occ_closed,
                                           ref.occ_open):
        raise RuntimeError(
            "the seeded run's converged occupations do not carry the "
            "declared character (the frozen-assignment disease: the "
            "solution drifted from the seed's configuration) -- the "
            "basin census or the seed needs revisiting")
    return mol, mf


# the census arms: the gate's own guess cycle (default, then the two
# named starts) -- the census must see exactly what the gate can reach
ROHF_CENSUS_GUESSES = (("default", None), ("atom", "atom"),
                       ("huckel", "huckel"))


def occupation_postcheck(mol, mf, occ_closed, occ_open) -> bool:
    """The occupation post-check, CHARACTER-OF-OCCUPIED form: every
    converged occupied column must carry its block's declared span by
    MAJORITY weight (closed columns >= 1/2 on the closed spans' S
    projector; open columns >= 1/2 on the open spans').

    The well-posed form on ill-conditioned bases: it never ranks all
    columns.  The ranking-consistency predecessor (the hook
    re-applied to the converged orbitals must reproduce the run's own
    occupations) silently degrades to noise comparison exactly where
    the hook's ranking is garbage-blind -- the metric-amplified
    virtual tail carries span weights ~1 too (the fresh-SCF disease's
    own finding), so any re-ranking compares noise to noise.  The
    character form reads only the OCCUPIED columns, whose weights are
    physical in every stored record.

    Mandatory wherever an accelerator runs under a forcing or
    following hook (the census's ``occ_ok`` column;
    :func:`atomterm.basis.run_series`'s per-run gate;
    :func:`run_rohf`'s seeded mode), and on every fresh seeded run:
    energy agreement alone cannot certify the configuration -- the
    second-order accelerator freezes occupations from its early
    cycles, and a converged solution can sit in a WRONG
    configuration's basin.

    Parameters
    ----------
    mol : Mole
        The atom (the spans resolve against it).
    mf : the converged SCF object.
    occ_closed, occ_open : the declared AO span patterns
        (:class:`ReferenceSpec` semantics -- label prefixes or
        ``"l0".."l9"`` whole-momentum spans).

    Returns
    -------
    bool
        True = the converged occupied columns carry the declared
        character (the block sizes must also match the molecule's
        electron counts exactly).

    Notes
    -----
    Measured instances: a d-forced quartet converging onto the
    d²s¹ solution with the (n+1)s MO occupied at hook-rank 0 (the
    census logs); a seeded target missing its anchor by ~1.5 Eh
    when its source's closed rank stole virtuals -- the census logs
    and the wayfinder record carry both.
    """
    p_closed = pattern_projector(mol, occ_closed)
    p_open = pattern_projector(mol, occ_open)
    occ = np.asarray(mf.mo_occ)
    n_closed = (mol.nelectron - mol.spin) // 2
    n_open = mol.spin
    closed = np.flatnonzero(occ == 2.0)
    opens = np.flatnonzero(occ == 1.0)
    if len(closed) != n_closed or len(opens) != n_open:
        return False
    for i in closed:
        w = float(np.real(mf.mo_coeff[:, i].conj() @ p_closed
                          @ mf.mo_coeff[:, i]))
        if w < 0.5:
            return False
    for i in opens:
        w = float(np.real(mf.mo_coeff[:, i].conj() @ p_open
                          @ mf.mo_coeff[:, i]))
        if w < 0.5:
            return False
    return True


# the seeded mode's conditioning policy -- the ONE owner of the
# threshold (Stage 1a's verdict: conservative, removes nothing on the
# twins, the guard for genuinely singular inputs)
_SOURCE_THRESH = 1e-8


def projected_source(cfg, ref: ReferenceSpec, chkfile: str,
                     source_basis: str = "6-31g",
                     thresh: float = _SOURCE_THRESH):
    """The on-demand seeded-mode source: the declared configuration's
    forced-occupation ROHF on a SMALL basis, its occupied columns
    projected into ``cfg``'s basis -- the configuration carrier the
    MOM hook follows (:func:`run_rohf`'s seeded mode; the certified
    repair for the fresh-SCF disease).

    Parameters
    ----------
    cfg : AtomConfig
        The TARGET atom (element, basis, charge, spin -- the source
        runs the same configuration on ``source_basis``).
    ref : ReferenceSpec
        The declared occupation patterns (the projection's honesty
        meter on the target) and the source's own forced hooks
        (``source_closed``/``source_open`` when the source basis
        needs its own vocabulary -- the contracted-source case;
        falling back to ``occ_closed``/``occ_open``), plus the
        acceleration policy the whole seeded flow runs under (the
        source inherits ``conv_accelerator``: the higher-spin member
        sources are the plain-DIIS oscillation class -- the same
        measured verdict that pins ``SeriesSpec.accel``'s default --
        while the heads' default-DIIS sources are untouched).
        ``seed_basis``/``occ_e_ref`` fields are ignored here -- no
        recursion, no gate: the small basis converges tame under the
        declared policy, and :func:`occupation_postcheck` is the
        honesty meter.
    source_basis : str
        The small basis (the default is the certified source).
    chkfile : str
        The source SCF's cache path (REQUIRED -- the caller owns
        layout; :func:`run_rohf` derives it beside its own chk; the
        cache convention applies: a cached source is trusted).
    thresh : float
        The projector's conditioning cut (:func:`windows.
        project_orbitals`; the conservative default is this module's
        declared policy).

    Returns
    -------
    (c_occ_b, occ, report)
        The projected occupied columns (S-orthonormal in the target
        basis), their occupations, and the projection's disclosure
        dict (per-column capture et al.).

    Raises
    ------
    RuntimeError
        The source's own post-check failed (the small-basis SCF did
        not land the declared configuration -- the seed would carry
        a wrong configuration; fix the source, never mask it).
    """
    from dataclasses import replace as dc_replace
    from atomterm.system import build_mol

    src_cfg = dc_replace(cfg, basis=source_basis)
    src_ref = ReferenceSpec(
        occ_closed=ref.source_closed or ref.occ_closed,
        occ_open=ref.source_open or ref.occ_open,
        conv_accelerator=ref.conv_accelerator)
    mol_s, mf_s = run_rohf(src_cfg, chkfile, ref=src_ref)
    if not occupation_postcheck(mol_s, mf_s, ref.occ_closed,
                                ref.occ_open):
        raise RuntimeError(
            f"the {source_basis} source did not land the declared "
            f"configuration (the post-check failed on "
            f"{cfg.symbol}/{source_basis} spin={cfg.spin}) -- the seed "
            "would carry a wrong configuration")
    occ = np.asarray(mf_s.mo_occ)
    cols = occ > 0
    mol_b = build_mol(cfg)
    c_b, report = project_orbitals(mol_b, mol_s,
                                   mf_s.mo_coeff[:, cols], thresh=thresh)
    return c_b, occ[cols], report


def rohf_basin_census(cfg, ref: ReferenceSpec,
                      guesses=ROHF_CENSUS_GUESSES) -> list[dict]:
    """The ROHF basin census: one character-forced SCF per starting
    guess, the basin table side by side -- the measurement behind the
    configuration-series gate (:class:`atomterm.basis.SeriesSpec`'s
    ``e_ref``).

    A degenerate open shell is multi-basin at fixed occupations (the
    metastable draws differ by tens of mEh); the census MEASURES the
    basins the gate's own guess cycle can reach, the declaration
    ANCHORS the gate on one row -- the deepest converged row is the
    physical configuration state (the variational pick; an
    unconverged arm is a divergence to widen the accelerator policy
    on, never a basin). Occupations stay character-forced throughout
    (``ref``); no checkpoint is written; the arms run under ``ref``'s
    declared acceleration policy (the production path's own).

    OCCUPATION SAFETY (``occ_ok``): the Newton/second-order path
    rotates orbitals with occupations FROZEN from its early cycles --
    a converged arm can therefore carry the wrong configuration's
    occupations.  Each arm carries the occupation post-check
    (:func:`occupation_postcheck`, the character-of-occupied form):
    a converged-but-``occ_ok=False`` arm is a WRONG CONFIGURATION,
    not a basin of the intended one -- the gate never anchors there.

    Parameters
    ----------
    cfg : AtomConfig
        The atom (the series configuration's own spin).
    ref : ReferenceSpec
        The configuration's occupation forcing (without it the arms
        measure the aufbau runaway, not the configuration's basins).
    guesses : (label, init_guess) pairs
        The arms; the default = the gate's own cycle.

    Returns
    -------
    list[dict]
        Per arm: ``label``, ``e_tot``, ``converged``, ``occ_ok``.
    """
    from atomterm.system import build_mol

    if ref is None or ref.occ_closed is None:
        raise ValueError(
            "the census needs the occupation forcing (ref= with "
            "occ_closed/occ_open) -- without it the arms measure "
            "the aufbau runaway, not the configuration's basins")
    mol = build_mol(cfg)
    out = []
    for label, guess in guesses:
        m = scf.ROHF(mol)
        m.get_occ = _forced_rohf_occ(mol, list(ref.occ_closed),
                                     list(ref.occ_open))
        m.conv_tol = 1e-12
        m.verbose = cfg.verbose
        if guess is not None:
            m.init_guess = guess
        mf = _accelerated(m, ref)
        out.append({"label": label, "e_tot": float(mf.e_tot),
                    "converged": bool(mf.converged),
                    "occ_ok": occupation_postcheck(
                        mol, mf, ref.occ_closed, ref.occ_open)})
    return out


def run_uhf(cfg, chkfile: str, ref: ReferenceSpec | None = None,
            dm0=None):
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
    dm0 : (2, nao, nao) ndarray, optional
        The caller-supplied spin-separated initial density
        (:func:`uhf_seed_density`'s form -- the projected source
        split high-spin). With ``dm0`` the occupation hooks stay
        UNATTACHED: the seed owns the configuration and plain aufbau
        maintains it -- the weight hook's ranking is the
        fresh-SCF disease's own vector on ill-conditioned bases (it
        re-ranks every cycle, seed or no seed), while aufbau is
        immune (the metric-amplified virtual tail sorts last).
        Seedless runs keep the declared hooks.

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
      axis: downstream core partitions consume one closed-shell core density
      and one MO set (the module docstring's declined-seam note).
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)

    def fresh_mf(guess: str | None = None):
        """One unconverged UHF object: the per-spin occupation hook
        attached (unless a seed density owns the configuration),
        chkfile/verbosity set, optional starting guess."""
        return _uhf_object(mol, chkfile, ref, guess, verbose=cfg.verbose,
                           forced=dm0 is None)

    gate = ref.occ_e_ref_uhf if ref is not None else None
    return _converged_cached(mol, fresh_mf, chkfile, ref, gate, cfg,
                             dm0=dm0)


def uhf_seed_density(cfg, ref: ReferenceSpec, chkfile: str) -> "np.ndarray":
    """The SEEDED UHF's spin-separated initial density ``(2, nao,
    nao)``: :func:`projected_source`'s occupied columns split
    HIGH-SPIN (the closed pairs feed both spins, the open shell
    alpha) -- the explicit seed for :func:`run_uhf`'s ``dm0`` on
    ill-conditioned bases, the fresh-SCF disease's UHF class (the
    twins' anchor rows: the pre-repair cached record sat 11.4 mEh
    above the censused basin). The source chk caches beside
    ``chkfile`` (the cache convention applies)."""
    c_occ, occ, _ = projected_source(
        cfg, ref, f"{chkfile}.{ref.seed_basis}.source.chk")
    occ = np.asarray(occ)
    occ_a = np.minimum(occ, 1.0)
    occ_b = occ - occ_a
    dm_a = (c_occ * occ_a) @ c_occ.conj().T
    dm_b = (c_occ * occ_b) @ c_occ.conj().T
    return np.stack([dm_a.real, dm_b.real]).astype(np.float64)


def _uhf_object(mol, chkfile, ref, guess=None, break_=None,
                verbose: int = 0, forced: bool = True):
    """One unconverged UHF object -- the ONE construction, shared by
    :func:`run_uhf` (the production path) and
    :func:`uhf_basin_census` (the measurement arms): the per-spin
    occupation hook attached (``forced=False`` leaves it off -- the
    seeded runs' aufbau form; the seed owns the configuration),
    chkfile/verbosity set, optional starting guess and spin-break
    amplitude (PySCF's ``init_guess_break`` -- the deliberate
    symmetry breaking of the start). ``chkfile=None`` disables the
    checkpoint (the census writes nothing)."""
    m = scf.UHF(mol)
    if forced and ref is not None and ref.occ_closed is not None:
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
        certificate any future UHF core admits on.

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
                 ref: "ReferenceSpec | None" = None, dm0=None):
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
    dm0 : (nao, nao) ndarray, optional
        The caller-supplied spin-summed initial density (the seeded
        form on ill-conditioned bases -- e.g.
        :func:`projected_source`'s density; the seed anchors the
        neighborhood so the character pick lands physical valence
        columns from cycle one). The hook stays armed: the
        fractional block IS the configuration.

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
    mf = _accelerated(mf, ref, dm0)
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
               occ_closed: tuple[str, ...], occ_open: tuple[str, ...],
               dominance: float = 0.5, picks: str = "gates"):
    """The folded reference: the valence ROHF self-consistent in the
    frozen confined-core field.

    The lane's reference class.  The core's confined core
    (:func:`atomterm.core.confined_core_orbitals` -- the span's own
    closed-shell Roothaan, the ion-like field) is FROZEN, and the
    valence electrons solve ROHF in its field, in the complement of
    the confined columns.  This is the reference the lane's
    construction assumes: a window cut on it carries the ion-field
    valence, not the neutral-field one.  The valence SCF seeds from
    the substrate's own valence, projected.  Deterministic given its
    inputs on a well-behaved landscape (the carbon fixture replays
    bit-stable); the uncontracted twins' folded landscape sits at a
    basin boundary -- execution micro-state selects the basin, so a
    draw is judged by the certificates (the converged-filling
    check below, then the term meters), never assumed.

    Mechanics: PySCF's ROHF kernel runs UNMODIFIED on the complement
    basis -- an orthonormal Euclidean space (``get_ovlp = I``,
    ``get_hcore`` = the projected frozen-field h, ``get_veff`` = the
    shim translating the complement density to the AO field and
    projecting the result back).  The complement localization is what
    makes the solve converge: a pinned core block inside the cycle
    would carry a NEVER-dying error vector (the confined core
    diagonalizes its own Fock, never the total one), so the block
    must exist in no cycle object at all.  The filling is
    CHARACTER-GATED and MOM-ANCHORED: the picks follow the previous
    cycle's columns (from the seed on) within dominance gates --
    the closed pairs among the closed-span-dominated columns
    (``w > 0.5``), the opens among the d-premium columns.  Neither
    character alone nor energy alone can pick here: the complement
    carries the span's unused tight radials, whose span character
    rivals the valence pair's at energies orders of magnitude
    higher, and a field whose aufbau order may invert relative to
    the neutral atom's.  Plain DIIS.

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
    c_span : (nao, M) S-orthonormal core-span columns -- the SAME
        realization the lane's core declares.
    n_pairs : the confined pairs (the core's chemistry).
    occ_closed, occ_open : AO span patterns
        (:class:`ReferenceSpec` vocabulary) -- the character spans
        of the open-shell premium.
    dominance : the closed-pair gate of the ``gates`` picks (the
        column's closed-span weight; inert for ``three-level``).  The
        signal is GAUGE-dependent -- the same field's C1 gauge carries
        a lower 4s weight than the D2h one -- so the cut is DECLARED
        per gauge, not universal: the default 0.5 serves the D2h
        gauges; a C1 gauge whose measured signal sits below it
        declares its own cut (0.4 on the uncontracted twins, the
        measured-honesty certificates judging the result).
    picks : the valence pick rule, DECLARED per atom-gauge -- both
        forms are measured-certified and neither is universal:
        ``"gates"`` (the default) = MOM-primary within dominance
        gates; ``"three-level"`` = the character class (rounded --
        a degenerate shell ties exactly), then the energy ascending,
        then the MOM anchor.  See the Notes for which atom certified
        which form and why the other form fails there.

    Returns
    -------
    (mol, mf)
        The molecule and the folded ROHF object (converged; the
        explicit ``e_tot``; ``mo_coeff`` = [the confined core | the
        valence complement], S-orthonormal; ``mo_energy``'s frozen
        block = the confined Fock's occupied eigenvalues).  The
        accessors are restored after the solve -- the object reads
        as a standard ROHF on the real mol.  The open pair's gauge
        is the solve's own choice -- the core and window
        constructions are gauge-blind to it.

    Raises
    ------
    ValueError
        If the complement cannot host the valence filling, or the
        patterns match no AO labels.
    RuntimeError
        If the valence SCF does not converge under plain DIIS, or the
        converged filling does not hold the declared electron count.

    Notes
    -----
    The pick rule is not universal -- measured on the uncontracted
    twins.  ``three-level`` is the certified form where the valence
    closed pair is a HYBRID (its closed-span weight sits below the
    tight radials' perfect s/p character, so only the energy
    separates the classes; the ion field's aufbau inversion is why
    energy alone must not pick).  ``gates`` is the certified form
    where the closed pair is NEARLY PURE (its weight ties the tight
    radials' class and the energy-ascending leg loses the tie-break
    to them; the seed declares the basin and the gates keep every
    pick inside it).  On the wrong landscape each form abandons the
    projected basin -- and the twins' folded landscape is
    draw-chaotic besides (the basin boundary sits inside the solve;
    the certificates judge every draw).
    """
    # the core family is this runner's only consumer (a fresh local
    # import, no module-level coupling)
    from atomterm.core import (closed_shell_core_energy,
                               confined_core_orbitals)
    from atomterm.windows import complete_orthonormal

    mol = mf.mol
    n_open = mol.spin
    n_closed_val = (mol.nelectron - 2 * n_pairs - mol.spin) // 2
    if picks not in ("gates", "three-level"):
        raise ValueError(
            f"picks {picks!r} not in ('gates', 'three-level') -- the "
            "declared valence pick rules")
    if n_closed_val < 0:
        raise ValueError(
            f"{n_pairs} confined pairs leave {mol.nelectron - 2 * n_pairs} "
            f"valence electrons but the sector wants {mol.spin} open -- "
            "the core over-freezes the atom")

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
    p_open = pattern_projector(mol, occ_open)
    p_closed = pattern_projector(mol, occ_closed)
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
            f"the substrate fills {len(occ_rows)} orbitals; the core's "
            f"{n_pairs} pairs + the valence {n_seed} need more -- wrong "
            "reference for this core")
    val_cols = occ_rows[n_pairs:n_pairs + n_seed]
    c_seed = (compl.conj().T @ s @ mf.mo_coeff[:, val_cols])
    c_seed = c_seed / np.linalg.norm(c_seed, axis=0)
    dm_seed = (c_seed * np.concatenate([np.full(n_closed_val, 2.0),
                                        np.ones(n_open)])) @ c_seed.conj().T
    prev_closed[0] = c_seed[:, :n_closed_val].copy()
    prev_open[0] = c_seed[:, n_closed_val:].copy()

    hc_compl = compl.conj().T @ h_eff @ compl
    eye = np.eye(compl.shape[1])

    mol_scf = mol
    if getattr(mol, "symmetry", False):
        # the symmetry layer symmetrizes get_ovlp against the
        # AO-adapted basis -- it would strangle the complement-space
        # shims.  The folded certificates are the SO(3) closures,
        # never the D2h bookkeeping: strip the symmetry for the SCF,
        # restore the real mol after (the MO columns are AO columns,
        # valid on both).
        mol_scf = gto.M(atom=mol.atom, basis=mol.basis, spin=mol.spin,
                        ecp=mol.ecp, charge=mol.charge, verbose=0)

    fm = scf.ROHF(mol_scf)
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
        # MOM-primary picks within the dominance gates -- the rule is
        # the docstring's; the seed declared the basin this anchor
        # keeps.  A cycle's qualification can dip transiently (the
        # certified solves recover); the CONVERGED state is what must
        # hold the filling, and run_folded verifies it post-kernel.
        occ = np.zeros(len(mo_energy))
        c_ao = compl @ mo_coeff
        w_open = np.einsum("ai,ab,bi->i", c_ao.conj(), p_open,
                           c_ao).real
        w_closed = np.einsum("ai,ab,bi->i", c_ao.conj(), p_closed,
                             c_ao).real
        mom_c = np.zeros(len(occ))
        mom_o = np.zeros(len(occ))
        if prev_open[0] is not None:
            mom_o = np.abs(prev_open[0].conj().T @ mo_coeff).sum(axis=0)
        if prev_closed[0] is not None:
            mom_c = np.abs(prev_closed[0].conj().T @ mo_coeff).sum(axis=0)
        if picks == "gates":
            qualified_c = [i for i in range(len(occ))
                           if w_closed[i] > dominance]
            closed = sorted(qualified_c,
                            key=lambda i: -mom_c[i])[:n_closed_val]
            occ[closed] = 2.0
            rest = [i for i in range(len(occ)) if occ[i] == 0.0]
            qualified_o = [i for i in rest
                           if _open_premium(w_open[i], w_closed[i]) > 1.0]
            opens = sorted(qualified_o,
                           key=lambda i: -mom_o[i])[:n_open]
        else:
            # the three-level key: the CHARACTER class (rounded -- a
            # degenerate shell ties exactly), then the ENERGY
            # ascending (the tight core-shell radials carry the same
            # pure s/p character as the valence pair -- only the
            # energy separates them), then the MOM anchor
            cls_c = np.round(w_closed, 8)
            closed = sorted(range(len(occ)),
                            key=lambda i: (-cls_c[i],
                                           round(float(mo_energy[i]), 6),
                                           -mom_c[i]))[:n_closed_val]
            occ[closed] = 2.0
            rest = [i for i in range(len(occ)) if occ[i] == 0.0]
            cls_o = np.round(_open_premium(w_open, w_closed), 8)
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
        v_dv = _closed_field(mf, d_v)
        e_one = float(np.einsum("ij,ji->", h_eff, d_v).real)
        e_coul = 0.5 * float(np.einsum("ij,ji->", v_dv, d_v).real)
        fm.scf_summary["e1"] = e_one
        fm.scf_summary["e2"] = e_coul
        return e_conf + e_one + e_coul, e_coul

    fm.energy_elec = _energy_elec

    fm.kernel()
    if not fm.converged:
        raise RuntimeError(
            f"the folded valence SCF did not converge "
            f"(E = {fm.e_tot:.8f}) -- the complement-local DIIS "
            "needs revisiting")
    # the filling certificate: a cycle's gate qualification may dip
    # transiently (the MOM picks ride through and recover), but the
    # CONVERGED state must hold the declared filling -- a short-filled
    # convergence is a lost basin, not a result (the electron count of
    # the final cycle's occupancy is the honesty check)
    n_val = int(round(float(np.sum(fm.mo_occ))))
    if n_val != 2 * n_closed_val + n_open:
        raise RuntimeError(
            f"the folded valence SCF converged on {n_val} valence "
            f"electrons; the declared filling wants "
            f"{2 * n_closed_val + n_open} -- the dominance gates lost "
            "the basin")
    # reassemble on the real substrate: the pinned block prepended,
    # the accessors restored (the consumers see a standard ROHF)
    fm.mo_coeff = np.hstack([c_conf, compl @ fm.mo_coeff])
    fm.mo_energy = np.concatenate([eps_conf, fm.mo_energy])
    fm.mo_occ = np.concatenate([np.full(n_pairs, 2.0), fm.mo_occ])
    fm.mol = mol
    for name in ("get_occ", "get_init_guess", "energy_elec",
                 "get_hcore", "get_ovlp", "get_veff"):
        del fm.__dict__[name]
    return mol, fm

