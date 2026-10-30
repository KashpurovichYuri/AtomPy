"""The SCF reference engine: ROHF with character-forced occupations and
the reference-basin gate, plus the fractional-occupation ensemble.

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

__all__ = ["ReferenceSpec", "run_rohf", "run_frac_occ"]


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
    conv_accelerator: str = "diis"

    def __post_init__(self):
        if (self.occ_closed is None) != (self.occ_open is None):
            raise ValueError("occ_closed and occ_open must be given together")
        if self.occ_e_ref is not None and self.occ_closed is None:
            raise ValueError(
                "occ_e_ref needs the occ_closed/occ_open patterns")
        if self.conv_accelerator not in ("diis", "newton", "diis+newton"):
            raise ValueError(
                f"unknown conv_accelerator: {self.conv_accelerator!r}")


# the runners --------------------------------------------------------

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
    the frozen-window relax rung close at 2.3e-8 there.
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

    def solve(m):
        """Run one SCF under the declared acceleration policy -- the
        ONE dispatch: "diis" plain, "newton" second-order from the
        start, "diis+newton" DIIS linearly first, then the Newton
        switch."""
        acc = ref.conv_accelerator if ref is not None else "diis"
        if acc == "newton":
            m = m.newton()
        elif acc == "diis+newton":
            m.run()
            m = m.newton()
        m.run()
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
    mf = solve(fresh_mf())
    if ref is not None and ref.occ_e_ref is not None:
        for guess in ("atom", "huckel"):
            if mf.converged and abs(mf.e_tot - ref.occ_e_ref) <= 1e-4:
                break
            mf = solve(fresh_mf(guess))
        else:
            if not (mf.converged and abs(mf.e_tot - ref.occ_e_ref) <= 1e-4):
                raise RuntimeError(
                    f"character-forced ROHF did not reach the reference "
                    f"basin E = {ref.occ_e_ref:.8f} (got {mf.e_tot:.8f}, "
                    f"converged={mf.converged}); the basin census or the "
                    f"guess list needs revisiting")
    return mol, mf


def run_frac_occ(cfg, chkfile: str, nelec_open: int, n_open: int,
                 open_idx: list[int] | None = None):
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
    mf.run()
    return mol, mf

