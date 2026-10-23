"""The SCF reference engine: ROHF with character-forced occupations and
the reference-basin gate, plus the fractional-occupation ensemble.

Improper convergence at SCF level has two measured faces (README 3.1):
the wrong configuration (aufbau runaway) and a metastable basin at
fixed occupations. The toolbox: character-forced filling (projector
weights, never energy order), the S-projector character metric for
near-degenerate bases, an energy gate with cycled starting guesses,
and trusted cached references.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import lib, scf

from atomterm.windows import ao_span_projector, span_rows

__all__ = ["ReferenceSpec", "run_rohf", "run_frac_occ_scf"]


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
    p_closed = ao_span_projector(
        s_mat, sum((span_rows(mol, p) for p in occ_closed), []))
    p_open = ao_span_projector(
        s_mat, sum((span_rows(mol, p) for p in occ_open), []))
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
    conv_accelerator : "diis" (default; the plain SCF), "newton"
        (second-order SOSCF from the start), or "diis+newton" (DIIS
        linearly first, then the Newton switch). On heavy uncontracted
        bases with near-degenerate f/d/p valence shells, plain DIIS
        with forced occupations can fail to converge from ANY start
        (U: E > +1e3, garbage blocks), while Newton lands the target
        basin deterministically from every guess tested.
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
        acc = ref.conv_accelerator if ref is not None else "diis"
        if acc == "newton":
            m = m.newton()
        return m

    def solve(m):
        """Run one SCF under the declared acceleration policy
        ("diis+newton": DIIS linearly first, then the Newton switch)."""
        acc = ref.conv_accelerator if ref is not None else "diis"
        if acc == "diis+newton":
            m.run()
            m2 = m.newton()
            m2.run()
            return m2
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

