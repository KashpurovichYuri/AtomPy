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
    "run_rohf",
    "run_frac_occ_scf",
    "select_shell_ao",
    "select_valence_ao",
    "shell_weights",
    "purified_s_core",
    "complete_orthonormal",
    "build_full_mo",
    "run_casci",
    "constrained_casscf",
    "core_virtual_only_mask",
    "no_core_window_exchange_mask",
    "run_casscf",
    "run_sa_casscf",
    "run_casscf_core_relax",
    "run_casscf_virtual_mix",
    "run_casscf_pyscf_start",
    "CASResult",
]


def run_rohf(cfg, chkfile: str):
    """ROHF for the configured atom, cached in ``chkfile``.

    If the chkfile exists, the converged MOs are loaded and no SCF is run.
    Returns ``(mol, mf)``.
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)
    try:
        _, rec = scf.chkfile.load_scf(chkfile)
        mf = scf.ROHF(mol)
        mf.__dict__.update(rec)
        mf.converged = True
        mf.verbose = cfg.verbose
        return mol, mf
    except (OSError, KeyError, ValueError):
        pass
    mf = scf.ROHF(mol)
    mf.chkfile = chkfile
    mf.verbose = cfg.verbose
    mf.conv_tol = 1e-12  # tight core density: keeps rotational symmetry clean
    mf.run()
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


def run_frac_occ_scf(cfg, chkfile: str, nelec_open: int, n_open: int,
                     open_idx: list[int] | None = None):
    """Average-of-configuration (fractional-occupation) HF, chkfile-cached.

    The open shell is occupied fractionally (``nelec_open / n_open`` per
    orbital, e.g. 2/3 per 2p AO for the carbon 2p^2 configuration), which
    makes the Fock operator exactly spherical. Unlike a symmetry-broken
    ROHF (whose core is polarized by the one occupied open-shell member),
    the converged core orbitals are pure l-eigenfunctions: this IS the
    rotationally invariant "scalar core" that spherically averaging a
    symmetry-broken core density reproduces to O(eps^2) in the tiny
    polarization weight eps -- while remaining a DIFFERENT physical
    reference (open-shell exchange with a smeared vs a localized valence),
    which is why the two cores' energies differ at the mEh scale.

    Control mechanism: ``get_occ`` is PySCF's occupation hook, called every
    SCF cycle to map canonical orbital energies to occupation numbers; the
    density and Fock are built from whatever it returns, so the SCF
    converges self-consistently in the ensemble density. The returned
    ``e_tot`` is the ensemble energy of the averaged configuration -- use
    this object as a CORE PROVIDER, never as a term energy.

    ``open_idx`` selects the fractional block by CHARACTER (the AO indices
    of the open shell): the ``n_open`` non-core orbitals with the largest
    weight on that shell's AOs (S metric) get the fractional occupation.
    This matters for transition metals: energy-ordered filling is
    unstable where 3d/4p quasi-degenerate -- on Sc the aufbau block runs
    away to {4s, 4p, ~2d} (a p-containing ensemble, not 3d1 4s2). Without
    ``open_idx`` the block is the HOMO block (energy order; fine for
    C/Si, where it selects the same orbitals as the projection would).

    Closed-shell RHF formalism with custom occupations; spin in ``cfg`` is
    ignored (the configuration average is an S-ensemble object — only the
    nelectron/spin PARITY of ``cfg`` is used, and it matches the atom).
    Returns ``(mol, mf)``.
    """
    from atomterm.system import build_mol

    # The caller's cfg is parity-consistent by construction (its ROHF ran
    # on it); any such spin works here because the fractional occupations
    # below define the density -- the RHF formalism never sees mol.spin.
    mol = build_mol(cfg)
    s_mat = None
    g_inv = None
    if open_idx is not None:
        assert len(open_idx) == n_open, (
            f"open_idx has {len(open_idx)} AOs but n_open={n_open}")
        s_mat = mol.intor_symmetric("int1e_ovlp")
        # AO-subspace projector metric: weight of MO i on the open shell =
        # (S c)_o^dag  S_oo^{-1}  (S c)_o  -- bounded by 1 for S-orthonormal
        # MOs (an S-projector), unlike a bare quadratic form c_o^dag S_oo c_o.
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
    # the occupation hook is part of the object's identity: set it BEFORE
    # the cache branch so a cache hit returns the same behavior
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
    l_k (S_l = S[np.ix_(rows_k, rows_k)]). S-orthonormal columns make every
    row of ``w`` sum to 1 -- this partition is valid because different-l
    subspaces are exactly orthogonal; grouping by (n, l) instead would NOT
    partition (2p and 3p AO subspaces overlap, weights double-count).

    Interpretation hierarchy: l-purity is NECESSARY but not SUFFICIENT for
    rotational closure -- a pure-p window can still be m-incomplete (radial
    mixtures lacking the full m-set), and only the commutator norms
    ||[H, Lz]||, ||[H, L^2]|| measure closure itself. Same metric quantity
    as the ``contamination`` of :func:`purified_s_core`: orbital weights
    live in the AO metric, not the Euclidean norm.
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

    Why this is needed (causal chain):

    1. Numerically converged SCF orbitals of an open-shell atom are NOT exact
       symmetry eigenfunctions: the reference density (e.g. |px py| for C) is
       anisotropic, the converged Fock inherits that anisotropy, and the
       canonical "s" MOs carry a small p/d tail (measured on C/cc-pVDZ:
       1s 1.9e-7, 2s 6.3e-6 weight deficit). The tail is numerical noise,
       not physics -- the core of a central-field atom is spherical.
    2. That tail breaks two preconditions at once:
       - orthonormality: <MO_1 | 2p_AO> ~ 1e-3 != 0, so a W built from such
         columns is not orthonormal, and second quantization (which assumes
         an orthonormal one-particle basis) is ill-defined on it;
       - rotational closure: L maps a pure solid harmonic within its (n,l)
         shell only; L acting on the p-tail leaves the span, so [H, Lz]
         fails at the ~1e-4 level instead of machine precision.
    3. Projection onto the s-AO subspace is the orbital-space analog of the
       spherically averaged density used in the folded path: it takes the
       l = 0 component exactly. The energy shifts only at second order in
       the tail amplitude (~1e-10 Eh) -- noise repair, not a model change.
    4. A Loewdin orthonormalization inside the remaining pair then restores
       orthonormality while staying inside the l = 0 subspace (rotations
       within l = 0 preserve purity).

    Returns ``(C, contamination)``: full-length (nao, n_core) orthonormal
    columns and the fraction ``1 - <c_s|S|c>`` of each canonical MO that was
    dropped (the SCF symmetry-breaking diagnostic). Note the MO s-weight
    lives in the AO metric (c^dag S c = 1), not the Euclidean norm.
    """
    import re

    labels = [line.split()[-1] for line in mol.ao_labels()]
    s_idx = [i for i, lbl in enumerate(labels) if re.fullmatch(r"\d+s", lbl)]

    c_can = mf.mo_coeff[:, :n_core]
    c_s = c_can[s_idx, :]                                # drop the p/d tail
    # MO normalization lives in the AO metric (c^dag S c = 1), so the
    # s-subspace weight is c_s^dag S_ss c_s -- not the Euclidean norm.
    s_ov = mf.get_ovlp()[np.ix_(s_idx, s_idx)]
    weight = np.einsum("ij,ij->j", c_s, s_ov @ c_s)
    contamination = np.clip(1.0 - weight, 0.0, None)     # fp noise can dip < 0
    c_s = c_s / np.sqrt(weight)                          # renormalize the rest

    gram = c_s.conj().T @ s_ov @ c_s                     # Loewdin in S metric
    ev, u = np.linalg.eigh(gram)
    c_s = c_s @ (u / np.sqrt(ev)) @ u.T

    c = np.zeros((mol.nao_nr(), n_core), dtype=c_s.dtype)
    c[s_idx, :] = c_s
    return c, contamination


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

    Construction is block-safe in the AO metric S: (1) the core block is
    Loewdin-orthonormalized within its own span (span- and density-preserving
    -- the folded core density cannot change); (2) the active columns are
    projected against the core span and renormalized, because raw valence
    AOs overlap the canonical core's same-l members (adjacent radial shells:
    <2p_MO|3p_AO> can reach ~0.3). That correction is same-l and m-diagonal,
    so an m-complete shell window keeps both its l-purity and its m-set
    structure. (3) the S-orthonormal complement is appended through an exact
    projector -- the naive `I - P` filter is only exact when v is already
    S-orthonormal, and silently returns extra columns otherwise (observed:
    (nao, nao+6) on Si).
    """
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    core = mf.mo_coeff[:, core_idx] if c_core is None else c_core
    act = np.eye(n)[:, act_idx]

    gram = core.conj().T @ s @ core
    ev, u = np.linalg.eigh(gram)
    core = core @ (u / np.sqrt(ev)) @ u.T

    p_core = core @ core.conj().T @ s
    act = act - p_core @ act
    gram_a = act.conj().T @ s @ act
    ev_a, u_a = np.linalg.eigh(gram_a)
    act = act @ (u_a / np.sqrt(ev_a)) @ u_a

    v = np.hstack([core, act])
    if v.shape[1] == n:
        return v
    return complete_orthonormal(v, s)


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
    if nroots == 1:
        energies = np.atleast_1d(float(np.atleast_1d(np.asarray(mc.e_tot)).ravel()[0]))
    else:
        energies = _root_energies(mc)
    ci = [mc.ci] if nroots == 1 else list(mc.ci)
    _persist(chkfile, name, energies, ci)
    return CASResult(name, energies, ci, mc.mo_coeff, mc, _core_energy(mc))


def run_casci(mf, ncas, nelecas, mo_full, chkfile, nroots=1, verbose=0,
              name="casci") -> CASResult:
    """CASCI in the given window. ``name`` is the chkfile store key -- keep
    it distinct per variation so differently-folded runs never collide."""
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
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
                       weights=None, verbose=0, name="casscf") -> CASResult:
    """CASSCF with declaratively constrained orbital rotations -- the one
    deep driver every orbital-optimization rung goes through.

    mo_full : starting MOs (e.g. raw shell columns via :func:`build_full_mo`);
        ``None`` = PySCF's default start (canonical MOs, canonical window) --
        the start-variant axis.
    mask : an ``uniq_var_indices``-shaped callable (see the masks above, or
        build one per window from AO labels); ``None`` = free optimization.
    weights : state-average over these weights. PySCF's SA kernel returns
        only the ENSEMBLE energy, so per-root energies come from a CASCI
        pass in the converged window -- the SA solution's own orbitals,
        diagonalizing the same active Hamiltonian the ensemble solved.
    name : chkfile store key under ``atomterm/`` -- keep it distinct per
        variation so differently-constrained runs never collide.
    """
    import types

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    if mo_full is not None:
        mc.mo_coeff = mo_full
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
    """SA-CASSCF over the given state weights (primary sector).

    Note: PySCF's state-average kernel returns only the ENSEMBLE energy
    (scalar e_cas = the weighted mean). The per-root energies reported in
    the result come from a CASCI pass in the converged window -- the SA
    solution's own orbitals, diagonalizing the same active Hamiltonian the
    ensemble was built from.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              weights=list(weights), mask=mask,
                              name="sa_casscf", verbose=verbose)


def run_casscf_core_relax(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with the active-window span frozen ("CASCI + core relaxation").

    The mask (:func:`core_virtual_only_mask`) permits only
    inactive<->virtual rotations, so the core polarizes in response to the
    correlated active CI density but no orbital ever enters or leaves the
    active window: the state's identity (and the window's rotation closure)
    is preserved by construction, while real physics (core relaxation,
    typically ~mEh) is added on top of plain CASCI.

    Peculiarity: the relaxed core is slightly anisotropic. The analysis
    fold spherically averages it, so the rung's ``max|dE|`` vs PySCF is a
    measure of that anisotropy -- reported, not hidden.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=core_virtual_only_mask,
                              name="casscf_relax", verbose=verbose)


def run_casscf_virtual_mix(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with core<->window exchange blocked ("virtmix" rung).

    Mask (:func:`no_core_window_exchange_mask`): active<->virtual AND
    inactive<->virtual rotations allowed, active<->inactive blocked. The
    window may borrow radial/polarization character from the virtual space
    -- real dynamic correlation of the active pair (radial separation is
    energetically favorable) -- but no orbital ever trades places between
    core and window, so the state's identity cannot change.

    Physics note: within the enlarged window the CI re-solves, so the
    effect appears as radial correlation (natural occupations spread below
    closed-shell values, natural orbitals delocalized over radials) -- the
    active-side analog of the core's first-order polarization by virtuals.
    Compare: rung4 (core-relax) < virtmix < free CASSCF brackets the three
    contributions.
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
