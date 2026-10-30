"""The final analysis layer: one pass over a finished stage that
yields the campaign's canonical measured object.

Position in the architecture: the engines finish a stage (casscf
returns a :class:`~atomterm.casscf.CASResult`); the consumers declare
stages and render facts (cases: report text; treatments/cells: the
flat json payload) -- none of them computes physics. :func:`stage_facts`
is the seam between the two: the ONE composition point that folds
(:mod:`atomterm.operators`), labels (:mod:`atomterm.analysis`), and
measures the window (:mod:`atomterm.windows`) over a stage's own
outputs. The guarantee: one pass, one :class:`StageFacts`, and every
stage type reports the same fields.

Conventions
-----------
- The fold uses the stage's OWN core columns (``cas.mo_coeff[:, :ncore]``):
  orbital optimization may legally reorganize the inactive span, and a
  fold different from the one PySCF's energies used breaks
  ``E_tot = E_core + e_cas`` at O(Eh).
- The reconstruction (``eigvalsh(H) + core_energy``) is the ANALYSIS
  ruler: meters, labels, and window facts are computed on it, while the
  tables' E columns carry the run's own totals -- the two are never
  compared in a table (README · Partition bookkeeping).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy as np

from atomterm.analysis import (AnalysisResult, analyze, configuration_summary,
                               shell_occupations)
from atomterm.windows import shell_weights
from atomterm.operators import build_operators, core_anisotropy, window_closure

__all__ = ["StageFacts", "stage_facts"]

# S-orthonormality tolerance of the divergence guard: above this the
# optimizer lost the AO metric and energies stop being variational
# (observed on Si from non-block-safe window construction, ticket 09).
ORTHO_TOL = 1e-6


# the facts container ------------------------------------------------

@dataclass
class StageFacts:
    """Everything a finished stage reports, presentation-free.

    Computed by :func:`stage_facts`; the readout properties below the
    fields are thin pass-throughs (renderers never reach into the
    analysis object), and :meth:`payload` flattens the object into
    the persisted dict whose keys are committed data -- renaming one
    is a payload break.

    Attributes
    ----------
    tag : str
        Stage label, carried through to reports and stores.
    analysis : AnalysisResult
        The term-symbol analysis (states, lz_spectrum, commutators).
    energies : np.ndarray
        The reconstruction -- the analysis ruler
        (``eigvalsh(H) + core_energy``).
    energies_pyscf : np.ndarray
        The run's own totals, in the run's root order.
    max_de_vs_pyscf : float
        Max |reconstruction - run| over the common root prefix (the
        three-contributions reading: :func:`stage_facts` Notes).
    n_dev : float
        Particle-number sanity: |diag(N) - (na+nb)|.
    ortho_err : float
        Max |C^dag S C - I| over the stage's full MO matrix.
    window_names : list[str]
        The l names of the window's shell composition.
    window_rows : np.ndarray
        Per-orbital l-weights (rows sum to 1).
    win_lz, win_lperp : float
        The one-particle closure defects ||[P, Lz]|| and ||[P, L_perp]||
        of the window span (:func:`atomterm.operators.window_closure`;
        0 = the exact closure certificate).
    configurations : list[str]
        Valence-configuration label per root.
    shell_electrons : dict
        The ground root's electron count per l -- the ELECTRON scale
        next to the SPAN weights of ``window_rows``; partitions
        na + nb.
    core_anisotropy : float or None
        The discarded anisotropy of the fold core -- measured only
        when ``anisotropy_ref`` was given, else absent.
    converged : bool
        The optimizer's own kernel report; CASCI stages and trusted
        cached stages default True and never flag.
    """

    tag: str
    analysis: AnalysisResult
    energies: np.ndarray
    energies_pyscf: np.ndarray
    max_de_vs_pyscf: float
    n_dev: float
    ortho_err: float
    window_names: list[str]
    window_rows: np.ndarray
    win_lz: float
    win_lperp: float
    configurations: list[str]
    shell_electrons: dict
    core_anisotropy: float | None = None
    converged: bool = True

    # --- shared readouts (thin pass-throughs: renderers never reach
    # --- into the analysis object)
    @property
    def ground_term(self) -> str:
        return self.analysis.states[0].term

    @property
    def ground_l2(self) -> float:
        return float(self.analysis.states[0].l2)

    @property
    def ground_lz2(self) -> float:
        return float(self.analysis.states[0].lz2)

    @property
    def comm_norm_h_lz(self) -> float:
        return float(self.analysis.comm_norm_h_lz)

    @property
    def diverged(self) -> bool:
        return self.ortho_err > ORTHO_TOL

    @property
    def e_ground(self) -> float:
        return float(self.energies[0])

    @property
    def window(self) -> dict[str, float]:
        """Mean l-weights of the window columns (the l-purity gauge),
        rounded to 1e-6 -- the values persist verbatim in payloads, and
        the rounding keeps regen diffs stable."""
        return {n: round(float(v), 6)
                for n, v in zip(self.window_names, self.window_rows.mean(axis=0))}

    @property
    def terms(self) -> dict[str, int]:
        return dict(Counter(s.term for s in self.analysis.states))

    def payload(self, kind: str) -> dict:
        """Flat json-friendly dict -- the payload the stores persist.

        Parameters
        ----------
        kind : the method axis value ("casci" / "casscf" /
        "sa-casscf") --
            recorded verbatim as ``"kind"``; ``core_anisotropy``
            appears whenever it was measured.
        """
        out = {
            "tag": self.tag,
            "kind": kind,
            "n_det": self.analysis.n_det,
            "diverged": self.diverged,
            "converged": self.converged,
            "ortho_err": self.ortho_err,
            "e_ground": self.e_ground,
            "energies_pyscf": np.asarray(self.energies_pyscf).tolist(),
            "max_de_vs_pyscf": self.max_de_vs_pyscf,
            "ground_term": self.ground_term,
            "ground_l2": self.ground_l2,
            "ground_lz2": self.ground_lz2,
            "terms": self.terms,
            "comm_norm_h_lz": self.analysis.comm_norm_h_lz,
            "comm_norm_h_l2": self.analysis.comm_norm_h_l2,
            "win_lz": self.win_lz,
            "win_lperp": self.win_lperp,
            "window": self.window,
            "electrons": self.shell_electrons,
            "n_dev": self.n_dev,
            "configurations": self.configurations,
        }
        if self.core_anisotropy is not None:
            out["core_anisotropy"] = self.core_anisotropy
        return out


# the one pass -------------------------------------------------------

def stage_facts(mol, mf, cas, c_act: np.ndarray, ncore: int, na: int, nb: int,
                *, tag: str, spherical_core: bool = True,
                anisotropy_ref=None) -> StageFacts:
    """Fold and analyze one finished stage.

    Parameters
    ----------
    mol, mf : the atom's molecule and reference object (integrals are
        geometry-defined; ``mf`` supplies hcore and the overlap).
    cas : the stage's :class:`~atomterm.casscf.CASResult` (mo_coeff,
        energies, ci, core_energy).
    c_act : the window columns to analyze -- raw shell AOs for CASCI
        stages, ``cas.mo_coeff[:, ncore:ncore+ncas]`` for optimized
        stages.
    ncore, na, nb : core size and spin sector of the stage.
    tag : stage label used throughout reports and chkfile stores.
    spherical_core : fold with the spherically averaged core density
        (the physical reference) or the raw one (the raw-fold
        detector).
    anisotropy_ref : ``(mf_ref, c_core_or_None)``, optional -- measure
        the core-density anisotropy of that reference's core and report
        it (the raw-fold detector column); ``None`` leaves the
        field absent.

    Returns
    -------
    StageFacts
        The one-pass facts; ``converged`` is read off the live mc
        object when the stage carries one (the mutation discipline),
        True otherwise.

    Notes
    -----
    **The three contributions of max|dE| vs PySCF.** Up to three
    effects separate the reconstruction from the run's totals: the
    core-field anisotropy of the fold, the window DEFINITION
    (``c_act``), and bookkeeping -- which must be exactly 0 (the fold
    convention is machine-exact against PySCF's ``get_h1eff``, 1e-15;
    the pin lives in test_operators). The window-definition term is not small in
    general: where the valence shell overlaps the core, the raw-AO
    analysis window shifts the CI ground substantially against the
    run's projected columns -- on such atoms max|dE| measures the
    window definition, not the fold.

    **Which CI labels the stage.** ``cas.ci`` drives the configuration
    labels and the shell electrons; both consumers map it through the
    PySCF-ordering bridge internally. On SA stages ``cas.ci`` is the
    SA kernel's CI (the constrained_casscf contract) -- the labels
    describe the ensemble-optimized CI, not the exact per-root pass.

    Examples
    --------
    On carbon, raw 2p window, Ms=0 sector (the textbook p^2 table;
    values pinned by tests/integration/test_stages.py)::

        cas = run_casci(mf, 3, (1, 1), mo, chk, nroots=9)
        facts = stage_facts(mol, mf, cas, c_act, ncore, 1, 1, tag="c")
        facts.terms
        {'3P': 3, '1D': 5, '1S': 1}
        facts.window["p"]          # the window is exactly l-pure
        1.0
        facts.win_lz, facts.n_dev  # closure + particle number exact
        (6.5e-16, 0.0)
    """
    c_core = cas.mo_coeff[:, :ncore]
    ops = build_operators(mol, mf, c_act, ncore, na, nb,
                          spherical_core=spherical_core, c_core=c_core)
    res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag,
                  sz=ops["Sz"], lz2=ops["Lz2"])
    e_mine = np.linalg.eigvalsh(ops["H"]) + cas.core_energy
    e_pyscf = np.sort(np.atleast_1d(np.asarray(cas.energies)))[: len(e_mine)]
    w_names, w = shell_weights(mol, c_act, mf.get_ovlp())
    closure = window_closure(mol, c_act)
    mo = cas.mo_coeff
    aniso = None
    if anisotropy_ref is not None:
        mf_ref, c_aniso = anisotropy_ref
        aniso = round(core_anisotropy(mol, mf_ref, ncore, c_core=c_aniso), 12)
    return StageFacts(
        tag=tag,
        analysis=res,
        energies=e_mine,
        energies_pyscf=np.asarray(cas.energies),
        max_de_vs_pyscf=float(np.max(np.abs(e_mine[: len(e_pyscf)] - e_pyscf))),
        n_dev=float(np.max(np.abs(np.real(np.diag(ops["N"])) - (na + nb)))),
        ortho_err=float(np.abs(mo.conj().T @ mf.get_ovlp() @ mo
                               - np.eye(mo.shape[1])).max()),
        window_names=w_names,
        window_rows=w,
        win_lz=closure["lz"],
        win_lperp=closure["lperp"],
        shell_electrons=shell_occupations(mol, mf.get_ovlp(), c_act,
                                          cas.ci[0], na, nb),
        configurations=configuration_summary(
            mol, mf.get_ovlp(), cas.ci[: len(e_mine)], c_act, na, nb),
        core_anisotropy=aniso,
        # the optimizer's own report, read off the live mc object (the
        # mutation discipline); CASCI stages and cache hits carry no mc
        converged=bool(getattr(getattr(cas, "mc", None), "converged", True)),
    )
