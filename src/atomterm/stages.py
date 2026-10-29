"""One analysis pass over a finished stage: the facts every renderer
reports, independent of presentation.

``stage_facts`` is the single owner of the fold-and-analyze bookkeeping;
callers only render (cases: report text; treatments/cells: the flat json
payload). One pass yields one :class:`StageFacts`, so every stage type
reports the same fields.

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


@dataclass
class StageFacts:
    """Everything a finished stage reports, presentation-free.

    Values are computed by :func:`stage_facts`; the properties below the
    fields are the thin readouts renderers share, and :meth:`payload`
    flattens the object into the persisted rungs/cells dict -- whose
    keys are committed data, so renaming one is a payload break.
    """

    tag: str
    analysis: AnalysisResult  # states, lz_spectrum, commutators
    energies: np.ndarray  # the reconstruction (analysis ruler): eigvalsh(H) + core_energy
    energies_pyscf: np.ndarray  # the run's own totals, in the run's root order
    max_de_vs_pyscf: float  # max |reconstruction - run| over the common root prefix
    n_dev: float  # particle-number sanity: |diag(N) - (na+nb)|
    ortho_err: float  # max |C^dag S C - I| over the stage's full MO matrix
    window_names: list[str]
    window_rows: np.ndarray  # per-orbital l-weights (rows sum to 1)
    win_lz: float  # one-particle closure defect ||[P, Lz]|| of the window span
    win_lperp: float  # ||[P, L_perp]|| (window_closure; 0 = exact certificate)
    configurations: list[str]
    # the ground root's electron count per l (analysis.
    # shell_occupations) -- the ELECTRON scale next to the SPAN
    # weights of window_rows; partitions na+nb
    shell_electrons: dict
    core_anisotropy: float | None = None  # CASCI reference rungs only
    # converged: the optimizer's own kernel report; CASCI stages and
    # trusted cached stages default True and never flag
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
        """Flat json-friendly dict -- the persisted rungs/cells payload.

        Parameters
        ----------
        kind : the method axis value ("casci" / "casscf" / "sa") --
            recorded verbatim as ``"kind"``; CASCI kinds also carry
            ``core_anisotropy`` when it was measured.
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
        (the physical reference) or the raw one (the anisotropy
        detector rung).
    anisotropy_ref : ``(mf_ref, c_core_or_None)``, optional -- measure
        the core-density anisotropy of that reference's core and report
        it (the CASCI rungs' detector column); ``None`` leaves the
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
    scripts/diag_fold.py). The window-definition term is not small in
    general: where the valence shell overlaps the core (Ti 4s|3s
    overlap ~0.1), the raw-AO analysis window shifts the CI ground by
    ~85 mEh against the run's projected columns -- on such atoms
    max|dE| measures the window definition, not the fold.

    **Which CI labels the stage.** ``cas.ci`` drives the configuration
    labels and the shell electrons; both consumers map it through the
    PySCF-ordering bridge internally. On SA stages ``cas.ci`` is the
    SA kernel's CI (the constrained_casscf contract) -- the labels
    describe the ensemble-optimized CI, not the exact per-root pass.
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
