"""One analysis pass over a finished calculation stage: the facts every
rung/stage reports, independent of presentation.

Single owner of the fold-and-analyze bookkeeping that used to live in both
``pipeline.stage_analysis`` and ``ladder._rung`` (identical physics, two
divergent dressings -- and ``core_anisotropy`` bolted on outside the fold).
``stage_facts`` computes; callers render (pipeline: report text; ladder: the
flat json payload). New rung types cannot forget a field anymore.

The fold uses the stage's OWN core columns (``cas.mo_coeff[:, :ncore]``):
orbital optimization may legally reorganize the inactive span, and a
different fold than the one PySCF's energies used breaks
``E_tot = E_core + e_cas`` at O(Eh) (fact 12).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from atomterm.analysis import AnalysisResult, analyze, configuration_summary
from atomterm.drivers import shell_weights
from atomterm.operators import build_operators, core_anisotropy, window_closure

__all__ = ["StageFacts", "stage_facts"]

# S-orthonormality tolerance of the divergence guard: above this the
# optimizer lost the AO metric and energies stop being variational
# (observed on Si from non-block-safe window construction, ticket 09).
ORTHO_TOL = 1e-6


@dataclass
class StageFacts:
    """Everything a finished stage reports, presentation-free."""

    tag: str
    analysis: AnalysisResult  # states, lz_spectrum, commutators
    energies: np.ndarray  # our reconstruction: eigvalsh(H) + core_energy
    energies_pyscf: np.ndarray
    max_de_vs_pyscf: float
    n_dev: float  # particle-number sanity: |diag(N) - (na+nb)|
    ortho_err: float  # max |C^dag S C - I| over the stage's full MO matrix
    window_names: list[str]
    window_rows: np.ndarray  # per-orbital l-weights (rows sum to 1)
    win_lz: float   # one-particle closure defects of the window SPAN
    win_lperp: float  # (operators.window_closure; 0 = exact certificate)
    configurations: list[str]

    # the three analysis fields every renderer/payload reads, exposed so
    # callers need not reach into the analysis object (interface leak
    # flagged by the design survey; kept as thin properties)
    @property
    def ground_term(self) -> str:
        return self.analysis.states[0].term

    @property
    def ground_l2(self) -> float:
        return float(self.analysis.states[0].l2)

    @property
    def comm_norm_h_lz(self) -> float:
        return float(self.analysis.comm_norm_h_lz)
    core_anisotropy: float | None = None  # CASCI reference rungs only

    @property
    def diverged(self) -> bool:
        return self.ortho_err > ORTHO_TOL

    @property
    def e_ground(self) -> float:
        return float(self.energies[0])

    @property
    def ground_term(self) -> str:
        return self.analysis.states[0].term

    @property
    def ground_l2(self) -> float:
        return self.analysis.states[0].l2

    @property
    def window(self) -> dict[str, float]:
        """Mean l-weights of the window columns (the l-purity gauge)."""
        return {n: round(float(v), 6)
                for n, v in zip(self.window_names, self.window_rows.mean(axis=0))}

    @property
    def terms(self) -> dict[str, int]:
        return dict(Counter(s.term for s in self.analysis.states))

    def payload(self, kind: str) -> dict:
        """Flat json-friendly dict (the ladder's persisted rung format)."""
        out = {
            "tag": self.tag,
            "kind": kind,
            "n_det": self.analysis.n_det,
            "diverged": self.diverged,
            "ortho_err": self.ortho_err,
            "e_ground": self.e_ground,
            "energies_pyscf": np.asarray(self.energies_pyscf).tolist(),
            "max_de_vs_pyscf": self.max_de_vs_pyscf,
            "ground_term": self.ground_term,
            "ground_l2": self.ground_l2,
            "terms": self.terms,
            "comm_norm_h_lz": self.analysis.comm_norm_h_lz,
            "comm_norm_h_l2": self.analysis.comm_norm_h_l2,
            "win_lz": self.win_lz,
            "win_lperp": self.win_lperp,
            "window": self.window,
            "n_dev": self.n_dev,
            "configurations": self.configurations,
        }
        if self.core_anisotropy is not None:
            out["core_anisotropy"] = self.core_anisotropy
        return out


def stage_facts(mol, mf, cas, c_act: np.ndarray, ncore: int, na: int, nb: int,
                *, tag: str, spherical_core: bool = True,
                anisotropy_ref=None) -> StageFacts:
    """Fold + analyze one finished stage.

    mol/mf : the atom's molecule and its reference mf object (integrals
        are geometry-defined; ``mf`` supplies hcore and the overlap).
    cas : the stage's CASResult (mo_coeff, energies, ci, core_energy).
    c_act : window columns to analyze -- raw shell AOs for CASCI stages,
        ``cas.mo_coeff[:, ncore:ncore+ncas]`` for optimized stages. NOTE
        (ticket 25): the raw-AO choice is a deliberate window DEFINITION;
        where the valence shell overlaps the core (Ti 4s|3s: overlap
        ~0.1) it shifts the CI ground by ~85 mEh vs the run's projected
        columns -- that offset, plus the mEh-class averaging effect, is
        the whole content of max_de_vs_pyscf. The fold CONVENTION itself
        is machine-exact vs PySCF's get_h1eff (1e-15, diag_fold.py).
    tag : stage label used throughout reports and chkfile stores.
    spherical_core : fold with the spherically averaged core density (the
        physical reference) or the raw one (the anisotropy detector rung).
    anisotropy_ref : ``(mf_ref, c_core_or_None)`` -- reference object whose
        core-density anisotropy should be reported (CASCI reference rungs;
        e.g. the fractional-occupation mf with its own core columns).
    """
    c_core = cas.mo_coeff[:, :ncore]
    ops = build_operators(mol, mf, c_act, ncore, na, nb,
                          spherical_core=spherical_core, c_core=c_core)
    res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag)
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
        configurations=configuration_summary(
            mol, mf.get_ovlp(), cas.ci[: len(e_mine)], c_act, na, nb),
        core_anisotropy=aniso,
    )
