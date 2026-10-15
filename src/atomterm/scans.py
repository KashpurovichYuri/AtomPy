"""Frozen-orbital scans: the radial-correlation mechanism and the
flatness of the member direction, measured WITHOUT any optimizer
(ticket 16, Phase 1).

Everything here is analysis-only: a scan window is an explicitly
constructed set of columns -- per member, a rotation of one radial
family into another -- folded with the spherically averaged reference
core (the casci_avg convention) and diagonalized in the determinant
basis. Because the construction is explicit, every property (closure,
energy decomposition) is attributable to the construction, not to an
optimizer's opinion.

Two families of scans:

- m-SYMMETRIC (one angle for all members): every such window is EXACTLY
  rotation-closed -- the members share their radial vector, which IS the
  closure condition (see GLOSSARY "m-completeness"). E(theta) is the
  clean radial-correlation curve at [H, Lz] = 0, and its one-body /
  two-body split shows the mechanism: the repulsion term drops faster
  than the (nucleus + core) attraction rises.
- PER-MEMBER (independent angles): closure fails at first order in the
  angle differences while the energy stays flat -- the direct,
  optimizer-free demonstration of fact 16's
  commutator = O(eps), energy = O(eps^2) at a stationary point.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from atomterm.detbasis import (
    commutator,
    determinant_basis,
    frobenius_norm,
    slater_condon_matrix,
)
from atomterm.operators import (
    active_space_integrals,
    lz_operator,
    spin_block_one_body,
    spin_block_two_body,
)

__all__ = ["ScanPoint", "scan_window", "scan_state", "radial_scan"]


@dataclass
class ScanPoint:
    """Ground state of one frozen window: energy parts + the Lz meter."""

    theta: float                 # the scan variable (radians)
    e_active: float              # ground eigenvalue of the folded H
    e_one: float                 # <H1> (nucleus + core field)
    e_two: float                 # <H2> (active-active repulsion)
    comm_h_lz: float             # ||[H, Lz]||_F in this window

    @property
    def e_total_minus_one_two(self) -> float:
        """Bookkeeping: e_active must equal e_one + e_two exactly."""
        return abs(self.e_active - (self.e_one + self.e_two))


def scan_window(mol, mf, thetas, fam_a, fam_b) -> np.ndarray:
    """Explicit window columns ``cos(t_m) A_m + sin(t_m) B_m``.

    fam_a / fam_b : AO rows of the two radial families (same l, same
    member ordering, e.g. the three 2p AOs and the three 3p AOs).
    thetas : one mixing angle per member. Columns are S-normalized;
    different members are exactly S-orthogonal, so the window is a
    valid orthonormal one-particle set for any angles.
    """
    thetas = np.atleast_1d(thetas)
    assert len(fam_a) == len(fam_b) == len(thetas)
    n = mol.nao_nr()
    s = mf.get_ovlp()
    cols = []
    for th, a, b in zip(thetas, fam_a, fam_b):
        c = np.zeros(n)
        c[a] = np.cos(th)
        c[b] = c[b] + np.sin(th)
        c = c / float(np.real(np.sqrt(c.conj().T @ s @ c)))
        cols.append(c)
    return np.stack(cols, axis=1)


def scan_state(mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int,
               theta: float = 0.0) -> ScanPoint:
    """Diagonalize the folded H over one frozen window; decompose the
    ground state into one-body and two-body expectations.

    The one-body part is the FOLDED h1 (nucleus + averaged core field),
    so e_one rising / e_two falling along the scan is exactly the in-out
    correlation mechanism.
    """
    h1, g = active_space_integrals(mol, mf, c_act, ncore, spherical_core=True)
    dets = determinant_basis(c_act.shape[1], na, nb)
    h1m = slater_condon_matrix(dets, one_body=spin_block_one_body(h1))
    h2m = slater_condon_matrix(dets, two_body=spin_block_two_body(g))
    h = np.real(h1m + h2m)
    evals, evecs = np.linalg.eigh(h)
    v0 = evecs[:, 0]
    lz = slater_condon_matrix(
        dets, one_body=spin_block_one_body(np.asarray(lz_operator(mol, c_act))))
    return ScanPoint(
        theta=float(theta),
        e_active=float(evals[0]),
        e_one=float(np.real(v0.conj().T @ h1m @ v0)),
        e_two=float(np.real(v0.conj().T @ h2m @ v0)),
        comm_h_lz=float(frobenius_norm(commutator(h, lz))),
    )


def radial_scan(mol, mf, ncore: int, na: int, nb: int, fam_a, fam_b,
                angles) -> list[ScanPoint]:
    """m-symmetric scan: every member at the same angle (exactly closed
    windows throughout). Compare with per-member scans built directly
    from :func:`scan_window` + :func:`scan_state`."""
    out = []
    for th in angles:
        c_act = scan_window(mol, mf, [th] * len(fam_a), fam_a, fam_b)
        out.append(scan_state(mol, mf, c_act, ncore, na, nb, theta=th))
    return out


def member_tilt_map(mol, mf, ncore: int, na: int, nb: int, fam_a, fam_b,
                    theta: float, delta: float, extra_cols=None) -> dict:
    """Tilt subsets of members by ``delta`` around the symmetric point
    ``theta`` -- the flatness measurement (fact 18).

    Returns ``{"tilted <indices>": ScanPoint}``. ``extra_cols``: columns
    appended unscanned (e.g. Sc's untouched 4s slot; different-l columns
    are exactly orthogonal to the scanned family).

    Reading: on p2 atoms any SINGLE-member tilt costs exactly 0 Eh (the
    term's M_L degeneracy leaves one component untouched and the CI finds
    it) while the commutator grows linearly; tilting two or more members
    costs real energy. On Sc the structure inverts: the symmetric
    direction gains nothing, member-asymmetric tilts DESCEND (a genuine
    instability of the degenerate 2D term -- see fact 18). NOTE [H, Lz]
    is blind to tilts of the m=0 member and to equal tilts of a full
    +-m pair: Lz mixes only within +-m pairs, so those deformations leak
    only through L^2 (trivial in single-irrep maximal-Ms sectors).
    """
    n = len(fam_a)
    out = {}
    for count in range(1, n + 1):
        for subset in _subsets(n, count):
            thetas = [theta + (delta if i in subset else 0.0) for i in range(n)]
            c_act = scan_window(mol, mf, thetas, fam_a, fam_b)
            if extra_cols is not None:
                c_act = np.hstack([c_act, extra_cols])
            key = "tilt " + (",".join(map(str, subset)) if subset else "none")
            out[key] = scan_state(mol, mf, c_act, ncore, na, nb, theta=theta)
    return out


def _subsets(n, count):
    from itertools import combinations

    return combinations(range(n), count)
