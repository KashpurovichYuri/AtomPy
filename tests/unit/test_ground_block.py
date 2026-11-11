"""Pure-numpy pins of the operators layer: the ground-block rule.

The block is the term's eigenspace MEASURED on an operator set -- the
selection rule's tolerance structure is the physics: the L^2 identity
carries the term (perturbed, not broken, by a residual anisotropy),
the energy window spans a fraction of the gap to the first state
OUTSIDE the L^2 band, and the exactly-degenerate case is untouched.
"""

import numpy as np
import pytest

from atomterm.operators import ground_block

pytestmark = pytest.mark.fast


def _ops(ev: list[float], l2: list[float]) -> dict:
    """A synthetic operator set: diagonal H at ``ev``, diagonal L^2``
    at ``l2`` (any eigenvector basis of both = the identity)."""
    return {"H": np.diag(ev).real, "L2": np.diag(l2).real}


def test_ground_block_spans_the_anisotropy_split():
    """A residual core-field anisotropy splits the multiplet's
    energies AND perturbs its <L^2> (the states stop being exact
    eigenstates of a non-scalar H): the block still collects the WHOLE
    term while the split stays below the window (one percent of the
    different-term gap; measured motivation: a 2D block split by
    7e-4 Eh against a 9e-2 term gap -- 0.8% -- with <L^2> drifting
    ~1e-3; the roundoff-fixed window saw one root and the NO arm
    built an m-unequal density from it)."""
    ev = ([0.0, 2e-4, 3e-4, 7e-4, 7e-4]       # the split block
          + [9e-2, 9e-2])                     # the next term
    l2 = [6.0 - 1.1e-3 * (i / 4) for i in range(5)] + [12.0, 12.0]
    e0, nroots, w, _ = ground_block(_ops(ev, l2))
    assert nroots == 5                          # the whole term
    assert abs(e0) < 1e-12
    assert np.allclose(w[:5], 0.2) and np.allclose(w[5:], 0.0)


def test_ground_block_excludes_same_l2_excited_and_keeps_exact():
    """The exactly-degenerate block is collected unchanged, and a
    same-<L^2> excited state OUTSIDE the energy window (beyond one
    percent of the different-term gap) stays out -- the window never
    grows past the term's neighborhood."""
    ev = [0.0, 0.0, 0.0, 1e-2, 1.0]            # 3-core + same-L2 + far
    l2 = [2.0, 2.0, 2.0, 2.0, 6.0]
    _, nroots, w, _ = ground_block(_ops(ev, l2))
    assert nroots == 3
    assert np.allclose(w[:3], 1 / 3) and np.allclose(w[3:], 0.0)


def test_ground_block_roundoff_when_no_other_term_near():
    """Every near state sharing the band: the window falls back to
    roundoff (1e-6) -- the certified closed-substrate behavior, the
    block exactly the degenerate eigenspace."""
    ev = [0.0, 1e-10, 1e-10, 5e-1]
    l2 = [2.0, 2.0, 2.0, 20.0]
    _, nroots, _, _ = ground_block(_ops(ev, l2))
    assert nroots == 3
