"""The two NEVPT2 adapters a state can be priced through.

``price_pyscf`` -- PySCF's own kernel: the external anchor (the free
routes' pricer; no production code in the path).
``price_classes`` -- the engine's own recipe with the Koopmans class
breakdown (``atomterm.perturb.nevpt2_rows_on``): the pinned route's
pricer, and the instrument that reads a class anomaly (a norm or an
energy per class).

Both consume the same :class:`~vv.routes.State`.  Their sums agree to
the kernel's own wobble (pinned in ``tests/unit/test_perturb.py`` and
here on a priced state).

THE PARTITION.  Both run the engine's prelude -- PySCF's
``canonicalize(cas_natorb=True)`` -- and that prelude RE-PARTITIONS
empty active orbitals: a natural orbital occupied below
``FRAC_OCC_THRESHOLD`` is canonicalized together with the virtuals,
sorted by orbital energy, and the lowest virtuals take its active
slots.  The priced partition is then not the solved window, and E(2)
depends on how many slots were empty rather than on which orbitals
the window named.  A frozen carbon 3P ground (its p_z empty) reads
-265.8 mEh; beside a non-coupling s orbital -165.1 mEh -- 100 mEh
from a spectator (with the window priced as solved: -264.5 and
-264.0).  The default here IS PySCF's convention (the production
engine's twin); :func:`resorted` counts the re-partitioned orbitals
and :func:`partition_as_solved` prices the window as solved.

Layering
--------
The adapters take a ready state; :meth:`vv.routes.Solution.price`
chooses one.  Strongly contracted, both.
"""

from contextlib import contextmanager

import numpy as np
from pyscf import fci, mrpt
from pyscf.mcscf import casci as _pyscf_casci

from atomterm.perturb import nevpt2_rows_on

# PySCF's cas_natorb canonicalizes the natural orbitals below this
# occupation together with the virtuals (see THE PARTITION above)
FRAC_OCC_THRESHOLD = float(_pyscf_casci.FRAC_OCC_THRESHOLD)


@contextmanager
def partition_as_solved():
    """Price the window AS SOLVED inside the block: PySCF's
    re-partitioning of empty active orbitals is suspended, and every
    solved active orbital -- occupied or empty -- takes the Dyall
    partition.  The twin of the default reading (:func:`resorted`
    says when they differ)."""
    saved = _pyscf_casci.FRAC_OCC_THRESHOLD
    _pyscf_casci.FRAC_OCC_THRESHOLD = 0.0
    try:
        yield
    finally:
        _pyscf_casci.FRAC_OCC_THRESHOLD = saved


def price_pyscf(state) -> float:
    """E(2) of the state by ``mrpt.NEVPT(mc, root).kernel()`` -- Eh
    (PySCF's partition convention; see THE PARTITION in the module
    docstring).

    On a multi-root readout the kernel writes its canonicalized
    vector back into the shared CI list; the list is detached first
    so the other roots' states stay untouched.
    """
    pt = mrpt.NEVPT(state.mc, root=state.root)
    if pt.fcisolver.nroots > 1:
        pt.ci = list(pt.ci)
    return float(pt.kernel())


def occupations(state) -> np.ndarray:
    """The state's natural occupations (the active 1-RDM's
    eigenvalues, descending; their sum is the active electron count)."""
    dm1 = fci.direct_spin1.make_rdm1(state.ci, state.ncas, state.nelecas)
    return np.sort(np.linalg.eigvalsh(np.asarray(dm1).real))[::-1]


def resorted(state) -> int:
    """How many of the state's active orbitals PySCF's NEVPT prelude
    re-partitions with the virtuals (occupation below
    ``FRAC_OCC_THRESHOLD``)."""
    return int(np.sum(occupations(state) < FRAC_OCC_THRESHOLD))


def price_classes(state):
    """E(2) of the state as one :class:`atomterm.perturb.PT2Row` --
    the eight Koopmans classes (norm and energy each), their sum the
    correction."""
    return nevpt2_rows_on(state.mc, roots=[state.root])[0]
