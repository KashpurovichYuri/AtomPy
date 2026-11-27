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

Layering
--------
The adapters take a ready state; :meth:`vv.routes.Solution.price`
chooses one.  Strongly contracted, both.
"""

from pyscf import mrpt

from atomterm.perturb import nevpt2_rows_on


def price_pyscf(state) -> float:
    """E(2) of the state by ``mrpt.NEVPT(mc, root).kernel()`` -- Eh.

    On a multi-root readout the kernel writes its canonicalized
    vector back into the shared CI list; the list is detached first
    so the other roots' states stay untouched.
    """
    pt = mrpt.NEVPT(state.mc, root=state.root)
    if pt.fcisolver.nroots > 1:
        pt.ci = list(pt.ci)
    return float(pt.kernel())


def price_classes(state):
    """E(2) of the state as one :class:`atomterm.perturb.PT2Row` --
    the eight Koopmans classes (norm and energy each), their sum the
    correction."""
    return nevpt2_rows_on(state.mc, roots=[state.root])[0]
