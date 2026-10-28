"""Unit-tier test helpers (importable as ``helpers`` from the unit test
files -- a unique basename; no conftest shadowing)."""

import numpy as np


def random_integrals(rng, norb: int, kind: str):
    """(h, g) random integrals for the construction pins, in the
    detbasis convention.

    kind="complex": raw complex draws -- g is non-particle-symmetric,
    so the two-body string contraction terms are exercised.
    kind="hermitian": real draws with h and g hermitianized (for
    matrix-level hermiticity checks).
    """
    if kind == "complex":
        h = rng.standard_normal((2 * norb, 2 * norb)) \
            + 1j * rng.standard_normal((2 * norb, 2 * norb))
        g = rng.standard_normal((2 * norb,) * 4) \
            + 1j * rng.standard_normal((2 * norb,) * 4)
    elif kind == "hermitian":
        a = rng.normal(size=(2 * norb, 2 * norb))
        h = 0.5 * (a + a.T)
        g = rng.normal(size=(2 * norb,) * 4)
        # Hermitian two-body operator: <pq|g|rs>* = <rs|g|pq>
        g = 0.5 * (g + g.transpose(2, 3, 0, 1).conj())
    else:
        raise ValueError(f"kind must be 'complex' or 'hermitian', got {kind!r}")
    return h, g
