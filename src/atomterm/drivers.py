"""Compatibility facade over the split engine modules.

The former monolithic drivers module is now four concerns:

- :mod:`atomterm.reference` -- ROHF, the reference-basin gate, the
  fractional-occupation ensemble;
- :mod:`atomterm.windows` -- window construction, selection, and the
  Window registry;
- :mod:`atomterm.core` -- the purified core substrate and the fold
  energy;
- :mod:`atomterm.casscf` -- CASCI / CASSCF / SA-CASSCF.

This module re-exports the historic surface so existing imports keep
working; NEW code imports the owning module directly.
"""

from atomterm.reference import (  # noqa: F401
    ReferenceSpec,
    run_frac_occ,
    run_rohf,
)
from atomterm.windows import (  # noqa: F401
    window_avas,
    build_full_mo,
    complete_orthonormal,
    lm_sectors,
    occupied_core_idx,
    primitive_chain,
    radial_basis,
    mask_radial_rotations,
    select_primitives,
    select_shell_ao,
    select_valence_ao,
    shell_weights,
    span_overlap,
)
from atomterm.core import (  # noqa: F401
    closed_shell_core_energy,
    purified_core_orbitals,
)
from atomterm.casscf import (  # noqa: F401
    CASResult,
    constrained_casscf,
    mask_frozen_window,
    mask_no_core_window_mix,
    run_casci,
    run_casscf,
    run_sa_casscf,
)

__all__ = [
    "ReferenceSpec",
    "run_rohf",
    "run_frac_occ",
    "select_shell_ao",
    "select_valence_ao",
    "primitive_chain",
    "select_primitives",
    "occupied_core_idx",
    "shell_weights",
    "purified_core_orbitals",
    "closed_shell_core_energy",
    "complete_orthonormal",
    "build_full_mo",
    "lm_sectors",
    "radial_basis",
    "mask_radial_rotations",
    "run_casci",
    "constrained_casscf",
    "mask_frozen_window",
    "mask_no_core_window_mix",
    "run_casscf",
    "run_sa_casscf",
    "window_avas",
    "span_overlap",
    "CASResult",
]
