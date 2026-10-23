"""Compatibility shim -- the module is now :mod:`atomterm.cases` (the
per-atom case declaration and launcher). Existing imports keep working."""

from atomterm.cases import (  # noqa: F401
    CaseSpec,
    main,
    run_atom,
    run_carbon,
    run_case,
)

__all__ = ["CaseSpec", "run_case", "run_atom", "run_carbon", "main"]
