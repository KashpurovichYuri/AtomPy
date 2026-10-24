"""Compat shim: the ladder lives in :mod:`atomterm.treatments` (ticket
27 -- the registry module is "treatments"; "ladder" is the historic
name). Existing imports keep working; NEW code imports the owner.
"""

from atomterm.treatments import (  # noqa: F401
    Rung,
    run_carbon_ladder,
    run_ladder,
)

__all__ = ["Rung", "run_carbon_ladder", "run_ladder"]
