"""
This is the core package of PQAnalysis.

It consists of all classes and functions that are used to handle 
atomic systems, atoms, elements and cells. Additionally it contains
some functions that are used to compute properties of atomic systems
within the api module, but can be used directly from the core package.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:  # pragma: no cover
    from .api import distance
    from .atom import Atom, Atoms, CustomElement, Element, Elements
    from .cell import Cell, Cells
    from .exceptions import (
        AtomError,
        CellError,
        ElementNotFoundError,
        ResidueError,
        ResidueWarning,
    )
    from .residue import QMResidue, Residue, Residues

_EXPORTS = {
    "ElementNotFoundError": ".exceptions",
    "ResidueError": ".exceptions",
    "ResidueWarning": ".exceptions",
    "AtomError": ".exceptions",
    "CellError": ".exceptions",
    "Cell": ".cell",
    "Cells": ".cell",
    "Atom": ".atom",
    "Atoms": ".atom",
    "Element": ".atom",
    "Elements": ".atom",
    "CustomElement": ".atom",
    "Residue": ".residue",
    "Residues": ".residue",
    "QMResidue": ".residue",
    "distance": ".api",
}

__all__ = list(_EXPORTS)



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
