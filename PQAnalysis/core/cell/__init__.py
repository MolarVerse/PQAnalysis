"""
This package contains classes and functions to handle Cell objects.
"""

from cell import Cell

from beartype.typing import NewType, Annotated
from beartype.vale import Is

#: A type hint for a list of cells
Cells = NewType(
    "Cells",
    Annotated[list,
              Is[lambda list: all(isinstance(atom, Cell) for atom in list)]]
)
