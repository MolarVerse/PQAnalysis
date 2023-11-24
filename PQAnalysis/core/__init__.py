from .atom import Atom
from .cell import Cell
from .atomicSystem import AtomicSystem

from beartype.vale import Is
from typing import Annotated

Atoms = Annotated[list, Is[lambda list: all(
    isinstance(atom, Atom) for atom in list)]]
