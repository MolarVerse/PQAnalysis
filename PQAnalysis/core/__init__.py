from .common import distance
from .exceptions import ElementNotFoundError, AtomicSystemPositionsError, AtomicSystemMassError
from .atom import Atom
from .cell.cell import Cell
from .atomicSystem import AtomicSystem

from beartype.vale import Is
from typing import Annotated

Atoms = Annotated[list, Is[lambda list: all(
    isinstance(atom, Atom) for atom in list)]]
