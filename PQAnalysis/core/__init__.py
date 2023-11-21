from .atom import Atom

from beartype.vale import Is
from typing import Annotated

Atoms = Annotated[list, Is[lambda list: all(
    isinstance(atom, Atom) for atom in list)]]
