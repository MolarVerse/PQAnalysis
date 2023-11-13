"""
A module containing a Mixin Class with the standard properties of an atomic system (i.e. standard getter and setter methods).

...

Classes
-------
_StandardPropertiesMixin
    A mixin class containing the standard properties of an atomic system (i.e. standard getter and setter methods).
"""

from beartype.typing import List

from ..atom import Atom
from ..cell import Cell
from ...types import Np1DNumberArray, Np2DNumberArray


class _StandardPropertiesMixin:
    """
    A mixin class containing the standard properties of an atomic system (i.e. standard getter and setter methods).
    """
    @property
    def atoms(self) -> List[Atom]:
        """
        Returns the atoms in the system.

        Returns
        -------
        List[Atom]
            The atoms in the system.
        """
        return self._atoms

    @property
    def pos(self) -> Np2DNumberArray:
        """
        Returns the positions of the atoms in the system.

        Returns
        -------
        Np2DNumberArray
            The positions of the atoms in the system.
        """
        return self._pos

    @property
    def vel(self) -> Np2DNumberArray:
        """
        Returns the velocities of the atoms in the system.

        Returns
        -------
        Np2DNumberArray
            The velocities of the atoms in the system.
        """
        return self._vel

    @property
    def forces(self) -> Np2DNumberArray:
        """
        Returns the forces on the atoms in the system.

        Returns
        -------
        Np2DNumberArray
            The forces on the atoms in the system.
        """
        return self._forces

    @property
    def charges(self) -> Np1DNumberArray:
        """
        Returns the charges of the atoms in the system.

        Returns
        -------
        Np1DNumberArray
            The charges of the atoms in the system.
        """
        return self._charges

    @property
    def cell(self) -> Cell:
        """
        Returns the unit cell of the system.

        Returns
        -------
        Cell
            The unit cell of the system.
        """
        return self._cell

    @cell.setter
    def cell(self, cell: Cell) -> None:
        """
        Sets the unit cell of the system.

        Parameters
        ----------
        cell : Cell
            The unit cell of the system.
        """
        self._cell = cell
