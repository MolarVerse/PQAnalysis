"""
A module containing the AtomicSystem class

...

Classes
-------
AtomicSystem
    A class for storing atomic systems.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import List, Any
from multimethod import multimethod

from ._decorators import check_atoms_pos, check_atoms_has_mass
from ._properties import _PropertiesMixin
from ._standardProperties import _StandardPropertiesMixin
from ._indexing import _IndexingMixin
from ._positions import _PositionsMixin

from .. import Atom, Cell
from ...types import Np2DNumberArray, Np1DNumberArray, Np1DIntArray


class AtomicSystem(_PropertiesMixin, _StandardPropertiesMixin, _IndexingMixin, _PositionsMixin):
    """
    A class for storing atomic systems.

    Inherits from the Mixins: _PropertiesMixin, _StandardPropertiesMixin, _IndexingMixin, _PositionsMixin
        - The _StandardPropertiesMixin contains the standard properties of an atomic system (i.e. standard getter and setter methods).
        - The _PropertiesMixin contains special properties derived from the standard properties
        - The _IndexingMixin contains methods for indexing the atomic system
        - The _PositionsMixin contains methods for computing properties based on the positions of the atoms
    """

    def __init__(self,
                 atoms: List[Atom] | None = None,
                 pos: Np2DNumberArray | None = None,
                 vel: Np2DNumberArray | None = None,
                 forces: Np2DNumberArray | None = None,
                 charges: Np1DNumberArray | None = None,
                 cell: Cell = Cell()
                 ) -> None:
        """
        Initializes the AtomicSystem with the given parameters.

        Parameters
        ----------
        atoms : List[Atom], optional
            A list of Atom objects, by default []
        pos : Np2DNumberArray, optional
            A 2d numpy.ndarray containing the positions of the atoms, by default np.zeros((0, 3)).
        vel : Np2DNumberArray, optional
            A 2d numpy.ndarray containing the velocities of the atoms, by default np.zeros((0, 3)).
        forces : Np2DNumberArray, optional
            A 2d numpy.ndarray containing the forces on the atoms, by default np.zeros((0, 3)).
        charges : Np1DNumberArray, optional
            A 1d numpy.ndarray containing the charges of the atoms, by default np.zeros(0).
        cell : Cell, optional
            The unit cell of the system. Defaults to a Cell with no periodic boundary conditions, by default Cell()
        """
        if atoms is None:
            atoms = []

        if pos is None:
            pos = np.zeros((0, 3))

        if vel is None:
            vel = np.zeros((0, 3))

        if forces is None:
            forces = np.zeros((0, 3))

        if charges is None:
            charges = np.zeros(0)

        self._atoms = atoms
        self._pos = pos
        self._vel = vel
        self._forces = forces
        self._charges = charges
        self._cell = cell

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the AtomicSystem is equal to another AtomicSystem.

        Parameters
        ----------
        other : AtomicSystem
            The other AtomicSystem to compare to.

        Returns
        -------
        bool
            Whether the AtomicSystem is equal to the other AtomicSystem.
        """
        if not isinstance(other, AtomicSystem):
            return False

        elif self.n_atoms != other.n_atoms:
            return False

        elif self.cell != other.cell:
            return False

        elif self.atoms != other.atoms:
            return False

        is_equal = True
        is_equal &= np.allclose(self.pos, other.pos)
        is_equal &= np.allclose(self.vel, other.vel)
        is_equal &= np.allclose(self.forces, other.forces)
        is_equal &= np.allclose(self.charges, other.charges)

        return is_equal

    # TODO: add possibility to index with atom list etc... similar to nearest neighbours
    @multimethod
    def __getitem__(self, key: Atom) -> AtomicSystem:
        """
        Returns a new AtomicSystem with the given key.

        Parameters
        ----------
        key : Atom
            The key as an atom to get the new AtomicSystem with only the matching atoms.

        Returns
        -------
        AtomicSystem
            The new AtomicSystem with the given key.
        """
        indices = np.argwhere(np.array(self.atoms) == key).flatten()

        return self.__getitem__(indices)

    @multimethod
    def __getitem__(self, key: int | slice | Np1DIntArray) -> AtomicSystem:
        """
        Returns a new AtomicSystem with the given key.

        Parameters
        ----------
        key : int | slice | Np1DIntArray
            The key to get the new AtomicSystem with.

        Returns
        -------
        AtomicSystem
            The new AtomicSystem with the given key.
        """

        if isinstance(key, int):
            keys = np.array([key])
        elif isinstance(key, slice):
            keys = np.array(range(self.n_atoms)[key])
        else:
            keys = np.array(key)

        if self.atoms != []:
            atoms = [self.atoms[key] for key in keys]
        else:
            atoms = None

        if np.shape(self.pos)[0] > 0:
            pos = self.pos[keys]
        else:
            pos = None

        if np.shape(self.vel)[0] > 0:
            vel = self.vel[keys]
        else:
            vel = None

        if np.shape(self.forces)[0] > 0:
            forces = self.forces[keys]
        else:
            forces = None

        if np.shape(self.charges)[0] > 0:
            charges = self.charges[keys]
        else:
            charges = None

        return AtomicSystem(atoms=atoms, pos=pos, vel=vel, forces=forces, charges=charges, cell=self.cell)
