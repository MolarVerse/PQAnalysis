"""
A module containing the AtomicSystem class and related functions.

...

Classes
-------
AtomicSystem
    A class for storing atomic systems.

Functions
---------
check_atoms_pos
    Decorator which checks that the number of atoms is equal to the number of positions.
check_atoms_has_mass
    Decorator which checks that all atoms have mass information.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import List, Any, Tuple, Iterable
from numbers import Real
from multimethod import multimethod

from .atom import Atom
from .cell import Cell
from ..types import Np2DNumberArray, Np1DNumberArray, Np2DIntArray, Np1DIntArray


def check_atoms_pos(func):
    """
    Decorator which checks that the number of atoms is equal to the number of positions.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Raises
    ------
    ValueError
        If the number of atoms is not equal to the number of positions.
    """
    def wrapper(*args, **kwargs):
        if args[0].pos.shape[0] != args[0].n_atoms:
            raise ValueError(
                "AtomicSystem contains a different number of atoms to positions.")
        return func(*args, **kwargs)
    return wrapper


def check_atoms_has_mass(func):
    """
    Decorator which checks that all atoms have mass information.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Raises
    ------
    ValueError
        If any atom does not have mass information.
    """
    def wrapper(*args, **kwargs):
        if not all([atom.mass is not None for atom in args[0].atoms]):
            raise ValueError(
                "AtomicSystem contains atoms without mass information.")
        return func(*args, **kwargs)
    return wrapper


class AtomicSystem:
    """
    A class for storing atomic systems.
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

    @property
    def PBC(self) -> bool:
        """
        Returns whether the system has periodic boundary conditions.

        Returns
        -------
        bool
            Whether the system has periodic boundary conditions. True if the system has a unit cell, False otherwise.
        """
        return self._cell != Cell()

    @property
    def n_atoms(self) -> int:
        """
        Returns the number of atoms in the system.

        Returns
        -------
        int
            The number of atoms in the system.
        """
        return len(self._atoms)

    @property
    @check_atoms_has_mass
    def atomic_masses(self) -> Np1DNumberArray:
        """
        Returns the masses of the atoms in the system.

        Returns
        -------
        Np1DNumberArray
            The masses of the atoms in the system.
        """
        return np.array([atom.mass for atom in self._atoms])

    @property
    def mass(self) -> Real:
        """
        Returns the total mass of the system.

        Returns
        -------
        Real
            The total mass of the system.
        """
        return np.sum(self.atomic_masses)

    @property
    @check_atoms_pos
    @check_atoms_has_mass
    def center_of_mass(self) -> Np1DNumberArray:
        """
        Returns the center of mass of the system.

        The decorated functions check that the number of atoms is equal to 
        the number of positions and that all atoms have mass information.

        Returns
        -------
        Np1DNumberArray
            The center of mass of the system.
        """
        # check if there are any atoms in the system otherwise self.mass would be 0
        if self.n_atoms == 0:
            return np.zeros(3)

        pos = self.cell.image(self.pos - self.pos[0]) + self.pos[0]

        pos = np.sum(
            pos * self.atomic_masses[:, None], axis=0) / self.mass

        pos = self.cell.image(pos)

        return pos

    @check_atoms_pos
    def nearest_neighbours(self, n: int = 1) -> Tuple[Np2DIntArray, Np2DNumberArray]:
        """
        Returns the n nearest neighbours of each atom in the system.

        Parameters
        ----------
        n : int, optional
            The number of nearest neighbours to return, by default 1

        Returns
        -------
        nearest_neighbours : Np2DIntArray
            The n nearest neighbours of each atom in the system.
        nearest_neighbours_distances : Np2DNumberArray
            The distances to the n nearest neighbours of each atom in the system.
        """
        nearest_neighbours = []
        nearest_neighbours_distances = []

        for atom_position in self.pos:
            delta_pos = self.pos - atom_position

            delta_pos = self.cell.image(delta_pos)

            distances = np.linalg.norm(delta_pos, axis=1)

            nearest_neighbours_atom = np.argsort(distances)[1:n+1]

            nearest_neighbours.append(nearest_neighbours_atom)
            nearest_neighbours_distances.append(
                distances[nearest_neighbours_atom])

        return np.array(nearest_neighbours), np.array(nearest_neighbours_distances)

    @property
    def combined_name(self) -> str:
        """
        Returns the combined name of the system.

        Therefore all atom.name are combined to a long string.

        Returns
        -------
        str
            The combined name of the system.
        """
        return ''.join([atom.name for atom in self.atoms])

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

        elif np.shape(self.pos) != np.shape(other.pos):
            return False

        elif np.shape(self.vel) != np.shape(other.vel):
            return False

        elif np.shape(self.forces) != np.shape(other.forces):
            return False

        elif np.shape(self.charges) != np.shape(other.charges):
            return False

        is_equal = True
        is_equal &= np.allclose(self.pos, other.pos)
        is_equal &= np.allclose(self.vel, other.vel)
        is_equal &= np.allclose(self.forces, other.forces)
        is_equal &= np.allclose(self.charges, other.charges)

        return is_equal

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
        indices = np.argwhere(np.array(self.atoms) == key)

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

    #######################
    #                     #
    # standard properties #
    #                     #
    #######################

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
