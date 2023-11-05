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

import numpy as np
import numpy.typing as npt

from beartype import beartype
from beartype.typing import List

from .atom import Atom
from .cell import Cell
from ..utils.mytypes import Numpy2DFloatArray, Numpy1DFloatArray


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


@beartype
class AtomicSystem:
    """
    A class for storing atomic systems.

    ...

    Attributes
    ----------
    atoms : list
        A list of Atom objects, by default []
    pos : Numpy2DFloatArray
        A 2d numpy.ndarray containing the positions of the atoms, by default np.zeros((0, 3)).
    vel : Numpy2DFloatArray
        A 2d numpy.ndarray containing the velocities of the atoms, by default np.zeros((0, 3)).
    forces : Numpy2DFloatArray
        A 2d numpy.ndarray containing the forces on the atoms, by default np.zeros((0, 3)).
    charges : Numpy1DFloatArray
        A 1d numpy.ndarray containing the charges of the atoms, by default np.zeros(0).
    cell : Cell
        The unit cell of the system, by default None.
    """

    def __init__(self,
                 atoms: List[Atom] = None,
                 pos: Numpy2DFloatArray = np.zeros((0, 3)),
                 vel: Numpy2DFloatArray = np.zeros((0, 3)),
                 forces: Numpy2DFloatArray = np.zeros((0, 3)),
                 charges: Numpy1DFloatArray = np.zeros(0),
                 cell: Cell | None = None
                 ):
        """
        Initializes the AtomicSystem with the given parameters.

        Parameters
        ----------
        atoms : List[Atom], optional
            A list of Atom objects, by default []
        pos : Numpy2DFloatArray, optional
            A 2d numpy.ndarray containing the positions of the atoms, by default np.zeros((0, 3)).
        vel : Numpy2DFloatArray, optional
            A 2d numpy.ndarray containing the velocities of the atoms, by default np.zeros((0, 3)).
        forces : Numpy2DFloatArray, optional
            A 2d numpy.ndarray containing the forces on the atoms, by default np.zeros((0, 3)).
        charges : Numpy1DFloatArray, optional
            A 1d numpy.ndarray containing the charges of the atoms, by default np.zeros(0).
        cell : Cell | None, optional
            The unit cell of the system, by default None.
        """
        if atoms is None:
            atoms = []

        self._atoms = atoms
        self._pos = pos
        self._vel = vel
        self._forces = forces
        self._charges = charges
        self._cell = cell

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
    def pos(self) -> Numpy2DFloatArray:
        """
        Returns the positions of the atoms in the system.

        Returns
        -------
        Numpy2DFloatArray
            The positions of the atoms in the system.
        """
        return self._pos

    @property
    def vel(self) -> Numpy2DFloatArray:
        """
        Returns the velocities of the atoms in the system.

        Returns
        -------
        Numpy2DFloatArray
            The velocities of the atoms in the system.
        """
        return self._vel

    @property
    def forces(self) -> Numpy2DFloatArray:
        """
        Returns the forces on the atoms in the system.

        Returns
        -------
        Numpy2DFloatArray
            The forces on the atoms in the system.
        """
        return self._forces

    @property
    def charges(self) -> Numpy1DFloatArray:
        """
        Returns the charges of the atoms in the system.

        Returns
        -------
        Numpy1DFloatArray
            The charges of the atoms in the system.
        """
        return self._charges

    @property
    def cell(self) -> Cell | None:
        """
        Returns the unit cell of the system.

        Returns
        -------
        Cell | None
            The unit cell of the system.
        """
        return self._cell

    @property
    def PBC(self) -> bool:
        """
        Returns whether the system has periodic boundary conditions.

        Returns
        -------
        bool
            Whether the system has periodic boundary conditions. True if the system has a unit cell, False otherwise.
        """
        return self._cell is not None

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
    def atomic_masses(self) -> Numpy1DFloatArray:
        """
        Returns the masses of the atoms in the system.

        Returns
        -------
        Numpy1DFloatArray
            The masses of the atoms in the system.
        """
        return np.array([atom.mass for atom in self._atoms])

    @property
    def mass(self) -> float:
        """
        Returns the total mass of the system.

        Returns
        -------
        float
            The total mass of the system.
        """
        return np.sum(self.atomic_masses)

    @property
    @check_atoms_pos
    @check_atoms_has_mass
    def center_of_mass(self) -> Numpy1DFloatArray:
        """
        Returns the center of mass of the system.

        The decorated functions check that the number of atoms is equal to 
        the number of positions and that all atoms have mass information.

        Returns
        -------
        Numpy1DFloatArray
            The center of mass of the system.
        """
        # check if there are any atoms in the system otherwise self.mass would be 0
        if self.n_atoms == 0:
            return np.zeros(3)

        # if sell is not None image the positions relative to the first atom
        if self.cell is not None:
            pos = self.cell.image(self.pos - self.pos[0]) + self.pos[0]
        else:
            pos = self.pos

        pos = np.sum(
            pos * self.atomic_masses[:, None], axis=0) / self.mass

        # back imaging
        if self.cell is not None:
            pos = self.cell.image(pos)

        return pos

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
