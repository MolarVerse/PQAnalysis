import numpy as np

from numbers import Real

from ..cell import Cell
from ...types import Np1DNumberArray, Np2DNumberArray


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
# @check_atoms_has_mass
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
# @check_atoms_pos
# @check_atoms_has_mass
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
