import numpy as np

from multimethod import multimethod
from beartype.typing import List, Tuple

# from ._decorators import check_atoms_pos
from ..atom import Atom
from ...types import Np2DIntArray, Np2DNumberArray, Np1DIntArray


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


# @check_atoms_pos
@multimethod
def nearest_neighbours(self,
                       n: int = 1,
                       atom_types: List[Atom] | List[str] | None = None,
                       element_types: List[Atom] | None = None,
                       ) -> Tuple[Np2DIntArray, Np2DNumberArray]:

    if atom_types is not None and element_types is not None:
        raise ValueError(
            "Cannot use both atom_types and element_types they are mutual exclusive.")

    elif atom_types is None and element_types is None:
        indices = None

    elif atom_types is not None and isinstance(atom_types[0], str):
        indices = self._indices_by_atom_type_names(atom_types)

    elif atom_types is not None and isinstance(atom_types[0], Atom):
        indices = self._indices_by_atom(atom_types)

    elif element_types is not None:
        indices = self._indices_by_element_types(element_types)

    return self.nearest_neighbours(n=n, indices=indices)


# @check_atoms_pos
@multimethod
def nearest_neighbours(self,
                       n: int = 1,
                       indices: Np1DIntArray | None = None
                       ) -> Tuple[Np2DIntArray, Np2DNumberArray]:
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
    indices : Np1DIntArray, optional
        The indices of the atoms to get the nearest neighbours of, by default None (all atoms)
    """

    if indices is None:
        indices = np.arange(self.n_atoms)

    nearest_neighbours = []
    nearest_neighbours_distances = []

    for atom_position in self.pos[indices]:
        delta_pos = self.pos - atom_position

        delta_pos = self.cell.image(delta_pos)

        distances = np.linalg.norm(delta_pos, axis=1)

        nearest_neighbours_atom = np.argsort(distances)[1:n+1]

        nearest_neighbours.append(nearest_neighbours_atom)
        nearest_neighbours_distances.append(
            distances[nearest_neighbours_atom])

    return np.array(nearest_neighbours), np.array(nearest_neighbours_distances)
