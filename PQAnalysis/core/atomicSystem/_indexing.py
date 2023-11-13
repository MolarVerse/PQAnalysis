import numpy as np

from beartype.typing import List

from ..atom import Atom, is_same_element_type
from ...types import Np1DIntArray


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


def _indices_by_atom_type_names(self, names: List[str]) -> Np1DIntArray:
    """
    Returns the indices of the atoms with the given atom type names.

    Parameters
    ----------
    names : List[str]
        The names of the atom types to get the indices of.

    Returns
    -------
    Np1DIntArray
        The indices of the atoms with the given atom type names.
    """
    indices = []
    for name in names:
        indices.append(np.argwhere(np.array(self.atoms) == name).flatten())
    return np.concatenate(indices)


def _indices_by_atom(self, atoms: List[Atom]) -> Np1DIntArray:
    """
    Returns the indices of the given atoms.

    Parameters
    ----------
    atoms : List[Atom]
        The atoms to get the indices of.

    Returns
    -------
    Np1DIntArray
        The indices of the given atoms.
    """
    indices = []
    for atom in atoms:
        indices.append(np.argwhere(np.array(self.atoms) == atom).flatten())
    return np.concatenate(indices)


def _indices_by_element_types(self, elements: List[Atom]) -> Np1DIntArray:
    """
    Returns the indices of the atoms with the given element types.

    Parameters
    ----------
    elements : List[Atom]
        The element types to get the indices of.

    Returns
    -------
    Np1DIntArray
        The indices of the atoms with the given element types.
    """
    indices = []
    for element in elements:
        bool_indices = np.array(
            [is_same_element_type(atom, element) for atom in self.atoms])
        indices.append(np.argwhere(bool_indices).flatten())
    return np.concatenate(indices)
