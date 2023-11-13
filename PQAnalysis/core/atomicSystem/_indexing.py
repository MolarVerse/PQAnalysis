"""
A module containing a Mixin Class with different helpers for the indexing of atoms.

...

Classes
-------
_IndexingMixin
    A mixin class containing methods for indexing the atomic system.
"""

import numpy as np

from beartype.typing import List

from ..atom import Atom, is_same_element_type
from ...types import Np1DIntArray


class _IndexingMixin:
    """
    A mixin class containing methods for indexing the atomic system.
    """

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
            bool_array = np.array(
                [atom.name == name for atom in self.atoms])
            indices.append(np.argwhere(bool_array).flatten())
        return np.sort(np.concatenate(indices))

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
        return np.sort(np.concatenate(indices))

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
        return np.sort(np.concatenate(indices))
