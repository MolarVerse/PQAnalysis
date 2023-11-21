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
from beartype.door import is_bearable

from ..atom import Atom, is_same_element_type
from ...types import Np1DIntArray


class _IndexingMixin:
    """
    A mixin class containing methods for indexing the atomic system.
    """

    def indices_from_atoms(self,
                           atoms: List[str] | List[Atom] | Np1DIntArray | None,
                           use_full_atom_info: bool = False
                           ) -> Np1DIntArray:
        """
        Returns the indices of the atoms with the given atom type names or atoms.

        Parameters
        ----------
        atoms : List[str] | List[Atom] | Np1DIntArray | None
            The atom type names or atoms to get the indices of. If None, all atoms are returned.

        Returns
        -------
        Np1DIntArray
            The indices of the atoms with the given atom type names or atoms.

        Raises
        ------
        ValueError
            If use_full_atom_info is True and atoms is not a List[Atom].
        """

        if atoms is None:
            return np.arange(self.n_atoms)

        elif not isinstance(atoms[0], Atom) and use_full_atom_info:
            raise ValueError(
                "use_full_atom_info can only be used with List[Atom]")

        elif isinstance(atoms[0], str):
            return self._indices_by_atom_type_names(atoms)

        elif isinstance(atoms[0], Atom) and use_full_atom_info:
            return self._indices_by_atom(atoms)

        elif isinstance(atoms[0], Atom) and not use_full_atom_info:
            return self._indices_by_element_types(atoms)

        # Note: here is is_bearable used instead of isinstance because
        #       isinstance(atoms, Np1DIntArray) returns False as
        #       Np1DIntArray is defined as a type alias and not a class.
        elif is_bearable(atoms, Np1DIntArray):
            return atoms

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
