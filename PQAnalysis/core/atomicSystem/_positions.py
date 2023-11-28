"""
A module containing a Mixin class for related information to the positions of an atomic system.

...

Classes
-------
_PositionsMixin
    A mixin class containing methods for related information to the positions of an atomic system.
"""

import numpy as np

from multimethod import multimethod
from beartype.typing import List, Tuple
from beartype.door import is_bearable

from ._decorators import check_atoms_pos

from .. import Atom, distance
from ...types import Np2DIntArray, Np2DNumberArray, Np1DIntArray, Np1DNumberArray


class _PositionsMixin:
    """
    A mixin class containing methods for related information to the positions of an atomic system.
    """

    @check_atoms_pos
    def _nearest_neighbours(self,
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
            distances = distance(atom_position, self.pos, self.cell)

            nearest_neighbours_atom = np.argsort(distances)[1:n+1]

            nearest_neighbours.append(nearest_neighbours_atom)
            nearest_neighbours_distances.append(
                distances[nearest_neighbours_atom])

        return np.array(nearest_neighbours), np.array(nearest_neighbours_distances)

    def nearest_neighbours(self,
                           n: int = 1,
                           atoms: List[Atom] | List[str] | Np1DIntArray | None = None,
                           use_full_atom_info: bool = False
                           ) -> Tuple[Np2DIntArray, Np2DNumberArray]:
        """
        Returns the n nearest neighbours of the given atoms in the system.

        It is possible to specify the atoms by their element type, by their name or by their index 
        or by the full Atom object. The parameter n specifies the number of closest nearest neighbours to return.
        With the parameter use_full_atom_info it is possible to specify if the atoms should be searched for as 
        only their element types or as their full atom object. (default: element types)

        for example:
            The object Atom('H1', 1) is equal to atom Atom('H') and Atom(1) if use_full_atom_info is False.
            The object Atom('H1', 1) is not equal to atom Atom('H') and Atom(1) if use_full_atom_info is True, 
            because also the atom_type name is compared.

        Parameters
        ----------
        n : int, optional
            The number of nearest neighbours to return, by default 1
        atoms : List[Atom] | List[str] | Np1DIntArray | None, optional
            The atoms to get the nearest neighbours of, by default None (all atoms)
        use_full_atom_info : bool, optional
            If the full atom object should be used to match the atoms, by default False

        Returns
        -------
        Tuple[Np2DIntArray, Np2DNumberArray]
            The n nearest neighbours of the given atoms in the system.
            The first array contains the indices of the nearest neighbours and the second array contains the distances to the nearest neighbours.
        """

        indices = self.indices_from_atoms(atoms, use_full_atom_info)

        return self._nearest_neighbours(n=n, indices=indices)
