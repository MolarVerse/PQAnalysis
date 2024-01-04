"""
A module containing a Mixin class for related information to the positions of an atomic system.

...

Classes
-------
_PositionsMixin
    A mixin class containing methods for related information to the positions of an atomic system.
"""

import numpy as np

from beartype.typing import Tuple

from ._decorators import check_atoms_pos
from .. import distance
from ...types import Np2DIntArray, Np2DNumberArray, Np1DIntArray, PositiveInt
from ...topology.selection import SelectionCompatible, Selection


class _PositionsMixin:
    """
    A mixin class containing methods for related information to the positions of an atomic system.
    """

    @check_atoms_pos
    def _nearest_neighbours(self,
                            n: PositiveInt = 1,
                            indices: Np1DIntArray | None = None
                            ) -> Tuple[Np2DIntArray, Np2DNumberArray]:
        """
        Returns the n nearest neighbours of each atom in the system.

        If no indices are given, the nearest neighbours of all atoms are returned.

        Parameters
        ----------
        n : PositiveInt, optional
            The number of nearest neighbours to return, by default 1
        indices : Np1DIntArray, optional
            The indices of the atoms to get the nearest neighbours of, by default None (all atoms)

        Returns
        -------
        nearest_neighbours : Np2DIntArray
            The n nearest neighbours of each atom in the system.
        nearest_neighbours_distances : Np2DNumberArray
            The distances to the n nearest neighbours of each atom in the system.
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
                           n: PositiveInt = 1,
                           selection: SelectionCompatible = None,
                           use_full_atom_info: bool = False
                           ) -> Tuple[Np2DIntArray, Np2DNumberArray]:
        """
        Returns the n nearest neighbours of the given atoms in the system.

        Examples
        --------
        >>> import numpy as np
        >>> from PQAnalysis import AtomicSystem
        >>> pos = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]])
        >>> system = AtomicSystem(pos=pos)
        >>> system.nearest_neighbours()
        (np.array([[1], [0], [1]]), np.array([[0.8660254], [0.8660254], [0.8660254]]))
        >>> system.nearest_neighbours(n=2)
        (np.array([[1, 2], [0, 2], [1, 0]]), np.array([[0.8660254, 1.73205081], [0.8660254, 1.73205081], [0.8660254, 1.73205081]]))
        >>> system.nearest_neighbours(selection=np.array([0]))
        (np.array([[1]]), np.array([[0.8660254]]))

        Parameters
        ----------
        n : PositiveInt, optional
            The number of nearest neighbours to return, by default 1
        selection : SelectionCompatible, optional
            Selection is either a selection object or any object that can be initialized via 'Selection(selection)', 
            default None (all atoms)
        use_full_atom_info : bool, optional
            If the full atom object should be used to match the atoms or only the element type, by default False

        Returns
        -------
        Tuple[Np2DIntArray, Np2DNumberArray]
            The n nearest neighbours of the given atoms in the system.
            The first array contains the indices of the nearest neighbours and the second array contains 
            the distances to the nearest neighbours.
        """

        indices = Selection(selection).select(
            self.topology, use_full_atom_info)

        return self._nearest_neighbours(n=n, indices=indices)
