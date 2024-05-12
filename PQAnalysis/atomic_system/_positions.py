"""
A module containing a Mixin class for related information to the positions of an atomic system.
"""

import numpy as np

from beartype.typing import Tuple

from PQAnalysis.core import distance
from PQAnalysis.topology import SelectionCompatible, Selection
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import (
    Np2DIntArray,
    Np2DNumberArray,
    Np1DIntArray,
    PositiveInt,
    Np1DNumberArray
)

from ._decorators import check_atoms_pos



class _PositionsMixin:

    """
    A mixin class containing methods for related 
    information to the positions of an atomic system.
    """

    @check_atoms_pos
    def _nearest_neighbours(
        self,
        n: PositiveInt = 1,
        indices: Np1DIntArray | None = None
    ) -> Tuple[Np2DIntArray,
        Np2DNumberArray]:
        """
        Returns the 'n' nearest neighbours of selected atoms in the system.

        If no indices are given, the 'n' nearest neighbours of all atoms are returned.
        With the parameter 'n' the number of nearest neighbours can be specified.

        Parameters
        ----------
        n : PositiveInt, optional
            The number of nearest neighbours to return, by default 1
        indices : Np1DIntArray, optional
            The indices of the atoms to get the nearest neighbours of,
            by default None (all atoms)

        Returns
        -------
        nearest_neighbours : Np2DIntArray
            The n nearest neighbour indices of each atom in the system.
        nearest_neighbours_distances : Np2DNumberArray
            The distances to the n nearest neighbours of each atom in the system.
        """

        indices = np.arange(self.n_atoms) if indices is None else indices

        nearest_neighbours = []
        nearest_neighbours_distances = []

        distances = distance(self.pos[indices], self.pos, self.cell)

        nearest_neighbours = np.argsort(distances, axis=-1)[:, 1:n + 1]
        nearest_neighbours_distances = np.take_along_axis(
            distances,
            nearest_neighbours,
            axis=-1
        )

        return nearest_neighbours, nearest_neighbours_distances

    @runtime_type_checking
    def nearest_neighbours(
        self,
        n: PositiveInt = 1,
        selection: SelectionCompatible = None,
        use_full_atom_info: bool = False
    ) -> Tuple[Np2DIntArray,
        Np2DNumberArray]:
        """
        Returns the n nearest neighbours of the given atoms in the system.

        If no selection of target atoms is given, the n nearest neighbours 
        of all atoms are returned. With the parameter 'n' the number of 
        nearest neighbours can be specified.

        Examples
        --------
        >>> import numpy as np
        >>> from PQAnalysis.atomic_system import AtomicSystem
        >>> from PQAnalysis.core import Atom

        >>> pos = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]])
        >>> atoms = [Atom('H'), Atom('H'), Atom('H')]
        >>> system = AtomicSystem(atoms=atoms, pos=pos)

        >>> system.nearest_neighbours()
        (array([[1], [0], [1]]), array([[0.8660254], [0.8660254], [0.8660254]]))

        >>> system.nearest_neighbours(n=2)
        (array([[1, 2],
               [0, 2],
               [1, 0]]), array([[0.8660254 , 1.73205081],
               [0.8660254 , 0.8660254 ],
               [0.8660254 , 1.73205081]]))

        >>> system.nearest_neighbours(selection=np.array([0]))
        (array([[1]]), array([[0.8660254]]))

        Parameters
        ----------
        n : PositiveInt, optional
            The number of nearest neighbours to return, by default 1
        selection : SelectionCompatible, optional
            Selection is either a selection object or any object that can 
            be initialized via 'Selection(selection)',
            default None (all atoms)
        use_full_atom_info : bool, optional
            If the full atom object should be used to match the 
            atoms or only the element type, 
            by default False

        Returns
        -------
        Tuple[Np2DIntArray, Np2DNumberArray]
            The n nearest neighbours of the given atoms in the system.
            The first array contains the indices of the nearest 
            neighbours and the second array contains the distances 
            to the nearest neighbours.
        """

        indices = Selection(selection
                            ).select(self.topology,
            use_full_atom_info)

        return self._nearest_neighbours(n=n, indices=indices)

    def image(self) -> None:
        """
        Images the positions of the system back into the cell.
        """
        self.pos = self.cell.image(self.pos)

    @runtime_type_checking
    def center(self, position: Np1DNumberArray, image: bool = True) -> None:
        """
        Center the positions of the system to a given position.

        Parameters
        ----------
        position : Np1DIntArray
            The position to recenter the system to.
        image : bool, optional
            If the system should be imaged back into the cell, by default True
        """
        self.pos -= position

        if image:
            self.image()
