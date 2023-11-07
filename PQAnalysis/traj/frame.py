"""
A module containing the Frame class.

...

Classes
-------
Frame
    A class for storing atomic systems with topology information.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import Any, List

from ..core.topology import Topology
from ..core.atomicSystem import AtomicSystem
from ..core.atom import Atom
from ..core.cell import Cell
from ..utils.mytypes import Numpy2DFloatArray, Numpy1DFloatArray


class Frame:
    """
    A class for storing atomic systems with topology information.

    ...

    Attributes
    ----------
    system : AtomicSystem
        The atomic system.
    topology : Topology
        The topology of the atomic system.
    """

    def __init__(self, system: AtomicSystem = AtomicSystem(), topology: Topology | None = None) -> None:
        """
        Initializes the Frame with the given parameters.

        Parameters
        ----------
        system : AtomicSystem, optional
            The atomic system, by default AtomicSystem()    
        topology : Topology, optional
            The topology of the atomic system, by default None
        """
        self.system = system
        self.topology = topology

    def compute_com_frame(self, group=None) -> Frame:
        """
        Computes a new Frame with the center of mass of the system or groups of atoms.  

        Parameters
        ----------
        group : int, optional
            group of atoms to compute the center of mass of, by default None (all atoms)

        Returns
        -------
        Frame
            A new Frame with the center of mass of the system or groups of atoms.

        Raises
        ------
        ValueError
            If the number of atoms in the selection is not a multiple of group.
        """
        if group is None:
            group = self.n_atoms

        elif self.n_atoms % group != 0:
            raise ValueError(
                'Number of atoms in selection is not a multiple of group.')

        pos = []
        names = []

        print(self.n_atoms)

        j = 0
        for i in range(0, self.n_atoms, group):
            atomic_system = AtomicSystem(
                atoms=self.atoms[i:i+group], pos=self.pos[i:i+group], cell=self.cell)

            print(atomic_system.center_of_mass)
            pos.append(atomic_system.center_of_mass)
            print(pos)
            names.append(atomic_system.combined_name)

            j += 1

        names = [Atom(name, use_guess_element=False) for name in names]

        print(pos)

        return Frame(AtomicSystem(pos=np.array(pos), atoms=names, cell=self.cell))

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the Frame is equal to another Frame.

        Parameters
        ----------
        other : Frame
            The other Frame to compare to.

        Returns
        -------
        bool
            Whether the Frame is equal to the other Frame.
        """
        if not isinstance(other, Frame):
            return False

        return self.system == other.system and self.topology == other.topology

    def __getitem__(self, key: int | slice) -> 'Frame':
        if self.topology is None:
            return Frame(system=self.system[key])
        else:
            return Frame(system=self.system[key], topology=self.topology[key])

    #########################
    #                       #
    #  Forwarded properties #
    #                       #
    #########################

    @property
    def PBC(self) -> bool:
        """
        Whether the system has periodic boundary conditions.

        Returns
        -------
        bool
            Whether the system has periodic boundary conditions.
        """
        return self.system.PBC

    @property
    def cell(self) -> Cell | None:
        """
        The unit cell of the system.

        Returns
        -------
        Cell | None
            The unit cell of the system.
        """
        return self.system.cell

    @cell.setter
    def cell(self, cell: Cell | None) -> Cell | None:
        """
        The unit cell of the system.

        Returns
        -------
        Cell | None
            The unit cell of the system.
        """
        self.system.cell = cell

    @property
    def n_atoms(self) -> int:
        """
        The number of atoms in the system.

        Returns
        -------
        int
            The number of atoms in the system.
        """
        return self.system.n_atoms

    @property
    def pos(self) -> Numpy2DFloatArray:
        """
        The positions of the atoms in the system.

        Returns
        -------
        Numpy2DFloatArray
            The positions of the atoms in the system.
        """
        return self.system.pos

    @property
    def atoms(self) -> List[Atom]:
        """
        The atoms in the system.

        Returns
        -------
        list[Atom]
            The atoms in the system.
        """
        return self.system.atoms
