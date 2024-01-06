"""
A module containing the AtomicSystem class

...

Classes
-------
AtomicSystem
    A class for storing atomic systems.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import Any

from ._properties import _PropertiesMixin
from ._standardProperties import _StandardPropertiesMixin
from ._positions import _PositionsMixin

from ..atom import Atom, Atoms
from ..cell import Cell
from ...types import Np2DNumberArray, Np1DNumberArray, Np1DIntArray
from ...topology import Topology


class AtomicSystem(_PropertiesMixin, _StandardPropertiesMixin, _PositionsMixin):
    """
    A class for storing atomic systems.

    Inherits from the Mixins: _PropertiesMixin, _StandardPropertiesMixin, _IndexingMixin, _PositionsMixin
        - The _StandardPropertiesMixin contains the standard properties of an atomic system (i.e. standard getter and setter methods).
        - The _PropertiesMixin contains special properties derived from the standard properties
        - The _IndexingMixin contains methods for indexing the atomic system
        - The _PositionsMixin contains methods for computing properties based on the positions of the atoms


    Examples
    --------
    >>> AtomicSystem(atoms=[Atom('C1'), Atom('C2')], pos=np.array([[0, 0, 0], [1, 0, 0]]))

    >>> AtomicSystem()

    >>> AtomicSystem(topology=Topology(atoms=[Atom('C1'), Atom('C2')]), pos=np.array([[0, 0, 0], [1, 0, 0]])
    """

    def __init__(self,
                 atoms: Atoms | None = None,
                 pos: Np2DNumberArray | None = None,
                 vel: Np2DNumberArray | None = None,
                 forces: Np2DNumberArray | None = None,
                 charges: Np1DNumberArray | None = None,
                 topology: Topology | None = None,
                 cell: Cell = Cell()
                 ) -> None:
        """
        Initializes the AtomicSystem with the given parameters.

        For the initialization of an AtomicSystem all parameters are optional. 
        If no value is given for a parameter, the default value is used which 
        is an empty list for atoms, an empty numpy.ndarray for pos, vel, forces
        and charges, a Topology() object for topology and a Cell() object for cell.

        If the shapes or lengths of the given parameters are not consistent, this will 
        only raise an error when a property or method is called that requires the
        given parameter. This is done to allow for the creation of an AtomicSystem
        with only a subset of the properties.

        One important restriction is that atoms and topology are mutually exclusive, 
        i.e. if atoms is given, topology cannot be given and vice versa (this is
        because the topology is derived from the atoms - if given).

        Parameters
        ----------
        atoms : Atoms, optional
            A list of Atom objects, by default []
        pos : Np2DNumberArray, optional
            A 2d numpy.ndarray containing the positions of the atoms, by default np.zeros((0, 3)).
        vel : Np2DNumberArray, optional
            A 2d numpy.ndarray containing the velocities of the atoms, by default np.zeros((0, 3)).
        forces : Np2DNumberArray, optional
            A 2d numpy.ndarray containing the forces on the atoms, by default np.zeros((0, 3)).
        charges : Np1DNumberArray, optional
            A 1d numpy.ndarray containing the charges of the atoms, by default np.zeros(0).
        topology : Topology, optional
            The topology of the system, by default Topology()
        cell : Cell, optional
            The unit cell of the system. Defaults to a Cell with no periodic boundary conditions, by default Cell()

        Raises
        ------
        ValueError
            If both atoms and topology are given.
        """
        if atoms is None and topology is None:
            topology = Topology()
        elif topology is None:
            topology = Topology(atoms=atoms)
        elif topology is not None and atoms is not None:
            raise ValueError(
                "Cannot initialize AtomicSystem with both atoms and topology arguments - they are mutually exclusive.")

        if pos is None:
            pos = np.zeros((0, 3))

        if vel is None:
            vel = np.zeros((0, 3))

        if forces is None:
            forces = np.zeros((0, 3))

        if charges is None:
            charges = np.zeros(0)

        self._topology = topology
        self._pos = pos
        self._vel = vel
        self._forces = forces
        self._charges = charges
        self._cell = cell

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the AtomicSystem is equal to another AtomicSystem.

        Parameters
        ----------
        other : AtomicSystem
            The other AtomicSystem to compare to.

        Returns
        -------
        bool
            Whether the AtomicSystem is equal to the other AtomicSystem.
        """
        if not isinstance(other, AtomicSystem):
            return False

        is_equal = True

        if self.n_atoms != other.n_atoms:
            return False

        is_equal &= self.topology == other.topology
        is_equal &= self.cell == other.cell
        is_equal &= np.allclose(self.pos, other.pos)
        is_equal &= np.allclose(self.vel, other.vel)
        is_equal &= np.allclose(self.forces, other.forces)
        is_equal &= np.allclose(self.charges, other.charges)

        return is_equal

    def __getitem__(self, key: Atom | int | slice | Np1DIntArray) -> AtomicSystem:
        """
        Returns a new AtomicSystem with the given key.

        Examples
        --------
        >>> system1 = AtomicSystem(atoms=[Atom('C1'), Atom('C2')], pos=np.array([[0, 0, 0], [1, 0, 0]]))
        >>> system1[0]
        AtomicSystem(atoms=[Atom('C1')], pos=np.array([[0, 0, 0]]))
        >>> system1[0:2]
        AtomicSystem(atoms=[Atom('C1'), Atom('C2')], pos=np.array([[0, 0, 0], [1, 0, 0]]))
        >>> system1[np.array([0, 1])]
        AtomicSystem(atoms=[Atom('C1'), Atom('C2')], pos=np.array([[0, 0, 0], [1, 0, 0]]))

        Parameters
        ----------
        key : Atom | int | slice | Np1DIntArray
            The key to get the new AtomicSystem with.

        Returns
        -------
        AtomicSystem
            The new AtomicSystem with the given key.
        """

        if isinstance(key, Atom):
            key = np.argwhere(np.array(self.atoms) == key).flatten()

        if isinstance(key, int):
            keys = np.array([key])
        elif isinstance(key, slice):
            keys = np.array(range(self.n_atoms)[key])
        else:
            keys = np.array(key)

        if np.shape(self.pos)[0] > 0:
            pos = self.pos[keys]
        else:
            pos = None

        if np.shape(self.vel)[0] > 0:
            vel = self.vel[keys]
        else:
            vel = None

        if np.shape(self.forces)[0] > 0:
            forces = self.forces[keys]
        else:
            forces = None

        if np.shape(self.charges)[0] > 0:
            charges = self.charges[keys]
        else:
            charges = None

        return AtomicSystem(pos=pos, vel=vel, forces=forces, charges=charges, cell=self.cell, topology=self.topology[keys])

    def __str__(self) -> str:
        """
        Returns the string representation of the AtomicSystem.

        Returns
        -------
        str
            The string representation of the AtomicSystem.
        """
        return f"AtomicSystem(topology=({self.topology}), cell=({self.cell}))"

    def __repr__(self) -> str:
        """
        Returns the string representation of the AtomicSystem.

        Returns
        -------
        str
            The string representation of the AtomicSystem.
        """
        return self.__str__()
