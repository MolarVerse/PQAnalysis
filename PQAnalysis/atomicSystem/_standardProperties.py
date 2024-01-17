"""
A module containing a Mixin Class with the standard properties of an atomic system (i.e. standard getter and setter methods).
"""

from PQAnalysis.core import Atoms, Cell
from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray
from PQAnalysis.topology import Topology


class _StandardPropertiesMixin:
    """
    A mixin class containing the standard properties of an atomic system (i.e. standard getter and setter methods).
    """
    @property
    def atoms(self) -> Atoms:
        """Atoms: The atoms in the system."""
        return self.topology.atoms

    @property
    def topology(self) -> Topology:
        """
        Topology: The topology of the system.

        In order to set the topology of the system, the number of atoms in the topology has to be equal to the number of atoms in the system.
        """
        return self._topology

    @topology.setter
    def topology(self, topology: Topology) -> None:
        if topology.n_atoms != self.n_atoms:
            raise ValueError(
                "The number of atoms already found in the AtomicSystem object have to be equal to the number of atoms in the new topology")

        self._topology = topology

    @property
    def pos(self) -> Np2DNumberArray:
        """Np2DNumberArray: The positions of the atoms in the system."""
        return self._pos

    @property
    def vel(self) -> Np2DNumberArray:
        """Np2DNumberArray: The velocities of the atoms in the system."""
        return self._vel

    @property
    def forces(self) -> Np2DNumberArray:
        """Np2DNumberArray: The forces acting on the atoms in the system."""
        return self._forces

    @property
    def charges(self) -> Np1DNumberArray:
        """Np1DNumberArray: The charges of the atoms in the system."""
        return self._charges

    @property
    def cell(self) -> Cell:
        """Cell: The unit cell of the system."""
        return self._cell

    @cell.setter
    def cell(self, cell: Cell) -> None:
        self._cell = cell
