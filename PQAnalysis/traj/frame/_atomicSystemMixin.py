"""
A module containing a Mixin Class with all methods and properties of an atomic system that are forwarded to the frame class.
"""

from beartype.typing import List

from PQAnalysis.core import Cell, Atom
from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray
from PQAnalysis.topology import Topology


class AtomicSystemMixin:
    @property
    def PBC(self) -> bool:
        """bool: Whether the system is periodic."""
        return self.system.PBC

    @property
    def cell(self) -> Cell:
        """Cell: The unit cell of the system."""
        return self.system.cell

    @cell.setter
    def cell(self, cell: Cell) -> None:
        self.system.cell = cell

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms in the system."""
        return self.system.n_atoms

    @property
    def pos(self) -> Np2DNumberArray:
        """Np2DNumberArray: The positions of the atoms in the system."""
        return self.system.pos

    @pos.setter
    def pos(self, pos: Np2DNumberArray) -> None:
        self.system.pos = pos

    @property
    def vel(self) -> Np2DNumberArray:
        """Np2DNumberArray: The velocities of the atoms in the system."""
        return self.system.vel

    @vel.setter
    def vel(self, vel: Np2DNumberArray) -> None:
        self.system.vel = vel

    @property
    def forces(self) -> Np2DNumberArray:
        """Np2DNumberArray: The forces of the atoms in the system."""
        return self.system.forces

    @forces.setter
    def forces(self, forces: Np2DNumberArray) -> None:
        self.system.forces = forces

    @property
    def charges(self) -> Np1DNumberArray:
        """Np1DNumberArray: The charges of the atoms in the system."""
        return self.system.charges

    @charges.setter
    def charges(self, charges: Np1DNumberArray) -> None:
        self.system.charges = charges

    @property
    def atoms(self) -> List[Atom]:
        """List[Atom]: The atoms in the system."""
        return self.system.atoms

    @property
    def topology(self) -> Topology:
        """Topology: The topology of the system."""
        return self.system.topology

    @topology.setter
    def topology(self, topology: Topology) -> None:
        self.system.topology = topology

    @property
    def center_of_mass(self) -> Np1DNumberArray:
        """Np1DNumberArray: The center of mass of the system."""
        return self.system.center_of_mass

    def image(self) -> None:
        """
        Images the positions of the atoms in the system to the unit cell
        """
        self.system.image()

    def center(self, position: Np1DNumberArray, image: bool = True) -> None:
        """
        Centers the system around a given position.

        Parameters
        ----------
        position : Np1DNumberArray, optional
            The position around which the system should be centered, by default None
        image : bool, optional
            Whether the positions should be imaged back into the cell, by default True
        """
        self.system.center(position, image=image)
