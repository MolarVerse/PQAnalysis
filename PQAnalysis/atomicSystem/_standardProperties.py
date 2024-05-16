"""
A module containing a Mixin Class with the standard properties of an atomic system (i.e. standard getter and setter methods).
"""

from ._decorators import check_atom_number_setters
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
        """
        Np2DNumberArray: The positions of the atoms in the system.

        In order to set the positions of the atoms in the system, the number of atoms in the system has to be equal to the number of positions.
        """
        return self._pos

    @pos.setter
    @check_atom_number_setters
    def pos(self, pos: Np2DNumberArray) -> None:
        self._pos = pos

    @property
    def has_pos(self) -> bool:
        """bool: A boolean indicating if the system has positions for all atoms."""
        return len(self.pos) == self.n_atoms

    @property
    def vel(self) -> Np2DNumberArray:
        """
        Np2DNumberArray: The velocities of the atoms in the system.

        In order to set the velocities of the atoms in the system, the number of atoms in the system has to be equal to the number of velocities.
        """
        return self._vel

    @vel.setter
    @check_atom_number_setters
    def vel(self, vel: Np2DNumberArray) -> None:
        self._vel = vel

    @property
    def has_vel(self) -> bool:
        """bool: A boolean indicating if the system has velocities for all atoms."""
        return len(self.vel) == self.n_atoms

    @property
    def forces(self) -> Np2DNumberArray:
        """
        Np2DNumberArray: The forces acting on the atoms in the system.

        In order to set the forces acting on the atoms in the system, the number of atoms in the system has to be equal to the number of forces.
        """
        return self._forces

    @forces.setter
    @check_atom_number_setters
    def forces(self, forces: Np2DNumberArray) -> None:
        self._forces = forces

    @property
    def has_forces(self) -> bool:
        """bool: A boolean indicating if the system has forces for all atoms."""
        return len(self.forces) == self.n_atoms

    @property
    def charges(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The charges of the atoms in the system.

        In order to set the charges of the atoms in the system, the number of atoms in the system has to be equal to the number of charges.
        """
        return self._charges

    @charges.setter
    @check_atom_number_setters
    def charges(self, charges: Np1DNumberArray) -> None:
        self._charges = charges

    @property
    def has_charges(self) -> bool:
        """bool: A boolean indicating if the system has charges for all atoms."""
        return len(self.charges) == self.n_atoms

    @property
    def cell(self) -> Cell:
        """Cell: The unit cell of the system."""
        return self._cell

    @cell.setter
    def cell(self, cell: Cell) -> None:
        self._cell = cell

    @property
    def energy(self) -> float:
        """float: The energy of the system."""
        return self._energy

    @energy.setter
    def energy(self, energy: float) -> None:
        self._energy = energy

    @property
    def has_energy(self) -> bool:
        """bool: A boolean indicating if the system has an energy."""
        return self._energy is not None

    @property
    def stress(self) -> Np2DNumberArray:
        """Np2DNumberArray: The stress tensor of the system."""
        return self._stress

    @stress.setter
    def stress(self, stress: Np2DNumberArray) -> None:
        self._stress = stress

    @property
    def has_stress(self) -> bool:
        """bool: A boolean indicating if the system has a stress tensor."""
        return self._stress is not None

    @property
    def virial(self) -> Np2DNumberArray:
        """Np2DNumberArray: The virial tensor of the system."""
        return self._virial

    @virial.setter
    def virial(self, virial: Np2DNumberArray) -> None:
        self._virial = virial

    @property
    def has_virial(self) -> bool:
        """bool: A boolean indicating if the system has a virial tensor."""
        return self._virial is not None
