"""
A module containing a Mixin Class with the standard properties
of an atomic system (i.e. standard getter and setter methods).
"""

from PQAnalysis.core import Atoms, Cell
from PQAnalysis.topology import Topology
from PQAnalysis.type_checking import runtime_type_checking_setter
from PQAnalysis.types import (
    Np1DNumberArray,
    Np2DNumberArray,
    Np3x3NumberArray,
    Real,
)

from ._decorators import check_atom_number_setters



class _StandardPropertiesMixin:

    """
    A mixin class containing the standard properties of an atomic 
    system (i.e. standard getter and setter methods).
    """

    @property
    def atoms(self) -> Atoms:
        """Atoms: The atoms in the system."""
        return self.topology.atoms

    @property
    def topology(self) -> Topology:
        """
        Topology: The topology of the system.

        In order to set the topology of the system, the number of atoms
        in the topology has to be equal to the number of atoms in the system.
        """
        return self._topology

    @topology.setter
    @runtime_type_checking_setter
    def topology(self, topology: Topology) -> None:
        if topology.n_atoms != self.n_atoms:
            raise ValueError(
                "The number of atoms already found in the AtomicSystem "
                "object have to be equal to the number of atoms in the new topology"
            )

        self._topology = topology

    @property
    def pos(self) -> Np2DNumberArray:
        """
        Np2DNumberArray: The positions of the atoms in the system.

        In order to set the positions of the atoms in the system, 
        the number of atoms in the system has to be equal to the 
        number of positions.
        """
        return self._pos

    @pos.setter
    @runtime_type_checking_setter
    @check_atom_number_setters
    def pos(self, pos: Np2DNumberArray) -> None:
        self._pos = pos

    def set_pos_no_checks(self, pos: Np2DNumberArray) -> None:
        """
        Set the positions of the atoms in the system without any checks.

        Parameters
        ----------
        pos : Np2DNumberArray
            The positions of the atoms in the system.
        """
        self._pos = pos

    @property
    def has_pos(self) -> bool:
        """bool: A boolean indicating if the system has positions for all atoms."""
        return len(self.pos) == self.n_atoms and self.n_atoms != 0

    @property
    def vel(self) -> Np2DNumberArray:
        """
        Np2DNumberArray: The velocities of the atoms in the system.

        In order to set the velocities of the atoms in the system, 
        the number of atoms in the system has to be equal to
        the number of velocities.
        """
        return self._vel

    @vel.setter
    @runtime_type_checking_setter
    @check_atom_number_setters
    def vel(self, vel: Np2DNumberArray) -> None:
        self._vel = vel

    def set_vel_no_checks(self, vel: Np2DNumberArray) -> None:
        """
        Set the velocities of the atoms in the system without any checks.

        Parameters
        ----------
        vel : Np2DNumberArray
            The velocities of the atoms in the system.
        """
        self._vel = vel

    @property
    def has_vel(self) -> bool:
        """bool: A boolean indicating if the system has velocities for all atoms."""
        return len(self.vel) == self.n_atoms and self.n_atoms != 0

    @property
    def forces(self) -> Np2DNumberArray:
        """
        Np2DNumberArray: The forces acting on the atoms in the system.

        In order to set the forces acting on the atoms in the system,
        the number of atoms in the system has to be equal to the number
        of forces.
        """
        return self._forces

    @forces.setter
    @runtime_type_checking_setter
    @check_atom_number_setters
    def forces(self, forces: Np2DNumberArray) -> None:
        self._forces = forces

    def set_forces_no_checks(self, forces: Np2DNumberArray) -> None:
        """
        Set the forces acting on the atoms in the system without any checks.

        Parameters
        ----------
        forces : Np2DNumberArray
            The forces acting on the atoms in the system.
        """
        self._forces = forces

    @property
    def has_forces(self) -> bool:
        """bool: A boolean indicating if the system has forces for all atoms."""
        return len(self.forces) == self.n_atoms and self.n_atoms != 0

    @property
    def charges(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The charges of the atoms in the system.

        In order to set the charges of the atoms in the system,
        the number of atoms in the system has to be equal to the
        number of charges.
        """
        return self._charges

    @charges.setter
    @runtime_type_checking_setter
    @check_atom_number_setters
    def charges(self, charges: Np1DNumberArray) -> None:
        self._charges = charges

    def set_charges_no_checks(self, charges: Np1DNumberArray) -> None:
        """
        Set the charges of the atoms in the system without any checks.

        Parameters
        ----------
        charges : Np1DNumberArray
            The charges of the atoms in the system.
        """
        self._charges = charges

    @property
    def has_charges(self) -> bool:
        """bool: A boolean indicating if the system has charges for all atoms."""
        return len(self.charges) == self.n_atoms and self.n_atoms != 0

    @property
    def cell(self) -> Cell:
        """Cell: The unit cell of the system."""
        return self._cell

    @cell.setter
    @runtime_type_checking_setter
    def cell(self, cell: Cell) -> None:
        self._cell = cell

    @property
    def energy(self) -> Real | None:
        """float: The energy of the system."""
        return self._energy

    @energy.setter
    @runtime_type_checking_setter
    def energy(self, energy: Real) -> None:
        self._energy = energy

    @property
    def has_energy(self) -> bool:
        """bool: A boolean indicating if the system has an energy."""
        return self._energy is not None

    @property
    def stress(self) -> Np3x3NumberArray | None:
        """Np2DNumberArray: The stress tensor of the system."""
        return self._stress

    @stress.setter
    @runtime_type_checking_setter
    def stress(self, stress: Np3x3NumberArray) -> None:
        self._stress = stress

    @property
    def has_stress(self) -> bool:
        """bool: A boolean indicating if the system has a stress tensor."""
        return self._stress is not None

    @property
    def virial(self) -> Np3x3NumberArray | None:
        """Np2DNumberArray: The virial tensor of the system."""
        return self._virial

    @virial.setter
    @runtime_type_checking_setter
    def virial(self, virial: Np3x3NumberArray) -> None:
        self._virial = virial

    @property
    def has_virial(self) -> bool:
        """bool: A boolean indicating if the system has a virial tensor."""
        return self._virial is not None
