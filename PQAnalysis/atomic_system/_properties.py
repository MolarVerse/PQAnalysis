"""
A module containing a Mixin Class with different properties derived
from the standard properties of an atomic system.
"""
from numbers import Real

import numpy as np

from beartype.typing import List

from PQAnalysis.core import Cell, CustomElement
from PQAnalysis.types import Np1DNumberArray

from ._decorators import check_atoms_has_mass, check_atoms_pos
from .exceptions import AtomicSystemError



class _PropertiesMixin:

    """
    A mixin class containing properties derived from the standard properties of an atomic system.
    """

    @property
    def pbc(self) -> bool:
        """bool: Whether the system has periodic boundary conditions."""
        return self._cell != Cell()

    @property
    def n_atoms(self) -> int:
        """
        int: The number of atoms in the system.

        Raises
        ------
        AtomicSystemError
            If the number of atoms in the topology, positions, velocities, 
            forces and charges are not equal (if they are not 0).

        """
        n_atoms = self._topology.n_atoms
        n_pos = len(self._pos)
        n_vel = len(self._vel)
        n_forces = len(self._forces)
        n_charges = len(self._charges)

        n_atoms_list = np.array([n_atoms, n_pos, n_vel, n_forces, n_charges])
        n_atoms_list = n_atoms_list[n_atoms_list != 0]

        if n_atoms_list.size == 0:
            return 0

        if not np.all(n_atoms_list == n_atoms_list[0]):
            self.logger.error(
                (
                    "The number of atoms (or atoms in the topology), "
                    "positions, velocities, forces and charges must be equal."
                ),
                exception=AtomicSystemError
            )

        return int(n_atoms_list[0])

    @property
    @check_atoms_has_mass
    def atomic_masses(self) -> Np1DNumberArray:
        """Np1DNumberArray: The atomic masses of the atoms in the system."""
        return np.array([atom.mass for atom in self.atoms])

    @property
    def mass(self) -> Real:
        """Real: The total mass of the system."""
        return np.sum(self.atomic_masses)

    @property
    @check_atoms_pos
    @check_atoms_has_mass
    def center_of_mass(self) -> Np1DNumberArray:
        """Np1DNumberArray: The center of mass of the system."""
        if self.n_atoms == 0:
            return np.zeros(3)

        relative_pos = self.cell.image(self.pos - self.pos[0]) + self.pos[0]

        weighted_average_pos = np.sum(
            relative_pos * self.atomic_masses[:, None], axis=0
        ) / self.mass

        center_of_mass = self.cell.image(weighted_average_pos)

        return center_of_mass

    @property
    def combined_name(self) -> str:
        """str: The combined name of the atoms in the system."""
        return ''.join(atom.name for atom in self.atoms)

    @property
    def unique_element_names(self) -> List[str]:
        """List[str]: The unique element names of the atoms in the system."""
        return list({atom.element_name for atom in self.atoms})

    @property
    def build_custom_element(self) -> CustomElement:
        """CustomElement: The custom element of the atoms in the system."""
        return CustomElement(self.combined_name, -1, self.mass)
