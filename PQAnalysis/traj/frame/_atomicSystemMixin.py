"""
A module containing a Mixin Class with all methods and properties of an atomic system that are forwarded to the frame class.
"""

from beartype.typing import List

from PQAnalysis.core import Cell, Atom
from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.topology import Topology
from PQAnalysis.atomicSystem import AtomicSystem


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

    def fit_atomic_system(self,
                          system: AtomicSystem,
                          number_of_additions: PositiveInt = 1,
                          max_iterations: PositiveInt = 100,
                          distance_cutoff: PositiveReal = 1.0,
                          max_displacement: PositiveReal | Np1DNumberArray = 0.1,
                          rotation_angle_step: PositiveInt = 10,
                          ) -> AtomicSystem | List[AtomicSystem]:
        """
        Fit the positions of the system to the positions of another system.

        First a random center of mass is chosen and a random displacement is applied to the system. Then the system is rotated in all possible ways and the distances between the atoms are checked. If the distances are larger than the distance cutoff, the system is fitted.

        Parameters
        ----------
        system : AtomicSystem
            The system that should be fitted into the positions of the AtomicSystem.
        number_of_additions : PositiveInt, optional
            The number of times the system should be added to the AtomicSystem, by default 1
        max_iterations : PositiveInt, optional
            The maximum number of iterations to try to fit the system into the positions of the AtomicSystem, by default 100
        distance_cutoff : PositiveReal, optional
            The distance cutoff for the fitting, by default 1.0
        max_displacement : PositiveReal | Np1DNumberArray, optional
            The maximum displacement percentage for the fitting, by default 0.1
        rotation_angle_step : PositiveInt, optional
            The angle step for the rotation of the system, by default 10

        Returns
        -------
        AtomicSystem
            The fitted AtomicSystem.

        Raises
        ------
        AtomicSystemError
            If the AtomicSystem has a vacuum cell.
        ValueError
            If the maximum displacement percentage is negative.
        AtomicSystemError
            If the system could not be fitted into the positions of the AtomicSystem within the maximum number of iterations.
        """
        return self.system.fit_atomic_system(
            system,
            number_of_additions=number_of_additions,
            max_iterations=max_iterations,
            distance_cutoff=distance_cutoff,
            max_displacement=max_displacement,
            rotation_angle_step=rotation_angle_step
        )
