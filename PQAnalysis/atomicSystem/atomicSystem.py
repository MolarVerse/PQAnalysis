"""
A module containing the AtomicSystem class
"""

from __future__ import annotations

import numpy as np
import itertools
import secrets

from scipy.spatial.transform import Rotation
from beartype.typing import Any, List

from ._properties import _PropertiesMixin
from ._standardProperties import _StandardPropertiesMixin
from ._positions import _PositionsMixin
from .exceptions import AtomicSystemError

from PQAnalysis.core import Atom, Atoms, Cell, distance
from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray, Np1DIntArray
from PQAnalysis.topology import Topology
from PQAnalysis.types import PositiveReal, PositiveInt


class AtomicSystem(_PropertiesMixin, _StandardPropertiesMixin, _PositionsMixin):
    """
    A class for storing atomic systems.

    It contains the standard properties of an atomic system (i.e. positions, velocities, forces, charges, topology and cell). The AtomicSystem class can be used as a container for the standard properties of any molecular/atomic system.

    Notes
    -----
    An atomic system does not have to containing any positions, velocities, forces and so forth. The only requirement is that the number of atoms in the topology is equal to the number of entries in the positions, velocities, forces and charges arrays. If e.g. only a system containing information of the velocities is needed, the positions, forces and charges arrays can be left empty (i.e. np.zeros((0, 3)) and np.zeros(0)). The same goes for the other properties. An empty cell can be created with Cell() and represents a system without periodic boundary conditions. (For more information see the documentation of :py:class:`~PQAnalysis.core.cell.cell.Cell`). As the topology is can be really and complex and most of the cases really specific to the system, here no further information is given. (For more information see the documentation of :py:class:`~PQAnalysis.topology.topology.Topology`). Furthermore for this reason if no specialization of the topology is needed, the atomic system can be initialized with only a list of atoms (see examples, and the documentation of :py:class:`~PQAnalysis.core.atom.atom.Atom`).


    Inherits from the Mixins: _PropertiesMixin, _StandardPropertiesMixin, _IndexingMixin, _PositionsMixin
        - The _StandardPropertiesMixin contains the standard properties of an atomic system (i.e. standard getter and setter methods).
        - The _PropertiesMixin contains special properties derived from the standard properties
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

        if topology is not None and atoms is not None:
            raise ValueError(
                "Cannot initialize AtomicSystem with both atoms and topology arguments - they are mutually exclusive."
            )

        if atoms is None and topology is None:
            topology = Topology()
        elif topology is None:
            topology = Topology(atoms=atoms)

        self._topology = topology
        self._pos = np.zeros((0, 3)) if pos is None else pos
        self._vel = np.zeros((0, 3)) if vel is None else vel
        self._forces = np.zeros((0, 3)) if forces is None else forces
        self._charges = np.zeros(0) if charges is None else charges
        self._cell = cell

    def fit_atomic_system(self,
                          system: AtomicSystem,
                          number_of_additions: PositiveInt = 1,
                          max_iterations: PositiveInt = 100,
                          distance_cutoff: PositiveReal = 1.0,
                          max_displacement: PositiveReal | Np1DNumberArray = 0.1,
                          rotation_angle_step: PositiveInt = 10,
                          ) -> List[AtomicSystem] | AtomicSystem:
        """
        Fit the positions of the system to the positions of another system.

        Parameters
        ----------
        system : AtomicSystem
            The system that should be fitted into the positions of the AtomicSystem.
        number_of_additions : PositiveInt, optional
            The number of times the system should be fitted into the positions of the AtomicSystem, by default 1
        max_iterations : PositiveInt, optional
            The maximum number of iterations to try to fit the system into the positions of the AtomicSystem, by default 100
        distance_cutoff : PositiveReal, optional
            The distance cutoff for the fitting, by default 1.0
        max_displacement : PositiveReal | Np1DNumberArray, optional
            The maximum displacement percentage for the fitting, by default 0.1
        rotation_angle_step : PositiveInt, optional
            The angle step for the rotation of the system, by default 10

        First a random center of mass is chosen and a random displacement is applied to the system. Then the system is rotated in all possible ways and the distances between the atoms are checked. If the distances are larger than the distance cutoff, the system is fitted.

        Returns
        -------
        List[AtomicSystem] | AtomicSystem
            The fitted AtomicSystem(s). If number_of_additions is 1, a single AtomicSystem is returned, otherwise a list of AtomicSystems is returned.

        Raises
        ------
        AtomicSystemError
            If the AtomicSystem has a vacuum cell.
        ValueError
            If the maximum displacement percentage is negative.
        AtomicSystemError
            If the system could not be fitted into the positions of the AtomicSystem within the maximum number of iterations.
        """

        systems = []

        for _ in range(number_of_additions):
            systems.append(
                self._fit_atomic_system(
                    system=system,
                    max_iterations=max_iterations,
                    distance_cutoff=distance_cutoff,
                    max_displacement=max_displacement,
                    rotation_angle_step=rotation_angle_step
                )
            )

        return systems if number_of_additions > 1 else systems[0]

    def _fit_atomic_system(self,
                           system: AtomicSystem,
                           max_iterations: PositiveInt = 100,
                           distance_cutoff: PositiveReal = 1.0,
                           max_displacement: PositiveReal | Np1DNumberArray = 0.1,
                           rotation_angle_step: PositiveInt = 10,
                           ) -> AtomicSystem:
        """
        Fit the positions of the system to the positions of another system.

        First a random center of mass is chosen and a random displacement is applied to the system. Then the system is rotated in all possible ways and the distances between the atoms are checked. If the distances are larger than the distance cutoff, the system is fitted.

        Parameters
        ----------
        system : AtomicSystem
            The system that should be fitted into the positions of the AtomicSystem.
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

        if self.cell.is_vacuum:
            raise AtomicSystemError(
                "Cannot fit into positions of a system with a vacuum cell."
            )

        if isinstance(max_displacement, float):
            max_displacement = np.array([max_displacement] * 3)

        if np.any(max_displacement < 0.0):
            raise ValueError(
                "The maximum displacement percentage must be a positive number."
            )

        iter_converged = None
        seed = secrets.randbits(128)
        rng = np.random.default_rng(seed=seed)

        for _iter in range(max_iterations):
            com = rng.random(3)
            com = com * self.cell.box_lengths - self.cell.box_lengths / 2

            rel_com_positions = system.pos - system.center_of_mass

            displacement = rng.random(3)
            displacement *= 2 * max_displacement
            displacement -= max_displacement

            new_pos = rel_com_positions + com + displacement

            rotation = Rotation.random()

            for x, y, z in itertools.product(range(0, 360, rotation_angle_step), repeat=3):
                rotation_angles = rotation.as_euler(
                    'xyz',
                    degrees=True
                )
                rotation_angles += np.array([x, y, z])
                rotation = Rotation.from_euler(
                    'xyz',
                    rotation_angles,
                    degrees=True
                )
                new_pos = rotation.apply(new_pos)

                distances = distance(self.pos, new_pos, self.cell)

                if np.all(distances > distance_cutoff):
                    iter_converged = _iter
                    break

            if iter_converged is not None:
                break

        if iter_converged is None:
            raise AtomicSystemError(
                "Could not fit the positions of the system. Try increasing the maximum number of iterations."
            )
        else:
            print(f"Fit converged after {_iter} iterations.")
            system = system.copy()
            system.pos = new_pos
            return system

    def copy(self) -> AtomicSystem:
        """
        Returns a copy of the AtomicSystem.

        Returns
        -------
        AtomicSystem
            A copy of the AtomicSystem.
        """
        return AtomicSystem(
            pos=self.pos,
            vel=self.vel,
            forces=self.forces,
            charges=self.charges,
            cell=self.cell,
            topology=self.topology
        )

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

        if self.n_atoms != other.n_atoms:
            return False

        if self.topology != other.topology:
            return False

        if self.cell != other.cell:
            return False

        if not np.allclose(self.pos, other.pos):
            return False

        if not np.allclose(self.vel, other.vel):
            return False

        if not np.allclose(self.forces, other.forces):
            return False

        if not np.allclose(self.charges, other.charges):
            return False

        return True

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
            key = np.array([key])

        keys = np.array(range(self.n_atoms))[key] if isinstance(
            key, slice) else np.array(key)

        pos = self.pos[keys] if np.shape(self.pos)[0] > 0 else None
        vel = self.vel[keys] if np.shape(self.vel)[0] > 0 else None
        forces = self.forces[keys] if np.shape(self.forces)[0] > 0 else None
        charges = self.charges[keys] if np.shape(self.charges)[0] > 0 else None

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
