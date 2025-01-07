"""
A module containing the AtomicSystem class
"""

import itertools
import logging
import sys
import numpy as np

from scipy.spatial.transform import Rotation
from beartype.typing import Any

# just for forwardref type hinting
from beartype.typing import List  # pylint: disable=unused-import

from PQAnalysis.core import Atom, Atoms, Cell, distance
from PQAnalysis.topology import Topology
from PQAnalysis.types import PositiveReal, PositiveInt, Bool
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.exceptions import PQNotImplementedError
from PQAnalysis.utils.random import get_random_seed
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.utils.math import allclose_vectorized
from PQAnalysis import __package_name__

from PQAnalysis.types import (
    Np2DNumberArray,
    Np1DNumberArray,
    Np1DIntArray,
    Np3x3NumberArray,
)

from ._properties import _PropertiesMixin
from ._standard_properties import _StandardPropertiesMixin
from ._positions import _PositionsMixin
from .exceptions import AtomicSystemError



class AtomicSystem(
    _PropertiesMixin,
    _StandardPropertiesMixin,
    _PositionsMixin,
):

    """
    A class for storing atomic systems.

    It contains the standard properties of an atomic system
    (i.e. positions, velocities, forces, charges, topology and cell).
    The AtomicSystem class can be used as a container for the 
    standard properties of any molecular/atomic system.

    Notes
    -----
    An atomic system does not have to contain any positions,
    velocities, forces and so on. The only requirement is that
    the number of atoms in the topology is equal to the number of
    entries in the positions, velocities, forces and charges 
    arrays. For example, if a system containing only information about
    the velocities is needed, the positions, forces and charges arrays
    can be left empty (i.e. np.zeros((0, 3)) and np.zeros(0)).
    The same goes for the other properties. An empty cell can be
    created with Cell() and represents a system without periodic
    boundary conditions. (For more information see the 
    documentation of :py:class:`~PQAnalysis.core.cell.cell.Cell`).
    Since, the topology can be really complex and most of the 
    cases really specific to the system, no further 
    information is given here. (For more information see the 
    documentation of :py:class:`~PQAnalysis.topology.topology.Topology`).
    For this reason, if no specialization of the topology
    is needed, the atomic system can be initialized with only a list
    of atoms (see examples, and the documentation of 
    :py:class:`~PQAnalysis.core.atom.atom.Atom`).


    Inherits from the Mixins: _PropertiesMixin, _StandardPropertiesMixin,
    _IndexingMixin, _PositionsMixin
        - The _StandardPropertiesMixin contains the standard properties of an atomic
        system (i.e. standard getter and setter methods).
        - The _PropertiesMixin contains special properties derived from the standard properties
        - The _PositionsMixin contains methods for computing properties based
        on the positions of the atoms


    Examples
    --------
    >>> import numpy as np
    >>> from PQAnalysis.core import Atom
    >>> from PQAnalysis.atomic_system import AtomicSystem

    >>> atoms = [Atom('C1', use_guess_element=False), Atom('C2', use_guess_element=False)]
    >>> AtomicSystem(atoms=atoms, pos=np.array([[0, 0, 0], [1, 0, 0]]))
    AtomicSystem(topology=(Topology with 2 atoms...), cell=(Cell()))

    >>> AtomicSystem()
    AtomicSystem(topology=(Topology with 0 atoms...), cell=(Cell()))

    >>> AtomicSystem(topology=Topology(atoms=[Atom('C'), Atom('C')]))
    AtomicSystem(topology=(Topology with 2 atoms...), cell=(Cell()))
    """

    logging.basicConfig(level=logging.INFO)
    fitting_logger = logging.getLogger(__name__ + ".fit_atomic_system")
    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(message)s"))
    fitting_logger.addHandler(handler)
    fitting_logger.propagate = False

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        atoms: Atoms | None = None,
        pos: Np2DNumberArray | None = None,
        vel: Np2DNumberArray | None = None,
        forces: Np2DNumberArray | None = None,
        charges: Np1DNumberArray | None = None,
        topology: Topology | None = None,
        energy: float | None = None,
        virial: Np3x3NumberArray | None = None,
        stress: Np3x3NumberArray | None = None,
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
        energy : float, optional
            The energy of the system, by default None
        virial : Np3x3NumberArray, optional
            The virial of the system, by default None
        stress : Np3x3NumberArray, optional
            The stress of the system, by default None
        cell : Cell, optional
            The unit cell of the system. Defaults to a Cell with no periodic
            boundary conditions, by default Cell()

        Raises
        ------
        ValueError
            If both atoms and topology are given.
        """

        if topology is not None and atoms is not None:
            self.logger.error(
                (
                    "Cannot initialize AtomicSystem with both atoms and topology "
                    "arguments - they are mutually exclusive."
                ),
                exception=AtomicSystemError
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
        self._energy = energy
        self._virial = virial
        self._stress = stress
        self._cell = cell

    #TODO: check why dynamic formatting does
    #      not work here for "AtomicSystem"
    #
    # @runtime_type_checking
    def fit_atomic_system(
        self,
        system: "AtomicSystem",
        number_of_additions: PositiveInt = 1,
        max_iterations: PositiveInt = 100,
        distance_cutoff: PositiveReal = 1.0,
        max_displacement: PositiveReal | Np1DNumberArray = 0.1,
        rotation_angle_step: PositiveInt = 10,
    ) -> "List[AtomicSystem] | AtomicSystem":
        """
        Fit the positions of the system to the positions of another system.

        Parameters
        ----------
        system : AtomicSystem
            The system that should be fitted into the positions of the AtomicSystem.
        number_of_additions : PositiveInt, optional
            The number of times the system should be fitted into the
            positions of the AtomicSystem, by default 1
        max_iterations : PositiveInt, optional
            The maximum number of iterations to try to fit the system
            into the positions of the AtomicSystem, by default 100
        distance_cutoff : PositiveReal, optional
            The distance cutoff for the fitting, by default 1.0
        max_displacement : PositiveReal | Np1DNumberArray, optional
            The maximum displacement percentage for the fitting, by default 0.1
        rotation_angle_step : PositiveInt, optional
            The angle step for the rotation of the system,
            by default 10

        First a random center of mass is chosen and a random displacement
        is applied to the system. Then the system is rotated in all possible
        ways and the distances between the atoms are checked. If the
        distances are larger than the distance cutoff, the system is fitted.

        Returns
        -------
        List[AtomicSystem] | AtomicSystem
            The fitted AtomicSystem(s). If number_of_additions is 1,
            a single AtomicSystem is returned, otherwise a list of
            AtomicSystems is returned.

        Raises
        ------
        AtomicSystemError
            If the AtomicSystem has a vacuum cell.
        ValueError
            If the maximum displacement percentage is negative.
        AtomicSystemError
            If the system could not be fitted into the positions
            of the AtomicSystem within the maximum number of iterations.
        """

        systems = []

        for i in range(number_of_additions):
            # concatenate the positions of this system with the
            # positions of all systems that have been fitted so far
            positions_to_fit_into = np.concatenate(
                [system.pos for system in systems] + [self.pos]
            )

            self.fitting_logger.info(
                (
                    f"Performing fitting for {i + 1}/{number_of_additions} "
                    "addition(s)."
                )
            )

            systems.append(
                self._fit_atomic_system(
                    positions_to_fit_into=positions_to_fit_into,
                    system=system,
                    max_iterations=max_iterations,
                    distance_cutoff=distance_cutoff,
                    max_displacement=max_displacement,
                    rotation_angle_step=rotation_angle_step
                )
            )

        return systems if number_of_additions > 1 else systems[0]

    def _fit_atomic_system(
        self,
        positions_to_fit_into: Np2DNumberArray,
        system: "AtomicSystem",
        max_iterations: PositiveInt = 100,
        distance_cutoff: PositiveReal = 1.0,
        max_displacement: PositiveReal | Np1DNumberArray = 0.1,
        rotation_angle_step: PositiveInt = 10,
    ) -> "AtomicSystem":
        """
        Fit the positions of the system to the positions of another system.

        First a random center of mass is chosen and a random displacement
        is applied to the system. Then the system is rotated in all possible
        ways and the distances between the atoms are checked. If the
        distances are larger than the distance cutoff, the system is fitted.

        Parameters
        ----------
        positions_to_fit_into : Np2DNumberArray
            The positions of the systems were the new system should be fitted into.
        system : AtomicSystem
            The system that should be fitted into the positions of the AtomicSystem.
        max_iterations : PositiveInt, optional
            The maximum number of iterations to try to fit the system into the
            positions of the AtomicSystem, by default 100
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
            If the system could not be fitted into the positions
            of the AtomicSystem within the maximum number of iterations.
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
        seed = get_random_seed()
        rng = np.random.default_rng(seed=seed)

        for _iter in range(max_iterations):
            com = rng.random(3)
            com = com * self.cell.box_lengths - self.cell.box_lengths / 2

            rel_com_positions = system.pos - system.center_of_mass

            displacement = rng.random(3)
            displacement *= 2 * max_displacement
            displacement -= max_displacement

            new_pos = rel_com_positions + com + displacement

            rotation = Rotation.random(random_state=rng)

            for x, y, z in itertools.product(range(0, 360, rotation_angle_step), repeat=3):
                rotation_angles = rotation.as_euler('xyz', degrees=True)
                rotation_angles += np.array([x, y, z])
                rotation = Rotation.from_euler(
                    'xyz',
                    rotation_angles,
                    degrees=True,
                )
                new_pos = rotation.apply(new_pos)

                distances = distance(positions_to_fit_into, new_pos, self.cell)

                if np.all(distances > distance_cutoff):
                    iter_converged = _iter
                    break

            if iter_converged is not None:
                break

        if iter_converged is None:
            raise AtomicSystemError(
                "Could not fit the positions of the system. "
                "Try increasing the maximum number of iterations."
            )

        self.fitting_logger.info(
            f"\tConverged after {_iter + 1} iterations.\n"
        )
        system = system.copy()
        system.pos = new_pos
        system.cell = self.cell
        system.image()
        return system

    @property
    def center_of_mass_residues(self) -> "AtomicSystem":
        """
        Computes the center of mass of the residues in the system.

        Returns
        -------
        AtomicSystem
            The center of mass of the residues in the system.
            
        Raises:
        -------
        AtomicSystemError
            If the number of residues in the system is not a 
            multiple of the number of atoms.
        PQNotImplementedError
            if system has forces, velocities or charges.
            
        TODO:
        -----
        Include also center of mass velocities, forces and so on...
        """
        if self.has_pos:
            residue_pos = np.zeros(
                (len(self.topology.residue_atom_indices), 3)
            )
        else:
            residue_pos = np.zeros((0, 3))

        residue_atoms = []

        if self.has_forces or self.has_vel or self.has_charges:
            self.logger.error(
                (
                    "Center of mass of residues not implemented for "
                    "systems with forces, velocities or charges."
                ),
                exception=PQNotImplementedError,
            )

        if len(self.topology.residue_ids) == 0:
            self.logger.error(
                "No residues in the system.",
                exception=AtomicSystemError,
            )

        if len(self.topology.residue_ids) == 1:
            return self.copy()

        for i, residue_indices in enumerate(self.topology.residue_atom_indices):

            residue_system = self[residue_indices]

            # check if residue_system has more than one atom otherwise return atom element
            if (
                residue_system.n_atoms != 1 and
                len(self.topology.residues) != 0
            ):
                custom_element = residue_system.build_custom_element
                custom_element.symbol = self.topology.residues[i].name
                custom_atom = Atom(custom_element)
            else:
                custom_atom = residue_system.atoms[0]

            residue_atoms.append(custom_atom)

            if residue_system.has_pos:
                residue_pos[i] = residue_system.center_of_mass

        topology = Topology(
            atoms=residue_atoms,
            residue_ids=self.topology.residue_ids_per_residue
        )

        return AtomicSystem(
            pos=residue_pos,
            cell=self.cell,
            topology=topology,
        )

    # TODO: refactor or discard this method
    def compute_com_atomic_system(self, group=None) -> "AtomicSystem":
        """
        Computes a new AtomicSystem with the center of mass of the system or groups of atoms.

        Parameters
        ----------
        group : int, optional
            group of atoms to compute the center of mass of, by default None (all atoms)

        Returns
        -------
        AtomicSystem
            A new AtomicSystem with the center of mass of the system or groups of atoms.

        Raises
        ------
        AtomicSystemError
            If the number of atoms in the selection is not a multiple of group.
        """
        if group is None:
            group = self.n_atoms

        elif self.n_atoms % group != 0:
            raise AtomicSystemError(
                'Number of atoms in selection is not a multiple of group.'
            )

        pos = []
        names = []

        for i in range(0, self.n_atoms, group):
            atomic_system = AtomicSystem(
                atoms=self.atoms[i:i + group],
                pos=self.pos[i:i + group],
                cell=self.cell
            )

            pos.append(atomic_system.center_of_mass)
            names.append(atomic_system.combined_name)

        names = [Atom(name, use_guess_element=False) for name in names]

        return AtomicSystem(pos=np.array(pos), atoms=names, cell=self.cell)

    def copy(self) -> "AtomicSystem":
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

    def __eq__(self, other: Any) -> Bool:
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

        return self.isclose(other)

    @runtime_type_checking
    def isclose(
        self,
        other: Any,
        rtol: PositiveReal = 1e-9,
        atol: PositiveReal = 0.0,
    ) -> Bool:
        """
        Checks whether the AtomicSystem is close to another AtomicSystem.

        Parameters
        ----------
        other : AtomicSystem
            The other AtomicSystem to compare to.
        rtol : PositiveReal, optional
            The relative tolerance, by default 1e-9
        atol : PositiveReal, optional
            The absolute tolerance, by default 0.0

        Returns
        -------
        bool
            Whether the AtomicSystem is close to the other AtomicSystem.
        """

        if not isinstance(other, AtomicSystem):
            return False

        if self.n_atoms != other.n_atoms:
            return False

        if self.topology != other.topology:
            return False

        if not self.cell.isclose(other.cell, rtol=rtol, atol=atol):
            return False

        if not allclose_vectorized(self.pos, other.pos, rtol=rtol, atol=atol):
            return False

        if not allclose_vectorized(self.vel, other.vel, rtol=rtol, atol=atol):
            return False

        if not allclose_vectorized(
            self.forces, other.forces, rtol=rtol, atol=atol
        ):
            return False

        if not allclose_vectorized(
            self.charges, other.charges, rtol=rtol, atol=atol
        ):
            return False

        return True

    def __getitem__(
        self, key: Atom | int | slice | Np1DIntArray
    ) -> "AtomicSystem":
        """
        Returns a new AtomicSystem with the given key.

        Examples
        --------
        >>> import numpy as np
        >>> from PQAnalysis.core import Atom, Residue
        >>> from PQAnalysis.atomic_system import AtomicSystem

        >>> atoms = [Atom('C'), Atom('C')]
        >>> pos = np.array([[0, 0, 0], [1, 0, 0]])
        >>> system = AtomicSystem(atoms=atoms, pos=pos)
        >>> system[0]
        AtomicSystem(topology=(Topology with 1 atoms...), cell=(Cell()))

        >>> system[0] == AtomicSystem(atoms=[Atom('C')], pos=np.array([[0, 0, 0]]))
        True

        >>> system[0:2]
        AtomicSystem(topology=(Topology with 2 atoms...), cell=(Cell()))

        >>> atoms = [Atom('C'), Atom('C')]
        >>> pos = np.array([[0, 0, 0], [1, 0, 0]])
        >>> system[0:2] == AtomicSystem(atoms=atoms, pos=pos)
        True

        >>> system[np.array([0, 1])]
        AtomicSystem(topology=(Topology with 2 atoms ...), cell=(Cell()))

        >>> atoms = [Atom('C'), Atom('C')]
        >>> pos = np.array([[0, 0, 0], [1, 0, 0]])
        >>> system[np.array([0, 1])] == AtomicSystem(atoms=atoms, pos=pos)
        True

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

        keys = np.array(range(self.n_atoms)
                        )[key] if isinstance(key, slice) else np.array(key)

        pos = self.pos[keys] if np.shape(self.pos)[0] > 0 else None
        vel = self.vel[keys] if np.shape(self.vel)[0] > 0 else None
        forces = self.forces[keys] if np.shape(self.forces)[0] > 0 else None
        charges = self.charges[keys] if np.shape(self.charges)[0] > 0 else None

        return AtomicSystem(
            pos=pos,
            vel=vel,
            forces=forces,
            charges=charges,
            cell=self.cell,
            topology=self.topology[keys]
        )

    def __len__(self) -> int:
        """
        Returns the number of atoms in the AtomicSystem.

        Returns
        -------
        int
            The number of atoms in the AtomicSystem.
        """
        return self.n_atoms

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

    __iter__ = None  # To avoid iteration over the AtomicSystem object
