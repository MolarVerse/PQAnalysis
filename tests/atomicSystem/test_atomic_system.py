"""
Tests for the AtomicSystem class.
"""

import pytest
import numpy as np

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, CustomElement, Element, Residue
from PQAnalysis.topology import Topology
from PQAnalysis.exceptions import PQTypeError, PQNotImplementedError
from PQAnalysis.type_checking import get_type_error_message

from PQAnalysis.atomic_system.exceptions import (
    AtomicSystemPositionsError,
    AtomicSystemMassError,
    AtomicSystemError,
)
from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray

from . import pytestmark  # pylint: disable=unused-import
from ..conftest import assert_logging_with_exception, assert_type_error_in_debug_mode



class TestAtomicSystem:

    """
    Tests for the AtomicSystem class.
    """

    def test__init__(self, caplog):
        """
        Test the __init__ method of the AtomicSystem class.

        Parameters
        ----------
        caplog : _pytest.logging.LogCaptureFixture
            The logging fixture.
        """
        system = AtomicSystem()
        assert system.atoms == []
        assert system.pos.shape == (0, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0, )
        assert system.cell == Cell()
        assert system.pbc is False
        assert system.n_atoms == 0
        assert system.atomic_masses.shape == (0, )
        assert np.isclose(system.mass, 0.0)
        assert np.allclose(system.center_of_mass, np.zeros(3))
        assert system.combined_name == ''

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]))
        assert system.atoms == []
        assert system.pos.shape == (2, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0, )
        assert system.cell == Cell()
        assert system.pbc is False
        assert system.n_atoms == 2
        assert np.allclose(system.pos, [[0, 0, 0], [1, 1, 1]])

        assert system.atomic_masses.shape == (0, )
        assert system.mass == 0

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            AtomicSystemPositionsError.message,
            AtomicSystemPositionsError,
            AtomicSystem.center_of_mass.__get__,
            system
        )

        assert system.combined_name == ''

        system = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        assert system.atoms == [Atom('C'), Atom('H')]
        assert system.pos.shape == (2, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0, )
        assert system.cell == Cell(0.75, 0.75, 0.75)
        assert system.pbc is True
        assert system.n_atoms == 2
        assert np.allclose(system.pos, [[0, 0, 0], [1, 1, 1]])
        assert np.allclose(system.atomic_masses, [12.0107, 1.00794])
        assert np.isclose(system.mass, 13.01864)
        assert np.allclose(
            system.center_of_mass, [0.01935571, 0.01935571, 0.01935571]
        )
        assert system.combined_name == 'CH'

        system = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]), atoms=[Atom('C'), Atom('H')]
        )

        assert np.allclose(
            system.center_of_mass,
            np.array([0.07742283, 0.07742283, 0.07742283])
        )

        system = AtomicSystem(
            atoms=[Atom('C', use_guess_element=False)],
            pos=np.array([[0, 0, 0]])
        )

        assert system.atoms == [Atom('C', use_guess_element=False)]

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            AtomicSystemMassError.message,
            AtomicSystemMassError,
            AtomicSystem.atomic_masses.__get__,
            system
        )

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            AtomicSystemMassError.message,
            AtomicSystemMassError,
            AtomicSystem.mass.__get__,
            system
        )

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            AtomicSystemMassError.message,
            AtomicSystemMassError,
            AtomicSystem.center_of_mass.__get__,
            system
        )

        atoms = [Atom('C1', use_guess_element=False), Atom('H')]
        topology = Topology(atoms=atoms)

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            (
                "Cannot initialize AtomicSystem with both atoms "
                "and topology arguments - they are mutually exclusive."
            ),
            AtomicSystemError,
            AtomicSystem,
            atoms=atoms,
            topology=topology
        )

        system = AtomicSystem(atoms=atoms)
        topology = Topology(atoms=atoms, residue_ids=np.array([0, 1]))
        assert system.topology != topology

        system.topology = topology
        assert system.topology == topology

        topology = Topology(atoms=[Atom('C1', use_guess_element=False)])
        with pytest.raises(ValueError) as exception:
            system.topology = topology
        assert str(
            exception.value
        ) == "The number of atoms already found in the AtomicSystem object have to be equal to the number of atoms in the new topology"

    def test__eq__(self):
        system1 = AtomicSystem()
        system2 = AtomicSystem()
        system3 = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]), atoms=[Atom('C'), Atom('H')]
        )
        system4 = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )
        system5 = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]), atoms=[Atom('C'), Atom('D')]
        )

        assert system1 == system2
        assert system1 != system3
        assert system1 != system4
        assert system3 != system4

        assert system1 != 1

        system1 = AtomicSystem(pos=np.array([[1, 1, 1]]))
        system2 = AtomicSystem(pos=np.array([[0, 0, 0]]))

        assert system1 != system2

        system1 = AtomicSystem(vel=np.array([[1, 1, 1]]))
        system2 = AtomicSystem(vel=np.array([[0, 0, 0]]))

        assert system1 != system2

        system1 = AtomicSystem(forces=np.array([[1, 1, 1]]))
        system2 = AtomicSystem(forces=np.array([[0, 0, 0]]))

        assert system1 != system2

        system1 = AtomicSystem(charges=np.array([1]))
        system2 = AtomicSystem(charges=np.array([0]))

        assert system1 != system2
        assert system3 != system5

    def test_isclose(self):
        system3 = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )
        system4 = AtomicSystem(
            pos=np.array([[0.00001, 0, 0], [1, 1.00001, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.750001, 0.75)
        )
        system5 = AtomicSystem(
            pos=np.array([[0, 0.01, 0], [1, 1.01, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.76, 0.75)
        )

        assert system3.isclose(system4, atol=1.0001e-5)
        assert system3.isclose(system4, rtol=1)
        assert not system3.isclose(system4, atol=1e-6)
        assert not system3.isclose(system4, rtol=1e-6)

        assert system3.isclose(system5, atol=1e-1)
        assert not system3.isclose(system5, rtol=1e-1)
        assert not system3.isclose(system5, atol=1e-3)
        assert not system3.isclose(system5, rtol=1e-2)

        system1 = AtomicSystem(
            pos=np.array([[0.0, 0, 0], [1, 100.1, 1]]),
            atoms=[Atom('C'), Atom('H')]
        )

        system2 = AtomicSystem(
            pos=np.array([[0.0, 0, 0], [1, 100.11, 1]]),
            atoms=[Atom('C'), Atom('H')]
        )

        assert system1.isclose(system2, atol=1.1e-2)
        assert system1.isclose(system2, rtol=1.1e-4)
        assert not system1.isclose(system2, atol=1.0e-3)
        assert not system1.isclose(system2, rtol=1.0e-5)

    def test__getitem__(self):
        system = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        assert system[0] == AtomicSystem(
            pos=np.array([[0, 0, 0]]),
            atoms=[Atom('C')],
            cell=Cell(0.75, 0.75, 0.75)
        )
        assert system[1] == AtomicSystem(
            pos=np.array([[1, 1, 1]]),
            atoms=[Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        with pytest.raises(IndexError):
            system[2]

        assert system[-1] == AtomicSystem(
            pos=np.array([[1, 1, 1]]),
            atoms=[Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        system = AtomicSystem(
            vel=np.array([[0, 0, 0], [1, 1, 1]]),
            forces=np.array([[0, 0, 0], [1, 1, 1]]),
            charges=np.array([0, 0]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        assert system[0] == AtomicSystem(
            vel=np.array([[0, 0, 0]]),
            forces=np.array([[0, 0, 0]]),
            charges=np.array([0]),
            atoms=[Atom('C')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        assert system[Atom(6)] == AtomicSystem(
            vel=np.array([[0, 0, 0]]),
            forces=np.array([[0, 0, 0]]),
            charges=np.array([0]),
            atoms=[Atom('C')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        system = AtomicSystem(
            vel=np.array([[0, 0, 0], [1, 1, 1]]),
            forces=np.array([[0, 0, 0], [1, 1, 1]]),
            charges=np.array([0, 0]),
            atoms=None,
            cell=Cell(0.75, 0.75, 0.75)
        )

        assert system[:] == system

    def test_n_atoms(self, caplog):
        system = AtomicSystem()
        assert system.n_atoms == 0

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]))
        assert system.n_atoms == 2

        system = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )
        assert system.n_atoms == 2

        system = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C')],
            cell=Cell(0.75, 0.75, 0.75)
        )

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            (
                "The number of atoms (or atoms in the topology), "
                "positions, velocities, forces and charges must be equal."
            ),
            AtomicSystemError,
            lambda system: system.n_atoms,
            system
        )

    def test__len__(self):
        system = AtomicSystem()
        assert len(system) == 0

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]))
        assert len(system) == 2

        system = AtomicSystem(
            pos=np.array([[0, 0, 0], [1, 1, 1]]),
            atoms=[Atom('C'), Atom('H')],
            cell=Cell(0.75, 0.75, 0.75)
        )
        assert len(system) == 2

    def test__str__(self):
        """
        Test the __str__ method of the AtomicSystem class.
        """
        pos = np.array([[0, 0, 0], [1, 1, 1]])
        vel = np.array([[0, 0, 0], [2, 2, 2]])
        charges = np.array([0, 1])
        forces = np.array([[0, 0, 0], [3, 3, 3]])
        cell = Cell(0.75, 0.75, 0.75)
        atoms = [Atom('C'), Atom('H')]

        system = AtomicSystem()
        assert str(
            system
        ) == f"AtomicSystem(topology=({system.topology}), cell=({system.cell}))"
        assert str(system) == repr(system)

        system = AtomicSystem(
            atoms=atoms,
            pos=pos,
            cell=cell,
            vel=vel,
            charges=charges,
            forces=forces
        )
        assert str(
            system
        ) == f"AtomicSystem(topology=({system.topology}), cell=({system.cell}))"
        assert str(system) == repr(system)

    def test_pos_setter(self, caplog):
        """
        Test the pos setter of the AtomicSystem class.

        Parameters
        ----------
        caplog : _pytest.logging.LogCaptureFixture
            The logging fixture.
        """

        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
        system.pos = np.array([[0, 0, 0], [1, 1, 1]])
        assert np.allclose(system.pos, np.array([[0, 0, 0], [1, 1, 1]]))

        def f(system, pos):
            system.pos = pos

        assert_logging_with_exception(
            caplog,
            "TypeChecking",
            "ERROR",
            get_type_error_message(
                "pos",
                np.array([0, 0, 0]),
                Np2DNumberArray,
            ),
            PQTypeError,
            f,
            system,
            np.array([0, 0, 0])
        )

        # should work without raising an exception
        assert_type_error_in_debug_mode(
            system.set_pos_no_checks, np.array([0, 0, 0])
        )

        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            (
                "The number of atoms in the AtomicSystem object have "
                "to be equal to the number of atoms in the new array "
                "in order to set the property."
            ),
            AtomicSystemError,
            f,
            system,
            np.array([[0, 0, 0]])
        )

        system = AtomicSystem(atoms=[Atom('C')], pos=np.array([[1, 1, 1]]))

        system.pos = np.array([[0, 0, 0]])
        assert np.allclose(system.pos, np.array([[0, 0, 0]]))

        system.set_pos_no_checks(np.array([[1, 1, 1]]))
        assert np.allclose(system.pos, np.array([[1, 1, 1]]))

    def test_vel_setter(self, caplog):
        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
        system.vel = np.array([[0, 0, 0], [1, 1, 1]])
        assert np.allclose(system.vel, np.array([[0, 0, 0], [1, 1, 1]]))

        def f(system, vel):
            system.vel = vel

        assert_logging_with_exception(
            caplog,
            "TypeChecking",
            "ERROR",
            get_type_error_message(
                "vel",
                np.array([0, 0, 0]),
                Np2DNumberArray,
            ),
            PQTypeError,
            f,
            system,
            np.array([0, 0, 0])
        )

        # should work without raising an exception
        assert_type_error_in_debug_mode(
            system.set_vel_no_checks, np.array([0, 0, 0])
        )

        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            (
                "The number of atoms in the AtomicSystem object have "
                "to be equal to the number of atoms in the new array "
                "in order to set the property."
            ),
            AtomicSystemError,
            f,
            system,
            np.array([[0, 0, 0]])
        )

        system = AtomicSystem(atoms=[Atom('C')], vel=np.array([[1, 1, 1]]))
        system.vel = np.array([[0, 0, 0]])
        assert np.allclose(system.vel, np.array([[0, 0, 0]]))

        system.set_vel_no_checks(np.array([[1, 1, 1]]))
        assert np.allclose(system.vel, np.array([[1, 1, 1]]))

    def test_force_setter(self, caplog):
        """
        Test the force setter of the AtomicSystem class.

        Parameters
        ----------
        caplog : _pytest.logging.LogCaptureFixture
            The logging fixture.
        """
        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
        system.forces = np.array([[0, 0, 0], [1, 1, 1]])
        assert np.allclose(system.forces, np.array([[0, 0, 0], [1, 1, 1]]))

        def f(system, forces):
            system.forces = forces

        assert_logging_with_exception(
            caplog,
            "TypeChecking",
            "ERROR",
            get_type_error_message(
                "forces",
                np.array([0, 0, 0]),
                Np2DNumberArray,
            ),
            PQTypeError,
            f,
            system,
            np.array([0, 0, 0])
        )

        # should work without raising an exception
        assert_type_error_in_debug_mode(
            system.set_forces_no_checks, np.array([0, 0, 0])
        )

        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            (
                "The number of atoms in the AtomicSystem object have "
                "to be equal to the number of atoms in the new array "
                "in order to set the property."
            ),
            AtomicSystemError,
            f,
            system,
            np.array([[0, 0, 0]])
        )

        system = AtomicSystem(atoms=[Atom('C')], forces=np.array([[1, 1, 1]]))
        system.forces = np.array([[0, 0, 0]])
        assert np.allclose(system.forces, np.array([[0, 0, 0]]))

        system.set_forces_no_checks(np.array([[1, 1, 1]]))
        assert np.allclose(system.forces, np.array([[1, 1, 1]]))

    def test_charge_setter(self, caplog):
        """
        Test the charge setter of the AtomicSystem class.

        Parameters
        ----------
        caplog : _pytest.logging.LogCaptureFixture
            The logging fixture.
        """
        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
        system.charges = np.array([0, 1])
        assert np.allclose(system.charges, np.array([0, 1]))

        def f(system, charges):
            system.charges = charges

        assert_logging_with_exception(
            caplog,
            "TypeChecking",
            "ERROR",
            get_type_error_message(
                "charges",
                np.array([[0, 0, 0]]),
                Np1DNumberArray,
            ),
            PQTypeError,
            f,
            system,
            np.array([[0, 0, 0]])
        )

        # should work without raising an exception
        assert_type_error_in_debug_mode(
            system.set_charges_no_checks, np.array([[0, 0, 0]])
        )

        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])

        assert_logging_with_exception(
            caplog,
            AtomicSystem.__qualname__,
            "ERROR",
            (
                "The number of atoms in the AtomicSystem object have "
                "to be equal to the number of atoms in the new array "
                "in order to set the property."
            ),
            AtomicSystemError,
            f,
            system,
            np.array([0, 0, 0])
        )

        system = AtomicSystem(atoms=[Atom('C')], charges=np.array([1]))
        system.charges = np.array([0])
        assert np.allclose(system.charges, np.array([0]))

        system.set_charges_no_checks(np.array([1]))
        assert np.allclose(system.charges, np.array([1]))

    def test_has_properties(self):
        """
        Test the has_properties method of the AtomicSystem class.
        """
        system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
        assert not system.has_pos
        assert not system.has_vel
        assert not system.has_forces
        assert not system.has_charges

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')], pos=np.array([[0, 0, 0], [1, 1, 1]])
        )
        assert system.has_pos
        assert not system.has_vel
        assert not system.has_forces
        assert not system.has_charges

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')], vel=np.array([[0, 0, 0], [1, 1, 1]])
        )
        assert not system.has_pos
        assert system.has_vel
        assert not system.has_forces
        assert not system.has_charges

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')],
            forces=np.array([[0, 0, 0], [1, 1, 1]])
        )
        assert not system.has_pos
        assert not system.has_vel
        assert system.has_forces
        assert not system.has_charges

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')], charges=np.array([0, 1])
        )
        assert not system.has_pos
        assert not system.has_vel
        assert not system.has_forces
        assert system.has_charges

    def test_energy(self):
        """
        Test the energy property of the AtomicSystem class.
        """
        system = AtomicSystem(energy=1.0)

        assert np.isclose(system.energy, 1.0)
        assert system.has_energy

        system.energy = 2.0
        assert np.isclose(system.energy, 2.0)
        assert system.has_energy

        system = AtomicSystem()
        assert not system.has_energy
        assert system.energy is None

    def test_virial(self):
        """
        Test the virial property of the AtomicSystem class.
        """
        system = AtomicSystem(
            virial=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        )

        assert np.allclose(
            system.virial, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        )
        assert system.has_virial

        system.virial = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        assert np.allclose(
            system.virial, np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        )
        assert system.has_virial

        system = AtomicSystem()
        assert not system.has_virial
        assert system.virial is None

    def test_stress(self):
        """
        Test the stress property of the AtomicSystem class.
        """
        system = AtomicSystem(
            stress=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        )

        assert np.allclose(
            system.stress, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        )
        assert system.has_stress

        system.stress = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        assert np.allclose(
            system.stress, np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        )
        assert system.has_stress

        system = AtomicSystem()
        assert not system.has_stress
        assert system.stress is None

    def test_center_of_mass_resiudes(self):
        system = AtomicSystem()

        with pytest.raises(AtomicSystemError) as exception:
            system.center_of_mass_residues  # pylint: disable=pointless-statement
        assert str(exception.value) == "No residues in the system."

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')],
            vel=np.array([[0, 0, 0], [1, 1, 1]]),
        )

        with pytest.raises(PQNotImplementedError) as exception:
            system.center_of_mass_residues  # pylint: disable=pointless-statement
        assert str(exception.value) == (
            "Center of mass of residues not implemented for "
            "systems with forces, velocities or charges."
        )

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')],
            forces=np.array([[0, 0, 0], [1, 1, 1]]),
        )

        with pytest.raises(PQNotImplementedError) as exception:
            system.center_of_mass_residues  # pylint: disable=pointless-statement
        assert str(exception.value) == (
            "Center of mass of residues not implemented for "
            "systems with forces, velocities or charges."
        )

        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H')],
            charges=np.array([0, 1]),
        )

        with pytest.raises(PQNotImplementedError) as exception:
            system.center_of_mass_residues  # pylint: disable=pointless-statement
        assert str(exception.value) == (
            "Center of mass of residues not implemented for "
            "systems with forces, velocities or charges."
        )

        topology = Topology(
            atoms=[Atom('C'), Atom('H')],
            residue_ids=np.array([0, 1]),
        )

        system = AtomicSystem(topology=topology)

        assert system.center_of_mass_residues == system

        topology = Topology(
            atoms=[
                Atom('C', use_guess_element=False),
                Atom('H', use_guess_element=False)
            ],
            residue_ids=np.array([0, 1]),
        )

        system = AtomicSystem(topology=topology)

        assert system.center_of_mass_residues == system

        reference_residues = [
            Residue(
                name="CH2",
                elements=[Element('C'), Element('H'), Element('H')],
                residue_id=1,
                total_charge=0,
                atom_types=np.array([0, 1, 1]),
                partial_charges=np.array([0, 0, 0])
            ),
            Residue(
                name="H",
                elements=[Element('H')],
                residue_id=2,
                total_charge=0,
                atom_types=np.array([1]),
                partial_charges=np.array([0])
            )
        ]

        topology = Topology(
            atoms=[Atom('C'), Atom('H'), Atom('H'), Atom('H')],
            residue_ids=np.array([1, 1, 1, 2]),
            reference_residues=reference_residues
        )

        system = AtomicSystem(
            topology=topology,
            pos=np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
        )

        com_system = system.center_of_mass_residues

        assert com_system.atoms == [
            Atom(CustomElement('CH2', -1, 14.027)), Atom(Element('H'))
        ]
        assert np.allclose(
            com_system.pos,
            np.array([[0.21557785, 0.21557785, 0.21557785], [3, 3, 3]]),
        )
