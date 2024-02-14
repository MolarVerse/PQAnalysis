import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.atomicSystem.exceptions import AtomicSystemPositionsError, AtomicSystemMassError
from PQAnalysis.core import Atom, Cell
from PQAnalysis.topology import Topology


class TestAtomicSystem:
    def test__init__(self):
        system = AtomicSystem()
        assert system.atoms == []
        assert system.pos.shape == (0, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0,)
        assert system.cell == Cell()
        assert system.PBC is False
        assert system.n_atoms == 0
        assert system.atomic_masses.shape == (0,)
        assert np.isclose(system.mass, 0.0)
        assert np.allclose(system.center_of_mass, np.zeros(3))
        assert system.combined_name == ''

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]))
        assert system.atoms == []
        assert system.pos.shape == (2, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0,)
        assert system.cell == Cell()
        assert system.PBC is False
        assert system.n_atoms == 2
        assert np.allclose(system.pos, [[0, 0, 0], [1, 1, 1]])

        assert system.atomic_masses.shape == (0,)
        assert system.mass == 0

        with pytest.raises(AtomicSystemPositionsError) as exception:
            system.center_of_mass
        assert str(
            exception.value) == AtomicSystemPositionsError.message

        assert system.combined_name == ''

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                              atoms=[Atom('C'), Atom('H')], cell=Cell(0.75, 0.75, 0.75))

        assert system.atoms == [Atom('C'), Atom('H')]
        assert system.pos.shape == (2, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0,)
        assert system.cell == Cell(0.75, 0.75, 0.75)
        assert system.PBC is True
        assert system.n_atoms == 2
        assert np.allclose(system.pos, [[0, 0, 0], [1, 1, 1]])
        assert np.allclose(system.atomic_masses, [12.0107, 1.00794])
        assert np.isclose(system.mass, 13.01864)
        assert np.allclose(system.center_of_mass, [
                           0.01935571, 0.01935571, 0.01935571])
        assert system.combined_name == 'CH'

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                              atoms=[Atom('C'), Atom('H')])

        assert np.allclose(system.center_of_mass, np.array(
            [0.07742283, 0.07742283, 0.07742283]))

        system = AtomicSystem(
            atoms=[Atom('C', use_guess_element=False)], pos=np.array([[0, 0, 0]]))

        assert system.atoms == [Atom('C', use_guess_element=False)]
        with pytest.raises(AtomicSystemMassError) as exception:
            system.atomic_masses
        assert str(
            exception.value) == AtomicSystemMassError.message
        with pytest.raises(AtomicSystemMassError) as exception:
            system.mass
        assert str(
            exception.value) == AtomicSystemMassError.message
        with pytest.raises(AtomicSystemMassError) as exception:
            system.center_of_mass
        assert str(
            exception.value) == AtomicSystemMassError.message

        atoms = [Atom('C1', use_guess_element=False), Atom('H')]
        topology = Topology(atoms=atoms)
        with pytest.raises(ValueError) as exception:
            AtomicSystem(atoms=atoms, topology=topology)
        assert str(
            exception.value) == "Cannot initialize AtomicSystem with both atoms and topology arguments - they are mutually exclusive."

        system = AtomicSystem(atoms=atoms)
        topology = Topology(atoms=atoms, residue_ids=np.array([0, 1]))
        assert system.topology != topology

        system.topology = topology
        assert system.topology == topology

        topology = Topology(atoms=[Atom('C1', use_guess_element=False)])
        with pytest.raises(ValueError) as exception:
            system.topology = topology
        assert str(
            exception.value) == "The number of atoms already found in the AtomicSystem object have to be equal to the number of atoms in the new topology"

    def test__eq__(self):
        system1 = AtomicSystem()
        system2 = AtomicSystem()
        system3 = AtomicSystem(pos=np.array(
            [[0, 0, 0], [1, 1, 1]]), atoms=[Atom('C'), Atom('H')])
        system4 = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                               atoms=[Atom('C'), Atom('H')], cell=Cell(0.75, 0.75, 0.75))
        system5 = AtomicSystem(pos=np.array(
            [[0, 0, 0], [1, 1, 1]]), atoms=[Atom('C'), Atom('D')])

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

    def test__getitem__(self):
        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                              atoms=[Atom('C'), Atom('H')], cell=Cell(0.75, 0.75, 0.75))

        assert system[0] == AtomicSystem(pos=np.array([[0, 0, 0]]),
                                         atoms=[Atom('C')], cell=Cell(0.75, 0.75, 0.75))
        assert system[1] == AtomicSystem(pos=np.array([[1, 1, 1]]),
                                         atoms=[Atom('H')], cell=Cell(0.75, 0.75, 0.75))

        with pytest.raises(IndexError):
            system[2]

        assert system[-1] == AtomicSystem(pos=np.array([[1, 1, 1]]),
                                          atoms=[Atom('H')], cell=Cell(0.75, 0.75, 0.75))

        system = AtomicSystem(vel=np.array([[0, 0, 0], [1, 1, 1]]), forces=np.array([[0, 0, 0], [1, 1, 1]]), charges=np.array([0, 0]),
                              atoms=[Atom('C'), Atom('H')], cell=Cell(0.75, 0.75, 0.75))

        assert system[0] == AtomicSystem(vel=np.array([[0, 0, 0]]), forces=np.array([[0, 0, 0]]), charges=np.array([0]),
                                         atoms=[Atom('C')], cell=Cell(0.75, 0.75, 0.75))

        assert system[Atom(6)] == AtomicSystem(vel=np.array([[0, 0, 0]]), forces=np.array([[0, 0, 0]]), charges=np.array([0]),
                                               atoms=[Atom('C')], cell=Cell(0.75, 0.75, 0.75))

        system = AtomicSystem(vel=np.array([[0, 0, 0], [1, 1, 1]]), forces=np.array([[0, 0, 0], [1, 1, 1]]), charges=np.array([0, 0]),
                              atoms=None, cell=Cell(0.75, 0.75, 0.75))

        assert system[:] == system

    def test_n_atoms(self):
        system = AtomicSystem()
        assert system.n_atoms == 0

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]))
        assert system.n_atoms == 2

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                              atoms=[Atom('C'), Atom('H')], cell=Cell(0.75, 0.75, 0.75))
        assert system.n_atoms == 2

        system = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                              atoms=[Atom('C')], cell=Cell(0.75, 0.75, 0.75))

        with pytest.raises(ValueError) as exception:
            system.n_atoms
        assert str(
            exception.value) == "The number of atoms (or atoms in the topology), positions, velocities, forces and charges must be equal."

    def test__str__(self):
        pos = np.array([[0, 0, 0], [1, 1, 1]])
        vel = np.array([[0, 0, 0], [2, 2, 2]])
        charges = np.array([0, 1])
        forces = np.array([[0, 0, 0], [3, 3, 3]])
        cell = Cell(0.75, 0.75, 0.75)
        atoms = [Atom('C'), Atom('H')]

        system = AtomicSystem()
        assert str(
            system) == f"AtomicSystem(topology=({system.topology}), cell=({system.cell}))"
        assert str(system) == repr(system)

        system = AtomicSystem(atoms=atoms, pos=pos, cell=cell,
                              vel=vel, charges=charges, forces=forces)
        assert str(
            system) == f"AtomicSystem(topology=({system.topology}), cell=({system.cell}))"
        assert str(system) == repr(system)