import pytest
import numpy as np

from PQAnalysis.core.atomicSystem import AtomicSystem
from PQAnalysis.core.atom import Atom
from PQAnalysis.core.cell import Cell


class TestAtomicSystem:
    def test__init__(self):
        system = AtomicSystem()
        assert system.atoms == []
        assert system.pos.shape == (0, 3)
        assert system.vel.shape == (0, 3)
        assert system.forces.shape == (0, 3)
        assert system.charges.shape == (0,)
        assert system.cell is None
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
        assert system.cell is None
        assert system.PBC is False
        assert system.n_atoms == 0
        assert np.allclose(system.pos, [[0, 0, 0], [1, 1, 1]])

        assert system.atomic_masses.shape == (0,)
        assert system.mass == 0

        with pytest.raises(ValueError) as exception:
            system.center_of_mass
        assert str(
            exception.value) == "AtomicSystem contains a different number of atoms to positions."

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
        with pytest.raises(ValueError) as exception:
            system.atomic_masses
        assert str(
            exception.value) == "AtomicSystem contains atoms without mass information."
        with pytest.raises(ValueError) as exception:
            system.mass
        assert str(
            exception.value) == "AtomicSystem contains atoms without mass information."
        with pytest.raises(ValueError) as exception:
            system.center_of_mass
        assert str(
            exception.value) == "AtomicSystem contains atoms without mass information."

    def test__eq__(self):
        system1 = AtomicSystem()
        system2 = AtomicSystem()
        system3 = AtomicSystem(pos=np.array(
            [[0, 0, 0], [1, 1, 1]]), atoms=[Atom('C'), Atom('H')])
        system4 = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]),
                               atoms=[Atom('C'), Atom('H')], cell=Cell(0.75, 0.75, 0.75))

        assert system1 == system2
        assert system1 != system3
        assert system1 != system4
        assert system3 != system4
        assert system4 == system4

        assert system1 != 1

        system1 = AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]]))
        system2 = AtomicSystem(pos=np.array([[0, 0, 0]]))

        assert system1 != system2

        system1 = AtomicSystem(vel=np.array([[0, 0, 0], [1, 1, 1]]))
        system2 = AtomicSystem(vel=np.array([[0, 0, 0]]))

        assert system1 != system2

        system1 = AtomicSystem(forces=np.array([[0, 0, 0], [1, 1, 1]]))
        system2 = AtomicSystem(forces=np.array([[0, 0, 0]]))

        assert system1 != system2

        system1 = AtomicSystem(charges=np.array([0, 1]))
        system2 = AtomicSystem(charges=np.array([0]))

        assert system1 != system2

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
