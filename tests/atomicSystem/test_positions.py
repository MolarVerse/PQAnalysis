import numpy as np
import pytest

from . import pytestmark

from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.atomicSystem.exceptions import AtomicSystemPositionsError
from PQAnalysis.core import Atom, Cell


class TestPositionsMixin:
    def test_nearest_neighbours(self):
        positions = np.array([[0, 0, 0], [10, 0, 0], [2, 1, 0], [1, 0, 0]])
        atoms = [Atom('C'), Atom('H1', 1), Atom(1), Atom('H1', 1)]

        system = AtomicSystem(pos=positions, atoms=atoms)

        indices, distances = system.nearest_neighbours()
        assert np.allclose(indices, [[3], [2], [3], [0]])
        assert np.allclose(
            distances, [[1.0], [np.sqrt(8*8+1*1)], [np.sqrt(2)], [1.0]])

        indices, distances = system._nearest_neighbours()
        assert np.allclose(indices, [[3], [2], [3], [0]])
        assert np.allclose(
            distances, [[1.0], [np.sqrt(8*8+1*1)], [np.sqrt(2)], [1.0]])

        indices, distances = system.nearest_neighbours(n=2)
        assert np.allclose(indices, [[3, 2], [2, 3], [3, 0], [0, 2]])
        assert np.allclose(
            distances, [[1.0, np.sqrt(5)], [np.sqrt(8*8+1), np.sqrt(9*9)], [np.sqrt(2), np.sqrt(5)], [1.0, np.sqrt(2)]])

        indices, distances = system.nearest_neighbours(
            n=1, selection=[Atom('H1', 1)], use_full_atom_info=True)
        assert np.allclose(indices, [[2], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1, selection=['H1'])
        assert np.allclose(indices, [[2], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1, selection=[Atom('H1', 1)], use_full_atom_info=False)
        assert np.allclose(indices, [[2], [3], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [np.sqrt(2.0)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1, selection=np.array([0, 2]))
        assert np.allclose(indices, [[3], [3]])
        assert np.allclose(distances, [[1.0], [np.sqrt(2)]])

        system = AtomicSystem(atoms=atoms)
        with pytest.raises(AtomicSystemPositionsError) as exception:
            system.nearest_neighbours(n=1)
        assert str(
            exception.value) == AtomicSystemPositionsError.message

    def test_image(self):
        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H1', 1)],
            pos=np.array([[0, 0, 0], [10, 0, 0]]),
            cell=Cell(8, 8, 8)
        )
        system.image()

        assert np.allclose(system.pos, np.array([[0, 0, 0], [2, 0, 0]]))

    def test_center(self):
        system = AtomicSystem(
            atoms=[Atom('C'), Atom('H1', 1)],
            pos=np.array([[0, 0, 0], [10, 0, 0]]),
            cell=Cell(8, 8, 8)
        )
        system.center(np.array([1, 0, 0]), image=False)

        assert np.allclose(system.pos, np.array([[-1, 0, 0], [9, 0, 0]]))

        system.center(np.array([1, 0, 0]), image=True)

        assert np.allclose(system.pos, np.array([[-2, 0, 0], [0, 0, 0]]))
