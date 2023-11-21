import numpy as np
import pytest

from PQAnalysis.core.atomicSystem.atomicSystem import AtomicSystem
from PQAnalysis.core.atom import Atom


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
            n=1, atoms=[Atom('H1', 1)], use_full_atom_info=True)
        assert np.allclose(indices, [[2], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1, atoms=['H1'])
        assert np.allclose(indices, [[2], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1, atoms=[Atom('H1', 1)], use_full_atom_info=False)
        assert np.allclose(indices, [[2], [3], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [np.sqrt(2.0)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1, atoms=np.array([0, 2]))
        assert np.allclose(indices, [[3], [3]])
        assert np.allclose(distances, [[1.0], [np.sqrt(2)]])

        with pytest.raises(ValueError) as exception:
            system.nearest_neighbours(
                n=1, atoms=['H1'], use_full_atom_info=True)
        assert str(
            exception.value) == "use_full_atom_info can only be used with List[Atom]"

        system = AtomicSystem(atoms=atoms)
        with pytest.raises(ValueError) as exception:
            system.nearest_neighbours(n=1)
        assert str(
            exception.value) == "AtomicSystem contains a different number of atoms to positions."
