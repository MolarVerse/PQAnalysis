import pytest
import numpy as np

from PQAnalysis.topology.shake_topology import ShakeTopologyGenerator
from PQAnalysis.core import Atom
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.traj import Trajectory
from PQAnalysis.exceptions import PQValueError

from . import pytestmark  # pylint: disable=unused-import

#pylint: disable=protected-access



class TestShakeTopologyGenerator:

    def test__init__(self):
        generator = ShakeTopologyGenerator()
        assert generator.selection_object is None
        assert generator._use_full_atom_info is False

        generator = ShakeTopologyGenerator([Atom('C'), Atom('H')], True)
        assert generator.selection_object == [Atom('C'), Atom('H')]
        assert generator._use_full_atom_info is True

        generator = ShakeTopologyGenerator(['C', 'H'])
        assert generator.selection_object == ['C', 'H']
        assert generator._use_full_atom_info is False

        generator = ShakeTopologyGenerator(np.array([0, 1]))
        assert np.allclose(generator.selection_object, [0, 1])
        assert generator._use_full_atom_info is False

    def test_generate_topology(self):
        atoms = [Atom('C'), Atom('H'), Atom('H'), Atom('O'), Atom('H')]
        pos = np.array(
            [[0.1, 0, 0], [1, 0, 0], [2.1, 0, 0], [3, 0, 0], [4, 0, 0]]
        )
        pos2 = np.array(
            [[0.5, 0, 0], [1, 0.5, 0], [2.5, 0, 0], [3, 0.5, 0], [4.5, 0, 0]]
        )
        system = AtomicSystem(pos=pos, atoms=atoms)
        system2 = AtomicSystem(pos=pos2, atoms=atoms)

        traj = Trajectory([system, system2])

        generator = ShakeTopologyGenerator(selection=[Atom('H')])
        generator.generate_topology(traj)
        indices, target_indices, distances = generator.indices, generator.target_indices, generator.distances

        assert np.allclose(indices, [1, 2, 4])
        assert np.allclose(target_indices, [0, 3, 3])
        assert np.allclose(distances, [0.80355339, 0.80355339, 1.29056942])

    def test_average_equivalents(self):
        atoms = [Atom('C'), Atom('H'), Atom('H'), Atom('O'), Atom('H')]
        pos = np.array(
            [[0.1, 0, 0], [1, 0, 0], [2.1, 0, 0], [3, 0, 0], [4, 0, 0]]
        )

        pos2 = np.array(
            [[0.5, 0, 0], [1, 0.5, 0], [2.5, 0, 0], [3, 0.5, 0], [4.5, 0, 0]]
        )

        system = AtomicSystem(pos=pos, atoms=atoms)
        system2 = AtomicSystem(pos=pos2, atoms=atoms)

        traj = Trajectory([system, system2])

        generator = ShakeTopologyGenerator(selection=[Atom('H')])
        generator.generate_topology(traj)
        generator.average_equivalents([np.array([1]), np.array([2, 4])])

        indices, target_indices, distances = generator.indices, generator.target_indices, generator.distances

        assert np.allclose(indices, [1, 2, 4])
        assert np.allclose(target_indices, [0, 3, 3])
        assert np.allclose(distances, [0.80355339, 1.047061405, 1.047061405])
        assert generator.line_comments is None

        generator.average_equivalents(
            [np.array([1]), np.array([2, 4])],
            comments=['test', 'test2'],
        )
        assert generator.line_comments == ['test', 'test2', 'test2']

    def test_add_comments(self):
        generator = ShakeTopologyGenerator()
        with pytest.raises(PQValueError) as exception:
            generator.add_comments(['test', 'test2'])
        assert str(
            exception.value
        ) == "The number of comments does not match the number of indices."

        generator.indices = np.array([1, 2, 3])

        with pytest.raises(PQValueError) as exception:
            generator.add_comments(['test', 'test2'])
        assert str(
            exception.value
        ) == "The number of comments does not match the number of indices."

        generator.add_comments(['test', 'test2', 'test3'])
        assert generator.line_comments == ['test', 'test2', 'test3']

    def test_write_topology(self, capsys):
        atoms = [Atom('C'), Atom('H'), Atom('H'), Atom('O'), Atom('H')]
        pos = np.array(
            [[0.1, 0, 0], [1, 0, 0], [2.1, 0, 0], [3, 0, 0], [4, 0, 0]]
        )

        pos2 = np.array(
            [[0.5, 0, 0], [1, 0.5, 0], [2.5, 0, 0], [3, 0.5, 0], [4.5, 0, 0]]
        )

        system = AtomicSystem(pos=pos, atoms=atoms)
        system2 = AtomicSystem(pos=pos2, atoms=atoms)

        traj = Trajectory([system, system2])

        generator = ShakeTopologyGenerator(selection=[Atom('H')])
        generator.generate_topology(traj)

        print()
        generator.write_topology()

        captured = capsys.readouterr()
        assert captured.out == """
SHAKE 3 2 0
    2     1   0.803553390593\t
    3     4   0.803553390593\t
    5     4   1.290569415042\t
END
"""

        print()

        generator.add_comments(['test', 'test2', 'test3'])
        generator.write_topology()

        captured = capsys.readouterr()
        assert captured.out == """
SHAKE 3 2 0
    2     1   0.803553390593\t # test
    3     4   0.803553390593\t # test2
    5     4   1.290569415042\t # test3
END
"""

        print()

        generator.average_equivalents(
            [np.array([1]), np.array([2, 4])],
            comments=['test', 'test2'],
        )
        generator.write_topology()

        captured = capsys.readouterr()
        assert captured.out == """
SHAKE 3 2 0
    2     1   0.803553390593\t # test
    3     4   1.047061402818\t # test2
    5     4   1.047061402818\t # test2
END
"""
