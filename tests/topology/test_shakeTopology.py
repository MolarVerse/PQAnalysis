import numpy as np

from PQAnalysis.topology.shakeTopology import ShakeTopologyGenerator
from PQAnalysis.core.atom import Atom
from PQAnalysis.core.atomicSystem import AtomicSystem
from PQAnalysis.traj.frame import Frame
from PQAnalysis.traj.trajectory import Trajectory


class TestShakeTopologyGenerator:
    def test__init__(self):
        generator = ShakeTopologyGenerator()
        assert generator._atoms is None
        assert generator._use_full_atom_info is False

        generator = ShakeTopologyGenerator([Atom('C'), Atom('H')], True)
        assert generator._atoms == [Atom('C'), Atom('H')]
        assert generator._use_full_atom_info is True

        generator = ShakeTopologyGenerator(['C', 'H'])
        assert generator._atoms == [Atom('C'), Atom('H')]
        assert generator._use_full_atom_info is False

        generator = ShakeTopologyGenerator(np.array([0, 1]))
        assert np.allclose(generator._atoms, [0, 1])
        assert generator._use_full_atom_info is False

    def test_generate_topology(self):
        atoms = [Atom('C'), Atom('H'), Atom('H'), Atom('O'), Atom('H')]
        pos = np.array([[0.1, 0, 0], [1, 0, 0], [2.1, 0, 0],
                        [3, 0, 0], [4, 0, 0]])
        pos2 = np.array([[0.5, 0, 0], [1, 0.5, 0], [2.5, 0, 0],
                         [3, 0.5, 0], [4.5, 0, 0]])
        system = AtomicSystem(pos=pos, atoms=atoms)
        system2 = AtomicSystem(pos=pos2, atoms=atoms)

        traj = Trajectory([Frame(system), Frame(system2)])

        generator = ShakeTopologyGenerator(atoms=[Atom('H')])
        indices, target_indices, distances = generator.generate_topology(traj)

        assert np.allclose(indices, [[1], [2], [4]])
        assert np.allclose(target_indices, [[0], [3], [3]])
        assert np.allclose(distances, [[1.0], [1.0], [1.0]])
