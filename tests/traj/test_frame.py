import numpy as np
import pytest

from PQAnalysis.traj.frame import Frame
from PQAnalysis.core.cell import Cell
from PQAnalysis.core.atomicSystem import AtomicSystem
from PQAnalysis.core.atom import Atom
from PQAnalysis.core.topology import Topology


class TestFrame:

    def test__init__(self):
        coords = np.array([[0, 0, 0], [1, 0, 0]])
        atoms = [Atom(atom) for atom in ['C', 'H']]
        cell = Cell(10, 10, 10)
        system = AtomicSystem(pos=coords, atoms=atoms, cell=cell)

        frame = Frame(system)

        assert frame.n_atoms == 2
        assert np.allclose(frame.pos, coords)
        assert frame.atoms == atoms
        assert frame.cell == Cell(10, 10, 10)
        assert frame.PBC

    def test_compute_com(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        cell = Cell(10, 10, 10)
        atoms = [Atom(atom) for atom in ['C', 'H', 'H']]
        frame = Frame(AtomicSystem(pos=coords, atoms=atoms, cell=cell))

        com_frame = frame.compute_com_frame()
        assert isinstance(com_frame, Frame)
        assert com_frame.n_atoms == 1
        assert np.allclose(com_frame.pos, [[0.21557785, 0, 0]])
        assert com_frame.system.combined_name == 'CHH'

        with pytest.raises(ValueError) as exception:
            frame.compute_com_frame(group=2)
        assert str(
            exception.value) == 'Number of atoms in selection is not a multiple of group.'

        com_frame = frame.compute_com_frame(group=1)
        assert isinstance(com_frame, Frame)
        assert com_frame.n_atoms == 3
        assert np.allclose(com_frame.pos, [[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        assert com_frame.atoms == [
            Atom(atom, use_guess_element=False) for atom in ['C', 'H', 'H']]

    def test__eq__(self):
        frame1 = Frame(AtomicSystem(pos=np.array([[0, 0, 0]])))
        frame2 = Frame(AtomicSystem(pos=np.array([[0, 0, 0]])))
        frame3 = Frame(AtomicSystem(pos=np.array([[0, 0, 0]])))
        frame4 = Frame(AtomicSystem(pos=np.array([[0, 0, 1]])))

        assert frame1 == frame2
        assert frame1 == frame3
        assert frame1 != frame4
        assert frame2 != frame4

        assert frame1 != 1

    def test__getitem__(self):
        frame = Frame(AtomicSystem(pos=np.array([[0, 0, 0], [1, 1, 1]])))
        assert np.allclose(frame[0].pos, [[0, 0, 0]])
        assert np.allclose(frame[1].pos, [[1, 1, 1]])
        assert np.allclose(frame[-1].pos, [[1, 1, 1]])

        with pytest.raises(IndexError) as exception:
            frame[2]
