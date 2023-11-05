import numpy as np
import pytest

from PQAnalysis.traj.frame import Frame, MolecularFrame
from PQAnalysis.core.cell import Cell
from PQAnalysis.atomicUnits.element import Elements
from PQAnalysis.coordinates.coordinates import Coordinates
from PQAnalysis.atomicUnits.molecule import Molecule


class TestFrame:

    def test__init__(self):
        coords = Coordinates([[0, 0, 0], [1, 0, 0]])
        elements = Elements(['C', 'H'])
        cell = Cell(10, 10, 10)

        frame = Frame(coords, elements)

        assert frame.n_atoms == 2
        assert frame.coordinates == coords
        assert frame.atoms == elements
        assert frame.cell is None

        coords.cell = cell
        frame = Frame(coords, elements)

        assert frame.n_atoms == 2
        assert frame.coordinates == coords
        assert frame.atoms == elements
        assert frame.cell == Cell(10, 10, 10)

        frame = Frame(coordinates=[[0, 0, 0]], atoms='C')
        assert frame.n_atoms == 1
        assert frame.coordinates == Coordinates([[0, 0, 0]])
        assert frame.atoms == Elements(['C'])
        assert frame.cell is None

        with pytest.raises(ValueError) as exception:
            Frame(coordinates=[[0, 0, 0], [1, 0, 0]], atoms=['C'])
        assert str(
            exception.value) == 'coordinates and atoms must have the same length.'

        frame = Frame(coordinates=[[0, 0, 0], [1, 0, 0]], atoms=[1, 2])
        assert frame.n_atoms == 2
        assert frame.coordinates == Coordinates([[0, 0, 0], [1, 0, 0]])
        assert frame.atoms == Elements(['H', 'He'])
        assert frame.cell is None

    def test_PBC(self):
        coords = Coordinates([[0, 0, 0], [1, 0, 0]])

        frame = Frame(coordinates=coords, atoms=['C', 'H'])
        assert not frame.PBC

        coords.cell = Cell(10, 10, 10)
        frame = Frame(coordinates=coords, atoms=['C', 'H'])
        assert frame.PBC

    def test__getitem__(self):
        frame = Frame(Coordinates([[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
                      Elements(['C', 'H', 'H']))

        assert frame[0].n_atoms == 1
        assert frame[0].coordinates == Coordinates([[0, 0, 0]])
        assert frame[0].atoms == Elements(['C'])
        assert frame[0].cell is None

    def test_compute_com(self):
        frame = Frame(Coordinates([[0, 0, 0], [1, 0, 0], [2, 0, 0]], Cell(10, 10, 10)),
                      ['C', 'H', 'H'])

        com_frame = frame.compute_com_frame()
        assert isinstance(com_frame, MolecularFrame)
        assert com_frame.n_molecules == 1
        assert com_frame.coordinates == Coordinates(
            [[0.21557785, 0, 0]], Cell(10, 10, 10))
        assert com_frame.molecules[0].name == 'chh'

        with pytest.raises(ValueError) as exception:
            frame.compute_com_frame(group=2)
        assert str(
            exception.value) == 'Number of atoms in selection is not a multiple of group.'

        com_frame = frame.compute_com_frame(group=1)
        assert isinstance(com_frame, MolecularFrame)
        assert com_frame.n_molecules == 3
        assert com_frame.coordinates == Coordinates(
            [[0, 0, 0], [1, 0, 0], [2, 0, 0]], Cell(10, 10, 10))
        assert com_frame.molecules[0].name == 'c'
        assert com_frame.molecules[1].name == 'h'
        assert com_frame.molecules[2].name == 'h'

    def test__eq__(self):
        frame1 = Frame(coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                       atoms=['C', 'H', 'H'])
        frame2 = Frame(coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                       atoms=['C', 'H', 'H'])
        assert frame1 == frame2
        frame2 = Frame(coordinates=Coordinates([[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
                       atoms=['C', 'H', 'H'])
        frame2.cell = Cell(10, 10, 10)
        assert frame1 != frame2

        frame2 = Frame(coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 1]],
                       atoms=['C', 'H', 'H'])
        assert frame1 != frame2

        frame2 = Frame(coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                       atoms=['C', 'H', 'O'])
        assert frame1 != frame2

        assert frame1 != 1

        frame2 = Frame(coordinates=[[0, 0, 0], [1, 0, 0]],
                       atoms=['C', 'H'])

        assert frame1 != frame2

    def test_is_combinable(self):
        frame1 = Frame([[0, 0, 0], [1, 0, 0], [2, 0, 0]], ['C', 'H', 'H'])
        frame2 = Frame([[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                       ['C', 'H', 'H'])
        frame3 = Frame([[0, 0, 0], [1, 0, 0]], ['C', 'H'])
        frame4 = Frame([[0, 0, 0], [1, 0, 0]], ['C', 'O'])

        assert frame1.is_combinable(frame2)
        assert not frame1.is_combinable(frame3)
        assert not frame1.is_combinable(frame4)
        assert not frame1.is_combinable(1)

        assert not frame3.is_combinable(frame4)


class TestMolecularFrame:
    def test__init__(self):
        frame = MolecularFrame(molecules=[Molecule(atoms=['C', 'H', 'H'],
                                                   coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])],
                               cell=Cell(10, 10, 10))

        assert frame.n_molecules == 1
        assert frame.molecules[0].name == 'chh'
        assert frame.cell == Cell(10, 10, 10)

        frame = MolecularFrame(molecules=[Molecule(atoms=['C', 'H', 'H'],
                                                   coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])])

        assert frame.n_molecules == 1
        assert frame.molecules[0].name == 'chh'
        assert frame.cell is None

        with pytest.raises(ValueError) as exception:
            MolecularFrame(molecules=[Molecule(atoms=['C', 'H', 'H'],
                                               coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])],
                           coordinates=Coordinates([[0, 0, 0], [1, 0, 0], [2, 0, 0]]))
