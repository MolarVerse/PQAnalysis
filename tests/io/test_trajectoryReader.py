import pytest
import numpy as np

from beartype.roar import BeartypeException

from PQAnalysis.io.trajectoryReader import TrajectoryReader, FrameReader
from PQAnalysis.core.cell import Cell
from PQAnalysis.traj.frame import Frame
from PQAnalysis.core.atomicSystem import AtomicSystem
from PQAnalysis.core.atom import Atom


class TestTrajectoryReader:
    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):
        with pytest.raises(FileNotFoundError) as exception:
            TrajectoryReader("tmp")
        assert str(exception.value) == "File tmp not found."

        open("tmp", "w")
        reader = TrajectoryReader("tmp")
        assert reader.filename == "tmp"
        assert reader.frames == []

    @pytest.mark.usefixtures("tmpdir")
    def test_read(self):

        file = open("tmp", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp")

        traj = reader.read()

        cell = Cell(1.0, 1.0, 1.0)
        atoms = [Atom(atom) for atom in ["h", "o"]]

        frame1 = Frame(system=AtomicSystem(
            atoms=atoms, pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), cell=cell))
        frame2 = Frame(system=AtomicSystem(
            atoms=atoms, pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]), cell=cell))

        assert traj[0] == frame1
        # NOTE: here cell is not none because of the consecutive reading of frames
        # Cell will be taken from the previous frame
        assert traj[1] == frame2


class TestFrameReader:

    def test__read_header_line__(self):
        reader = FrameReader()

        with pytest.raises(ValueError) as exception:
            reader.__read_header_line__("1 2.0 3.0")
        assert str(
            exception.value) == "Invalid file format in header line of Frame."

        n_atoms, cell = reader.__read_header_line__(
            "1 2.0 3.0 4.0 5.0 6.0 7.0")
        assert n_atoms == 1
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [5.0, 6.0, 7.0])

        n_atoms, cell = reader.__read_header_line__("2 2.0 3.0 4.0")
        assert n_atoms == 2
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [90.0, 90.0, 90.0])

        n_atoms, cell = reader.__read_header_line__("3")
        assert n_atoms == 3
        assert cell is None

    def test__read_xyz__(self):
        reader = FrameReader()

        with pytest.raises(ValueError) as exception:
            reader.__read_xyz__(
                ["", "", "h 1.0 2.0 3.0", "o 2.0 2.0"], n_atoms=2)
        assert str(
            exception.value) == "Invalid file format in xyz coordinates of Frame."

        xyz, atoms = reader.__read_xyz__(
            ["", "", "h 1.0 2.0 3.0", "o 2.0 2.0 2.0"], n_atoms=2)
        assert np.allclose(xyz, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert atoms == ["h", "o"]

    def test_read(self):
        reader = FrameReader()

        with pytest.raises(BeartypeException):
            reader.read(["tmp"])

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no 2.0 2.0 2.0")
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom) for atom in ["h", "o"]]
        assert np.allclose(frame.pos, [
                           [1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0")
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom, use_guess_element=False)
                               for atom in ["h", "o1"]]
        assert np.allclose(frame.pos, [
                           [1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
