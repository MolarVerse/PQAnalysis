import pytest
import numpy as np

from beartype.roar import BeartypeException

from PQAnalysis.io.frameReader import FrameReader
from PQAnalysis.core.cell import Cell
from PQAnalysis.core.atom import Atom
from PQAnalysis.exceptions import TrajectoryFormatError
from PQAnalysis.traj.formats import TrajectoryFormat


class TestFrameReader:

    def test__read_header_line(self):
        reader = FrameReader()

        with pytest.raises(ValueError) as exception:
            reader._read_header_line("1 2.0 3.0")
        assert str(
            exception.value) == "Invalid file format in header line of Frame."

        n_atoms, cell = reader._read_header_line(
            "1 2.0 3.0 4.0 5.0 6.0 7.0")
        assert n_atoms == 1
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [5.0, 6.0, 7.0])

        n_atoms, cell = reader._read_header_line("2 2.0 3.0 4.0")
        assert n_atoms == 2
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [90.0, 90.0, 90.0])

        n_atoms, cell = reader._read_header_line("3")
        assert n_atoms == 3
        assert cell is None

    def test__read_xyz(self):
        reader = FrameReader()

        with pytest.raises(ValueError) as exception:
            reader._read_xyz(
                ["", "", "h 1.0 2.0 3.0", "o 2.0 2.0"], n_atoms=2)
        assert str(
            exception.value) == "Invalid file format in xyz coordinates of Frame."

        xyz, atoms = reader._read_xyz(
            ["", "", "h 1.0 2.0 3.0", "o 2.0 2.0 2.0"], n_atoms=2)
        assert np.allclose(xyz, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert atoms == ["h", "o"]

    def test__read_scalar(self):
        reader = FrameReader()

        with pytest.raises(ValueError) as exception:
            reader._read_scalar(["", "", "h 1.0 2.0 3.0"], n_atoms=1)
        assert str(
            exception.value) == "Invalid file format in scalar values of Frame."

        scalar, atoms = reader._read_scalar(["", "", "h 1.0"], n_atoms=1)
        assert np.allclose(scalar, [1.0])
        assert atoms == ["h"]

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

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0", format="vel")
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom, use_guess_element=False)
                               for atom in ["h", "o1"]]
        assert np.allclose(frame.vel, [
                           [1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0", format="force")
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom, use_guess_element=False)
                               for atom in ["h", "o1"]]
        assert np.allclose(frame.forces, [
                           [1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0\no1 2.0", format="charge")
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom, use_guess_element=False)
                               for atom in ["h", "o1"]]
        assert np.allclose(frame.charges, [1.0, 2.0])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

    def test_read_invalid_format(self):
        reader = FrameReader()

        with pytest.raises(TrajectoryFormatError) as exception:
            reader.read("", format="invalid")
        assert str(
            exception.value) == f"""
'invalid' is not a valid TrajectoryFormat.
Possible values are: {TrajectoryFormat.member_repr()}
or their case insensitive string representation: {TrajectoryFormat.value_repr()}"""
