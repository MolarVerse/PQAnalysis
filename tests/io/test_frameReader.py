import pytest
import numpy as np

from PQAnalysis.io import _FrameReader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError
from PQAnalysis.core import Cell, Atom
from PQAnalysis.traj.exceptions import TrajectoryFormatError
from PQAnalysis.traj import TrajectoryFormat
from PQAnalysis.topology import Topology

from . import pytestmark



class TestFrameReader:

    def test__read_header_line(self):
        reader = _FrameReader()

        with pytest.raises(FrameReaderError) as exception:
            reader._read_header_line("1 2.0 3.0")
        assert str(
            exception.value
        ) == "Invalid file format in header line of Frame."

        n_atoms, cell = reader._read_header_line("1 2.0 3.0 4.0 5.0 6.0 7.0")
        assert n_atoms == 1
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [5.0, 6.0, 7.0])

        n_atoms, cell = reader._read_header_line("2 2.0 3.0 4.0")
        assert n_atoms == 2
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [90.0, 90.0, 90.0])

        n_atoms, cell = reader._read_header_line("3")
        assert n_atoms == 3
        assert cell == Cell()

    def test__read_xyz(self):
        reader = _FrameReader()

        with pytest.raises(FrameReaderError) as exception:
            reader._read_xyz(["", "", "h 1.0 2.0 3.0", "o 2.0 2.0"], n_atoms=2)
        assert str(
            exception.value
        ) == "Invalid file format in xyz coordinates of Frame."

        xyz, atoms = reader._read_xyz(
            ["", "", "h 1.0 2.0 3.0", "o 2.0 2.0 2.0"], n_atoms=2)
        assert np.allclose(xyz, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert atoms == ["h", "o"]

    def test__read_scalar(self):
        reader = _FrameReader()

        with pytest.raises(FrameReaderError) as exception:
            reader._read_scalar(["", "", "h 1.0 2.0 3.0"], n_atoms=1)
        assert str(
            exception.value
        ) == "Invalid file format in scalar values of Frame."

        scalar, atoms = reader._read_scalar(["", "", "h 1.0"], n_atoms=1)
        assert np.allclose(scalar, [1.0])
        assert atoms == ["h"]

    def test_read(self):
        reader = _FrameReader()

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no 2.0 2.0 2.0"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom) for atom in ["h", "o"]]
        assert np.allclose(frame.pos, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.pos, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0",
            traj_format="vel"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.vel, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0",
            traj_format="force"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.forces, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0\no1 2.0",
            traj_format="charge"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.charges, [1.0, 2.0])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

    def test_read_invalid_format(self):
        reader = _FrameReader()

        with pytest.raises(TrajectoryFormatError) as exception:
            reader.read("", traj_format="invalid")
        assert str(exception.value) == (
            "\n"
            "'invalid' is not a valid TrajectoryFormat.\n"
            f"Possible values are: {TrajectoryFormat.member_repr()} "
            "or their case insensitive string representation: "
            f"{TrajectoryFormat.value_repr()}"
        )

    def test__get_topology(self):
        reader = _FrameReader()

        topology = reader._get_topology(["h", "o"], None)
        assert topology.atoms == [Atom(atom) for atom in ["h", "o"]]

        topology = Topology()
        topology = reader._get_topology(["h", "o"], topology)
        assert topology == Topology()
