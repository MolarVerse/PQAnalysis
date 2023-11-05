import pytest
import sys

from PQAnalysis.io.trajectoryWriter import TrajectoryWriter, write_trajectory
from PQAnalysis.traj.frame import Frame
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.core.cell import Cell
from PQAnalysis.coordinates.coordinates import Coordinates
from PQAnalysis.atomicUnits.element import Elements


def test_write_trajectory(capsys):
    coordinates1 = Coordinates([[0, 0, 0], [0, 0, 1]], cell=Cell(10, 10, 10))
    atoms = Elements(['h', 'o'])
    coordinates2 = Coordinates([[0, 0, 0], [0, 0, 1]], cell=Cell(11, 10, 10))
    traj = Trajectory([Frame(coordinates=coordinates1, atoms=atoms), Frame(
        coordinates=coordinates2, atoms=atoms)])

    write_trajectory(traj, format="pimd-qmcf")

    captured = capsys.readouterr()
    assert captured.out == "2 10 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n2 11 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n"


class TestTrajectoryWriter:
    def test__init__(self):

        with pytest.raises(ValueError) as exception:
            TrajectoryWriter(format="notAFormat")
        assert str(
            exception.value) == "Invalid format. Has to be either \'pimd-qmcf\', \'qmcfc\' or \'None\'."

        writer = TrajectoryWriter()
        assert writer.file == sys.stdout
        assert writer.filename is None
        assert writer.mode == "a"
        assert writer.format == "pimd-qmcf"

        writer = TrajectoryWriter(format="qmcfc")
        assert writer.format == "qmcfc"

        writer = TrajectoryWriter(format="pimd-qmcf")
        assert writer.format == "pimd-qmcf"

    def test__write_header__(self, capsys):

        writer = TrajectoryWriter()
        writer.__write_header__(1, Cell(10, 10, 10))

        captured = capsys.readouterr()
        assert captured.out == "1 10 10 10 90 90 90\n\n"

        writer.__write_header__(1)
        captured = capsys.readouterr()
        assert captured.out == "1\n\n"

    def test__write_coordinates__(self, capsys):

        writer = TrajectoryWriter()
        writer.__write_coordinates__(
            atoms=Elements(["h", "o"]), xyz=[[0, 0, 0], [0, 0, 1]])

        captured = capsys.readouterr()
        assert captured.out == "h 0 0 0\no 0 0 1\n"

        writer.format = "qmcfc"
        writer.__write_coordinates__(
            atoms=Elements(["h", "o"]), xyz=[[0, 0, 0], [0, 0, 1]])

        captured = capsys.readouterr()
        assert captured.out == "X   0.0 0.0 0.0\nh 0 0 0\no 0 0 1\n"

    def test_write(self, capsys):

        coordinates1 = Coordinates(
            [[0, 0, 0], [0, 0, 1]], cell=Cell(10, 10, 10))
        atoms = Elements(['h', 'o'])
        coordinates2 = Coordinates(
            [[0, 0, 0], [0, 0, 1]], cell=Cell(11, 10, 10))
        traj = Trajectory([Frame(coordinates=coordinates1, atoms=atoms), Frame(
            coordinates=coordinates2, atoms=atoms)])
        writer = TrajectoryWriter()

        writer.write(traj)

        captured = capsys.readouterr()
        assert captured.out == "2 10 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n2 11 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n"
