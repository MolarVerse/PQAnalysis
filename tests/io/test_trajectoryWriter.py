import pytest
import sys
import numpy as np

from PQAnalysis.io.trajectoryWriter import TrajectoryWriter, write_trajectory
from PQAnalysis.traj.frame import Frame
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.traj.formats import TrajectoryFormat, MDEngineFormat
from PQAnalysis.core.cell import Cell
from PQAnalysis.core.atomicSystem import AtomicSystem
from PQAnalysis.core.atom import Atom
from PQAnalysis.exceptions import MDEngineFormatError

# TODO: here only one option is tested - think of a better way to test all options


def test_write_trajectory(capsys):
    atoms = [Atom(atom) for atom in ['h', 'o']]
    coordinates1 = np.array([[0, 0, 0], [0, 0, 1]])
    coordinates2 = np.array([[0, 0, 0], [0, 0, 1]])
    frame1 = Frame(AtomicSystem(
        atoms=atoms, pos=coordinates1, cell=Cell(10, 10, 10)))
    frame2 = Frame(AtomicSystem(
        atoms=atoms, pos=coordinates2, cell=Cell(11, 10, 10)))
    traj = Trajectory([frame1, frame2])

    write_trajectory(traj, format="pimd-qmcf")

    captured = capsys.readouterr()
    assert captured.out == "2 10 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n2 11 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n"


class TestTrajectoryWriter:
    def test__init__(self):

        with pytest.raises(MDEngineFormatError) as exception:
            TrajectoryWriter(format="notAFormat")
        assert str(
            exception.value) == f"""
'notaformat' is not a valid MDEngineFormat.
Possible values are: {MDEngineFormat.member_repr()}
or their case insensitive string representation: {MDEngineFormat.value_repr()}"""

        writer = TrajectoryWriter()
        assert writer.file == sys.stdout
        assert writer.filename is None
        assert writer.mode == "a"
        assert writer.format == MDEngineFormat.PIMD_QMCF

        writer = TrajectoryWriter(format="qmcfc")
        assert writer.format == MDEngineFormat.QMCFC

        writer = TrajectoryWriter(format="pimd-qmcf")
        assert writer.format == MDEngineFormat.PIMD_QMCF

    def test__write_header(self, capsys):

        writer = TrajectoryWriter()
        writer._write_header(1, Cell(10, 10, 10))

        captured = capsys.readouterr()
        assert captured.out == "1 10 10 10 90 90 90\n"

        writer._write_header(1)
        captured = capsys.readouterr()
        assert captured.out == "1\n"

    def test__write_comment(self, capsys):

        writer = TrajectoryWriter()
        writer._write_comment(Frame(AtomicSystem(
            atoms=[Atom(atom) for atom in ["h", "o"]], cell=Cell(10, 10, 10))))

        captured = capsys.readouterr()
        assert captured.out == "\n"

        forces = np.array([[1, 0, 3], [0, 2, 1]])
        writer._type = TrajectoryFormat.FORCE
        writer._write_comment(Frame(AtomicSystem(
            atoms=[Atom(atom) for atom in ["h", "o"]], cell=Cell(10, 10, 10), forces=forces)))

        captured = capsys.readouterr()
        assert captured.out == "sum of forces: 1 2 4\n"

    def test__write_xyz(self, capsys):

        writer = TrajectoryWriter()
        writer._write_xyz(
            atoms=[Atom(atom) for atom in ["h", "o"]], xyz=np.array([[0, 0, 0], [0, 0, 1]]))

        captured = capsys.readouterr()
        assert captured.out == "h 0 0 0\no 0 0 1\n"

        writer.format = "qmcfc"
        writer._write_xyz(
            atoms=[Atom(atom) for atom in ["h", "o"]], xyz=np.array([[0, 0, 0], [0, 0, 1]]))

        captured = capsys.readouterr()
        assert captured.out == "X   0.0 0.0 0.0\nh 0 0 0\no 0 0 1\n"

    def test__write_scalar(self, capsys):

        writer = TrajectoryWriter()
        writer._write_scalar(
            atoms=[Atom(atom) for atom in ["h", "o"]], scalar=np.array([1, 2]))

        captured = capsys.readouterr()
        assert captured.out == "h 1\no 2\n"

    def test_write(self, capsys):

        atoms = [Atom(atom) for atom in ['h', 'o']]
        coordinates1 = np.array([[0, 0, 0], [0, 0, 1]])
        coordinates2 = np.array([[0, 0, 0], [0, 0, 1]])

        frame1 = Frame(AtomicSystem(
            atoms=atoms, pos=coordinates1, cell=Cell(10, 10, 10)))
        frame2 = Frame(AtomicSystem(
            atoms=atoms, pos=coordinates2, cell=Cell(11, 10, 10)))

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        writer.write(traj)

        captured = capsys.readouterr()
        assert captured.out == "2 10 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n2 11 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n"

        frame1 = Frame(AtomicSystem(
            atoms=atoms, vel=coordinates1, cell=Cell(10, 10, 10)))
        frame2 = Frame(AtomicSystem(
            atoms=atoms, vel=coordinates2, cell=Cell(11, 10, 10)))

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        writer.write(traj, type="vel")

        captured = capsys.readouterr()
        assert captured.out == "2 10 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n2 11 10 10 90 90 90\n\nh 0 0 0\no 0 0 1\n"

        frame1 = Frame(AtomicSystem(
            atoms=atoms, forces=coordinates1, cell=Cell(10, 10, 10)))
        frame2 = Frame(AtomicSystem(
            atoms=atoms, forces=coordinates2, cell=Cell(11, 10, 10)))

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        writer.write(traj, type="force")

        captured = capsys.readouterr()
        assert captured.out == "2 10 10 10 90 90 90\nsum of forces: 0 0 1\nh 0 0 0\no 0 0 1\n2 11 10 10 90 90 90\nsum of forces: 0 0 1\nh 0 0 0\no 0 0 1\n"

        charges1 = np.array([1, 2])
        charges2 = np.array([3, 4])

        frame1 = Frame(AtomicSystem(
            atoms=atoms, charges=charges1, cell=Cell(10, 10, 10)))
        frame2 = Frame(AtomicSystem(
            atoms=atoms, charges=charges2, cell=Cell(11, 10, 10)))

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        writer.write(traj, type="charge")

        captured = capsys.readouterr()
        assert captured.out == "2 10 10 10 90 90 90\n\nh 1\no 2\n2 11 10 10 90 90 90\n\nh 3\no 4\n"
