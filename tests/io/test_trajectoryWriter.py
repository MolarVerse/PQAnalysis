import pytest
import sys
import numpy as np

from . import pytestmark

from PQAnalysis.io import TrajectoryWriter, write_trajectory, FileWritingMode
from PQAnalysis.traj import Trajectory, TrajectoryFormat, MDEngineFormat
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.traj.exceptions import MDEngineFormatError



def test_write_trajectory(capsys):
    atoms = [Atom(atom) for atom in ['h', 'o']]
    coordinates1 = np.array([[0, 0, 0], [0, 0, 1]])
    coordinates2 = np.array([[0, 0, 0], [0, 0, 1]])
    frame1 = AtomicSystem(atoms=atoms, pos=coordinates1, cell=Cell(10, 10, 10))
    frame2 = AtomicSystem(atoms=atoms, pos=coordinates2, cell=Cell(11, 10, 10))
    traj = Trajectory([frame1, frame2])

    print()
    write_trajectory(traj, engine_format="PQ")

    captured = capsys.readouterr()
    assert captured.out == """
2 10 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
2 11 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""



class TestTrajectoryWriter:

    def test__init__(self):

        with pytest.raises(MDEngineFormatError) as exception:
            TrajectoryWriter(engine_format="notAFormat")
        assert str(exception.value) == (
            "\n"
            "'notaformat' is not a valid MDEngineFormat.\n"
            f"Possible values are: {MDEngineFormat.member_repr()} "
            "or their case insensitive string representation: "
            f"{MDEngineFormat.value_repr()}"
        )

        writer = TrajectoryWriter()
        assert writer.file == sys.stdout
        assert writer.filename is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.format == MDEngineFormat.PQ

        writer = TrajectoryWriter(engine_format="qmcfc")
        assert writer.format == MDEngineFormat.QMCFC

        writer = TrajectoryWriter(engine_format="PQ")
        assert writer.format == MDEngineFormat.PQ

    def test__write_header(self, capsys):

        writer = TrajectoryWriter()
        assert writer.mode == FileWritingMode.WRITE
        writer._write_header(1, Cell(10, 10, 10))

        captured = capsys.readouterr()
        assert captured.out == "1 10 10 10 90 90 90\n"

        writer._write_header(1)
        captured = capsys.readouterr()
        assert captured.out == "1\n"

    def test__write_comment(self, capsys):

        writer = TrajectoryWriter()
        writer.type = TrajectoryFormat.XYZ
        writer._write_comment(
            AtomicSystem(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            cell=Cell(10,
            10,
            10)
            )
        )

        captured = capsys.readouterr()
        assert captured.out == "\n"

        forces = np.array([[1, 0, 3], [0, 2, 1]])
        writer.type = TrajectoryFormat.FORCE
        writer._write_comment(
            AtomicSystem(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            cell=Cell(10,
            10,
            10),
            forces=forces
            )
        )

        captured = capsys.readouterr()
        assert captured.out == "sum of forces: 1.000000e+00 2.000000e+00 4.000000e+00\n"

    def test__write_xyz(self, capsys):

        writer = TrajectoryWriter()
        writer.type = TrajectoryFormat.XYZ

        print()
        writer._write_xyz(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            xyz=np.array([[0,
            0,
            0],
            [0,
            0,
            1]])
        )

        captured = capsys.readouterr()
        assert captured.out == """
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        writer.format = "qmcfc"

        print()
        writer._write_xyz(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            xyz=np.array([[0,
            0,
            0],
            [0,
            0,
            1]])
        )

        captured = capsys.readouterr()
        assert captured.out == """
X   0.0 0.0 0.0
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

    def test__write_scalar(self, capsys):

        writer = TrajectoryWriter()
        writer._write_scalar(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            scalar=np.array([1,
            2])
        )

        captured = capsys.readouterr()
        assert captured.out == "h 1\no 2\n"

    def test_write(self, capsys):

        atoms = [Atom(atom) for atom in ['h', 'o']]
        coordinates1 = np.array([[0, 0, 0], [0, 0, 1]])
        coordinates2 = np.array([[0, 0, 0], [0, 0, 1]])

        frame1 = AtomicSystem(
            atoms=atoms,
            pos=coordinates1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            pos=coordinates2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()
        assert writer.mode == FileWritingMode.WRITE

        print()
        writer.write(traj)
        assert writer.mode == FileWritingMode.APPEND

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
2 11 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        writer = TrajectoryWriter()
        print()
        writer.write(frame1)

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        frame1 = AtomicSystem(
            atoms=atoms,
            vel=coordinates1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            vel=coordinates2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        print()
        writer.write(traj, traj_type="vel")

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00
o 0.000000000000e+00 0.000000000000e+00 1.000000000000e+00
2 11 10 10 90 90 90

h 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00
o 0.000000000000e+00 0.000000000000e+00 1.000000000000e+00
"""

        frame1 = AtomicSystem(
            atoms=atoms,
            forces=coordinates1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            forces=coordinates2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        print()
        writer.write(traj, traj_type="force")

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90
sum of forces: 0.000000e+00 0.000000e+00 1.000000e+00
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
2 11 10 10 90 90 90
sum of forces: 0.000000e+00 0.000000e+00 1.000000e+00
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        charges1 = np.array([1, 2])
        charges2 = np.array([3, 4])

        frame1 = AtomicSystem(
            atoms=atoms,
            charges=charges1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            charges=charges2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        print()
        writer.write(traj, traj_type="charge")

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h 1
o 2
2 11 10 10 90 90 90

h 3
o 4
"""
