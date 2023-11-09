import pytest
import numpy as np

from _pytest.capture import CaptureFixture

from PQAnalysis.io.boxWriter import BoxWriter, write_box
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.traj.frame import Frame
from PQAnalysis.core.cell import Cell
from PQAnalysis.core.atomicSystem import AtomicSystem
from PQAnalysis.core.atom import Atom


class TestBoxWriter:

    def test__init__(self):
        with pytest.raises(ValueError) as exception:
            BoxWriter(filename="tmp", format="r")
        assert str(
            exception.value) == "Invalid format. Has to be either \'vmd\', \'data\' or \'None\'."

        writer = BoxWriter(filename="tmp", format="vmd")
        assert writer.file is None
        assert writer.mode == "a"
        assert writer.filename == "tmp"
        assert writer.format == "vmd"

        writer = BoxWriter(filename="tmp", format="data")
        assert writer.file is None
        assert writer.mode == "a"
        assert writer.filename == "tmp"
        assert writer.format == "data"

        writer = BoxWriter(filename="tmp")
        assert writer.file is None
        assert writer.mode == "a"
        assert writer.filename == "tmp"
        assert writer.format == "data"

    atoms1 = [Atom("H")]
    pos1 = np.array([[0, 1, 2]])
    cell1 = Cell(10, 10, 10)

    def test__check_PBC__(self):
        system1 = AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=self.cell1)
        system2 = AtomicSystem(atoms=self.atoms1, pos=self.pos1)

        frame1 = Frame(system1)
        frame2 = Frame(system2)

        traj1 = Trajectory([frame1, frame1])
        traj2 = Trajectory([frame1, frame2])

        writer = BoxWriter()

        try:
            writer.__check_PBC__(traj1)
        except:
            assert False

        with pytest.raises(ValueError) as exception:
            writer.__check_PBC__(traj2)
        assert str(
            exception.value) == "At least on cell of the trajectory is None. Cannot write box file."

    def test_write_box_file(self, capsys: CaptureFixture):
        writer = BoxWriter()

        cell1 = Cell(10, 10, 10, 90, 90, 90)
        cell2 = Cell(10, 10, 11, 90, 90, 120)

        frame1 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell1))
        frame2 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell2))

        traj = Trajectory([frame1, frame2])

        writer.write_box_file(traj, reset_counter=True)

        captured = capsys.readouterr()

        assert captured.out == "1 10 10 10 90 90 90\n2 10 10 11 90 90 120\n"

    def test_write_vmd(self, capsys: CaptureFixture):
        writer = BoxWriter()

        cell1 = Cell(10, 10, 10, 90, 90, 90)
        cell2 = Cell(10, 10, 11, 90, 90, 90)
        frame1 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell1))
        frame2 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell2))

        traj = Trajectory([frame1, frame2])

        print("")
        writer.write_vmd(traj)

        captured = capsys.readouterr()

        assert captured.out == """
8
Box   10 10 10    90 90 90
X   -5.0 -5.0 -5.0
X   -5.0 -5.0 5.0
X   -5.0 5.0 -5.0
X   -5.0 5.0 5.0
X   5.0 -5.0 -5.0
X   5.0 -5.0 5.0
X   5.0 5.0 -5.0
X   5.0 5.0 5.0
8
Box   10 10 11    90 90 90
X   -5.0 -5.0 -5.5
X   -5.0 -5.0 5.5
X   -5.0 5.0 -5.5
X   -5.0 5.0 5.5
X   5.0 -5.0 -5.5
X   5.0 -5.0 5.5
X   5.0 5.0 -5.5
X   5.0 5.0 5.5
"""

    def test_write(self, capsys: CaptureFixture):
        writer = BoxWriter()

        cell1 = Cell(10, 10, 10, 90, 90, 90)
        cell2 = Cell(10, 10, 11, 90, 90, 90)
        frame1 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell1))
        frame2 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell2))

        traj = Trajectory([frame1, frame2])

        print("")
        writer.format = "data"
        writer.write(traj)
        print("")
        writer.format = "vmd"
        writer.write(traj)

        captured = capsys.readouterr()

        assert captured.out == """
1 10 10 10 90 90 90
2 10 10 11 90 90 90

8
Box   10 10 10    90 90 90
X   -5.0 -5.0 -5.0
X   -5.0 -5.0 5.0
X   -5.0 5.0 -5.0
X   -5.0 5.0 5.0
X   5.0 -5.0 -5.0
X   5.0 -5.0 5.0
X   5.0 5.0 -5.0
X   5.0 5.0 5.0
8
Box   10 10 11    90 90 90
X   -5.0 -5.0 -5.5
X   -5.0 -5.0 5.5
X   -5.0 5.0 -5.5
X   -5.0 5.0 5.5
X   5.0 -5.0 -5.5
X   5.0 -5.0 5.5
X   5.0 5.0 -5.5
X   5.0 5.0 5.5
"""

    def test_write_box(self, capsys: CaptureFixture):

        cell1 = Cell(10, 10, 10, 90, 90, 90)
        cell2 = Cell(10, 10, 11, 90, 90, 90)
        frame1 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell1))
        frame2 = Frame(AtomicSystem(
            atoms=self.atoms1, pos=self.pos1, cell=cell2))

        traj = Trajectory([frame1, frame2])

        print("")
        write_box(traj, format="data")
        print("")
        write_box(traj, format="vmd")

        captured = capsys.readouterr()

        assert captured.out == """
1 10 10 10 90 90 90
2 10 10 11 90 90 90

8
Box   10 10 10    90 90 90
X   -5.0 -5.0 -5.0
X   -5.0 -5.0 5.0
X   -5.0 5.0 -5.0
X   -5.0 5.0 5.0
X   5.0 -5.0 -5.0
X   5.0 -5.0 5.0
X   5.0 5.0 -5.0
X   5.0 5.0 5.0
8
Box   10 10 11    90 90 90
X   -5.0 -5.0 -5.5
X   -5.0 -5.0 5.5
X   -5.0 5.0 -5.5
X   -5.0 5.0 5.5
X   5.0 -5.0 -5.5
X   5.0 -5.0 5.5
X   5.0 5.0 -5.5
X   5.0 5.0 5.5
"""
