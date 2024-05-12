import pytest
import numpy as np

from . import pytestmark

from _pytest.capture import CaptureFixture

from PQAnalysis.io import BoxWriter, write_box, BoxFileFormat, FileWritingMode
from PQAnalysis.io.exceptions import BoxWriterError, BoxFileFormatError
from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem



class TestBoxWriter:

    def test__init__(self):
        with pytest.raises(BoxFileFormatError) as exception:
            BoxWriter(filename="tmp", output_format="r")
        assert str(exception.value) == (
            "\n"
            "'r' is not a valid BoxFileFormat.\n"
            "Possible values are: BoxFileFormat.VMD, BoxFileFormat.DATA "
            "or their case insensitive string representation: "
            "vmd, data"
        )

        writer = BoxWriter(filename="tmp", output_format="vmd")
        assert writer.file is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename == "tmp"
        assert writer.output_format == BoxFileFormat.VMD

        writer = BoxWriter(filename="tmp", output_format="data")
        assert writer.file is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename == "tmp"
        assert writer.output_format == BoxFileFormat.DATA

        writer = BoxWriter(filename="tmp")
        assert writer.file is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename == "tmp"
        assert writer.output_format == BoxFileFormat.DATA

    atoms1 = [Atom("H")]
    pos1 = np.array([[0, 1, 2]])
    cell1 = Cell(10, 10, 10)

    def test__check_PBC__(self):
        system1 = AtomicSystem(
            atoms=self.atoms1,
            pos=self.pos1,
            cell=self.cell1
        )
        system2 = AtomicSystem(atoms=self.atoms1, pos=self.pos1)

        frame1 = system1
        frame2 = system2

        traj1 = Trajectory([frame1, frame1])
        traj2 = Trajectory([frame1, frame2])

        writer = BoxWriter()

        try:
            writer.__check_pbc__(traj1)
        except:
            assert False

        with pytest.raises(BoxWriterError) as exception:
            writer.__check_pbc__(traj2)
        assert str(
            exception.value
        ) == "At least on cell of the trajectory is None. Cannot write box file."

    def test_write_box_file(self, capsys: CaptureFixture):
        writer = BoxWriter()

        cell1 = Cell(10, 10, 10, 90, 90, 90)
        cell2 = Cell(10, 10, 11, 90, 90, 120)

        frame1 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell1)
        frame2 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell2)

        traj = Trajectory([frame1, frame2])

        writer.write_box_file(traj, reset_counter=True)

        captured = capsys.readouterr()

        assert captured.out == "1 10 10 10 90 90 90\n2 10 10 11 90 90 120\n"

    def test_write_vmd(self, capsys: CaptureFixture):
        writer = BoxWriter()

        cell1 = Cell(10, 10, 10, 90, 90, 90)
        cell2 = Cell(10, 10, 11, 90, 90, 90)
        frame1 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell1)
        frame2 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell2)

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
        frame1 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell1)
        frame2 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell2)

        traj = Trajectory([frame1, frame2])

        print("")
        writer.output_format = "data"
        writer.write(traj)
        print("")
        writer.output_format = "vmd"
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
        frame1 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell1)
        frame2 = AtomicSystem(atoms=self.atoms1, pos=self.pos1, cell=cell2)

        traj = Trajectory([frame1, frame2])

        print("")
        write_box(traj, output_format="data")
        print("")
        write_box(traj, output_format="vmd")

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
