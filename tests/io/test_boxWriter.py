import pytest

from _pytest.capture import CaptureFixture

from PQAnalysis.io.boxWriter import BoxWriter, write_box
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.traj.frame import Frame
from PQAnalysis.pbc.cell import Cell


def test__init__():
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


def test__check_PBC__():
    traj1 = Trajectory([Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(
        10, 10, 10)), Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10))])

    traj2 = Trajectory([Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10)),
                        Frame(atoms=["H"], xyz=[[0, 0, 0]])])

    writer = BoxWriter()

    try:
        writer.__check_PBC__(traj1)
    except:
        assert False

    with pytest.raises(ValueError) as exception:
        writer.__check_PBC__(traj2)
    assert str(
        exception.value) == "At least on cell of the trajectory is None. Cannot write box file."


def test_write_box_file(capsys: CaptureFixture):
    writer = BoxWriter()

    traj = Trajectory([Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10, 90, 90, 90)), Frame(
        atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 11, 90, 90, 120))])

    writer.write_box_file(traj)

    captured = capsys.readouterr()

    assert captured.out == "1 10 10 10 90 90 90\n2 10 10 11 90 90 120\n"


def test_write_vmd(capsys: CaptureFixture):
    writer = BoxWriter()

    traj = Trajectory([Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10, 90, 90, 90)), Frame(
        atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 11, 90, 90, 90))])

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


def test_write(capsys: CaptureFixture):
    writer = BoxWriter()

    traj = Trajectory([Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10, 90, 90, 90)), Frame(
        atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 11, 90, 90, 90))])

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


def test_write_box(capsys: CaptureFixture):
    traj = Trajectory([Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10, 90, 90, 90)), Frame(
        atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 11, 90, 90, 90))])

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
