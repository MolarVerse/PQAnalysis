import pytest

from PQAnalysis.io.boxWriter import BoxWriter
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
