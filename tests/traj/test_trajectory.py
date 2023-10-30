import pytest

from PQAnalysis.traj.frame import Frame
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.pbc.cell import Cell


def test__init__():
    frames = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["C"], xyz=[[1, 1, 2]]), Frame(atoms=["O"], xyz=[[2, 1, 2]])]
    traj = Trajectory(frames)
    assert traj.frames == frames

    assert Trajectory().frames == []


def test_check_PBC():
    frames = [Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10)),
              Frame(atoms=["H"], xyz=[[0, 0, 0]])]
    traj = Trajectory(frames)
    assert traj.check_PBC() == False

    frames = [Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10)),
              Frame(atoms=["H"], xyz=[[0, 0, 0]], cell=Cell(10, 10, 10))]
    traj = Trajectory(frames)
    assert traj.check_PBC() == True


def test__len__():
    frames = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["C"], xyz=[[1, 1, 2]]), Frame(atoms=["O"], xyz=[[2, 1, 2]])]
    traj = Trajectory(frames)
    assert len(traj) == 3

    assert len(Trajectory()) == 0


def test__getitem__():
    frames = [Frame(atoms=["H"], xyz=[[0, 2, 3]]), Frame(
        atoms=["C"], xyz=[[1, 2, 3]]), Frame(atoms=["O"], xyz=[[2, 2, 3]])]
    traj = Trajectory(frames)
    assert traj[0] == frames[0]
    assert traj[1] == frames[1]
    assert traj[2] == frames[2]
    assert traj[-1] == frames[-1]

    with pytest.raises(IndexError) as exception:
        traj[3]

    with pytest.raises(IndexError) as exception:
        Trajectory()[0]
    assert str(exception.value) == "list index out of range"


def test__iter__():
    frames = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["C"], xyz=[[1, 1, 2]]), Frame(atoms=["O"], xyz=[[2, 1, 2]])]
    traj = Trajectory(frames)
    assert [frame for frame in traj] == frames

    assert [frame for frame in Trajectory()] == []


def test__contains__():
    frames = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["C"], xyz=[[1, 1, 2]]), Frame(atoms=["O"], xyz=[[2, 1, 2]])]
    traj = Trajectory(frames)
    assert Frame(atoms=["H"], xyz=[[0, 1, 2]]) in traj
    assert Frame(atoms=["C"], xyz=[[1, 1, 2]]) in traj
    assert Frame(atoms=["O"], xyz=[[2, 1, 2]]) in traj
    assert Frame(atoms=["N"], xyz=[[2, 1, 2]]) not in traj

    assert Frame(atoms=["H"], xyz=[[0, 1, 2]]) not in Trajectory()
    assert Frame(atoms=["C"], xyz=[[1, 1, 2]]) not in Trajectory()
    assert Frame(atoms=["O"], xyz=[[2, 1, 2]]) not in Trajectory()
    assert Frame(atoms=["N"], xyz=[[2, 1, 2]]) not in Trajectory()


def test__add__():
    frames1 = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["H"], xyz=[[1, 1, 2]]), Frame(atoms=["H"], xyz=[[2, 1, 2]])]
    frames2 = [Frame(atoms=["H"], xyz=[[3, 1, 2]]), Frame(
        atoms=["H"], xyz=[[4, 1, 2]]), Frame(atoms=["H"], xyz=[[5, 1, 2]])]
    traj1 = Trajectory(frames1)
    traj2 = Trajectory(frames2)
    assert traj1 + traj2 == Trajectory(frames1 + frames2)

    assert Trajectory() + Trajectory() == Trajectory()
    assert Trajectory() + traj1 == traj1
    assert traj1 + Trajectory() == traj1

    frames2 = [Frame(atoms=["O"], xyz=[[3, 1, 2]]), Frame(
        atoms=["O"], xyz=[[4, 1, 2]]), Frame(atoms=["O"], xyz=[[5, 1, 2]])]

    with pytest.raises(ValueError) as exception:
        traj1 + Trajectory(frames2)
    assert str(exception.value) == "Frames are not compatible."

    with pytest.raises(TypeError) as exception:
        traj1 + frames2
    assert str(exception.value) == "only Trajectory can be added to Trajectory."


def test__eq__():
    frames1 = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["C"], xyz=[[1, 1, 2]]), Frame(atoms=["O"], xyz=[[2, 1, 2]])]
    frames2 = [Frame(atoms=["N"], xyz=[[3, 1, 2]]), Frame(
        atoms=["He"], xyz=[[4, 1, 2]]), Frame(atoms=["Ne"], xyz=[[5, 1, 2]])]
    frames3 = [Frame(atoms=["H"], xyz=[[0, 1, 2]]), Frame(
        atoms=["C"], xyz=[[1, 1, 2]]), Frame(atoms=["O"], xyz=[[2, 1, 2]])]
    traj1 = Trajectory(frames1)
    traj2 = Trajectory(frames2)
    traj3 = Trajectory(frames3)
    assert traj1 == traj1
    assert traj1 == traj3
    assert traj2 == traj2
    assert traj1 != traj2
    assert traj2 != traj3

    assert Trajectory() == Trajectory()
    assert traj1 != Trajectory()
    assert Trajectory() != traj1

    assert Trajectory() != 1
