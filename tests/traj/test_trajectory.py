import pytest

from PQAnalysis.traj.frame import Frame
from PQAnalysis.traj.trajectory import Trajectory


def test__init__():
    frames = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
    traj = Trajectory(frames)
    assert traj.frames == frames

    assert Trajectory().frames == []


def test__len__():
    frames = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
    traj = Trajectory(frames)
    assert len(traj) == 3

    assert len(Trajectory()) == 0


def test__getitem__():
    frames = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
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
    frames = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
    traj = Trajectory(frames)
    assert [frame for frame in traj] == frames

    assert [frame for frame in Trajectory()] == []


def test__contains__():
    frames = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
    traj = Trajectory(frames)
    assert Frame(0, 0, 0) in traj
    assert Frame(1, 1, 1) in traj
    assert Frame(2, 2, 2) in traj
    assert Frame(3, 3, 3) not in traj

    assert Frame(0, 0, 0) not in Trajectory()
    assert Frame(1, 1, 1) not in Trajectory()
    assert Frame(2, 2, 2) not in Trajectory()
    assert Frame(3, 3, 3) not in Trajectory()


def test__add__():
    frames1 = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
    frames2 = [Frame(3, 3, 3), Frame(4, 4, 4), Frame(5, 5, 5)]
    traj1 = Trajectory(frames1)
    traj2 = Trajectory(frames2)
    assert traj1 + traj2 == Trajectory(frames1 + frames2)

    assert Trajectory() + Trajectory() == Trajectory()
    assert Trajectory() + traj1 == traj1
    assert traj1 + Trajectory() == traj1


def test__eq__():
    frames1 = [Frame(0, 0, 0), Frame(1, 1, 1), Frame(2, 2, 2)]
    frames2 = [Frame(3, 3, 3), Frame(4, 4, 4), Frame(5, 5, 5)]
    traj1 = Trajectory(frames1)
    traj2 = Trajectory(frames2)
    assert traj1 == traj1
    assert traj2 == traj2
    assert traj1 != traj2

    assert Trajectory() == Trajectory()
    assert traj1 != Trajectory()
    assert Trajectory() != traj1

    assert Trajectory() != 1
