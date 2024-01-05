import pytest
import numpy as np

from PQAnalysis.traj import Frame, Trajectory, TrajectoryError
from PQAnalysis.core import Cell, Atom, AtomicSystem
from PQAnalysis.topology import Topology


class TestTrajectory:
    atoms1 = [Atom("H")]
    atoms2 = [Atom("C")]
    atoms3 = [Atom("O")]
    system1 = AtomicSystem(atoms=atoms1, pos=np.array([[0, 1, 2]]))
    system2 = AtomicSystem(atoms=atoms2, pos=np.array([[1, 1, 2]]))
    system3 = AtomicSystem(atoms=atoms3, pos=np.array([[2, 1, 2]]))
    frame1 = Frame(system1)
    frame2 = Frame(system2)
    frame3 = Frame(system3)
    frames = [frame1, frame2, frame3]

    def test__init__(self):

        traj = Trajectory(self.frames)
        assert traj.frames == self.frames

        traj.frames = []
        assert traj.frames == []

        assert Trajectory().frames == []

    def test_check_PBC(self):
        traj = Trajectory(self.frames)
        assert traj.check_PBC() == False

        system1 = AtomicSystem(atoms=self.atoms1, pos=np.array(
            [[0, 1, 2]]), cell=Cell(10, 10, 10))
        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array(
            [[1, 1, 2]]), cell=Cell(10, 10, 10))
        frame1 = Frame(system1)
        frame2 = Frame(system2)
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_PBC() == True

        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array(
            [[1, 1, 2]]))
        frame1 = Frame(system1)
        frame2 = Frame(system2)
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_PBC() == False

        traj = Trajectory()
        assert traj.check_PBC() == False

    def test_check_vacuum(self):
        traj = Trajectory(self.frames)
        assert traj.check_vacuum() == True

        system1 = AtomicSystem(atoms=self.atoms1, pos=np.array(
            [[0, 1, 2]]), cell=Cell(10, 10, 10))
        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array(
            [[1, 1, 2]]), cell=Cell(10, 10, 10))
        frame1 = Frame(system1)
        frame2 = Frame(system2)
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_vacuum() == False

        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array(
            [[1, 1, 2]]))
        frame1 = Frame(system1)
        frame2 = Frame(system2)
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_vacuum() == False

    def test_box_volumes(self):
        traj = Trajectory(self.frames)
        assert traj.box_volumes[0] > 10**10
        assert traj.box_volumes[1] > 10**10
        assert traj.box_volumes[2] > 10**10

        system1 = AtomicSystem(atoms=self.atoms1, pos=np.array(
            [[0, 1, 2]]), cell=Cell(10, 10, 10))
        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array(
            [[1, 1, 2]]), cell=Cell(11, 11, 11))
        frame1 = Frame(system1)
        frame2 = Frame(system2)
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert np.allclose(traj.box_volumes, np.array([1000, 11*11*11]))

    def test__len__(self):
        traj = Trajectory(self.frames)
        assert len(traj) == 3

        assert len(Trajectory()) == 0

    def test__getitem__(self):
        traj = Trajectory(self.frames)
        assert traj[0] == self.frames[0]
        assert traj[1] == self.frames[1]
        assert traj[2] == self.frames[2]
        assert traj[-1] == self.frames[-1]

        assert traj[0:2] == Trajectory(self.frames[0:2])

        with pytest.raises(IndexError) as exception:
            traj[3]

        with pytest.raises(IndexError) as exception:
            Trajectory()[0]
        assert str(exception.value) == "list index out of range"

    def test__iter__(self):
        traj = Trajectory(self.frames)
        assert [frame for frame in traj] == self.frames

        assert [frame for frame in Trajectory()] == []

    def test__contains__(self):
        traj = Trajectory(self.frames[0:2])
        assert self.frame1 in traj
        assert self.frame2 in traj
        assert self.frame3 not in traj

    def test__add__(self):
        traj1 = Trajectory(self.frames)
        traj2 = Trajectory(self.frames)
        assert traj1 + traj2 == Trajectory(self.frames + self.frames)

        assert Trajectory() + Trajectory() == Trajectory()
        assert Trajectory() + traj1 == traj1
        assert traj1 + Trajectory() == traj1

    def test__eq__(self):
        frames1 = self.frames[0:2]
        frames2 = self.frames[1:3]
        frames3 = self.frames[0:2]
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

    def test_append(self):
        traj = Trajectory()
        print(len(traj))
        traj.append(self.frame1)
        print(len(traj))

        assert traj.frames == [self.frame1]

    def test_property_topology(self):
        frame1 = Frame()
        frame2 = Frame()

        traj = Trajectory([frame1, frame2])
        assert traj.topology == Topology()

        frame2 = Frame(system=AtomicSystem(atoms=[Atom("C")]))
        traj = Trajectory([frame1, frame2])
        with pytest.raises(TrajectoryError) as exception:
            traj.topology
        assert str(
            exception.value) == "All frames in the trajectory must have the same topology."

        traj = Trajectory([frame2, frame2])
        assert traj.topology == Topology(atoms=[Atom("C")])

        traj = Trajectory()
        assert traj.topology == Topology()
