"""
Unit tests for the Trajectory class.
"""

import sys
import pytest
import numpy as np

from ..conftest import assert_logging, assert_logging_with_exception

from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.topology import Topology
from PQAnalysis.exceptions import PQIndexError



class TestTrajectory:

    atoms1 = [Atom("H")]
    atoms2 = [Atom("C")]
    atoms3 = [Atom("O")]
    system1 = AtomicSystem(atoms=atoms1, pos=np.array([[0, 1, 2]]))
    system2 = AtomicSystem(atoms=atoms2, pos=np.array([[1, 1, 2]]))
    system3 = AtomicSystem(atoms=atoms3, pos=np.array([[2, 1, 2]]))
    frame1 = system1
    frame2 = system2
    frame3 = system3
    frames = [frame1, frame2, frame3]

    def test__init__(self):

        traj = Trajectory(self.frames)
        assert traj.frames == self.frames

        traj.frames = []
        assert traj.frames == []

        assert Trajectory().frames == []

        traj = Trajectory(self.frame1)
        assert traj.frames == [self.frame1]

        traj = Trajectory(self.system1)
        assert traj.frames == [self.system1]

    def test_check_PBC(self):
        traj = Trajectory(self.frames)
        assert traj.check_pbc() == False

        system1 = AtomicSystem(
            atoms=self.atoms1,
            pos=np.array([[0, 1, 2]]),
            cell=Cell(10, 10, 10)
        )
        system2 = AtomicSystem(
            atoms=self.atoms2,
            pos=np.array([[1, 1, 2]]),
            cell=Cell(10, 10, 10)
        )
        frame1 = system1
        frame2 = system2
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_pbc() == True

        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array([[1, 1, 2]]))
        frame1 = system1
        frame2 = system2
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_pbc() == False

        traj = Trajectory()
        assert traj.check_pbc() == False

    def test_check_vacuum(self):
        traj = Trajectory(self.frames)
        assert traj.check_vacuum() == True

        system1 = AtomicSystem(
            atoms=self.atoms1,
            pos=np.array([[0, 1, 2]]),
            cell=Cell(10, 10, 10)
        )
        system2 = AtomicSystem(
            atoms=self.atoms2,
            pos=np.array([[1, 1, 2]]),
            cell=Cell(10, 10, 10)
        )
        frame1 = system1
        frame2 = system2
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_vacuum() == False

        system2 = AtomicSystem(atoms=self.atoms2, pos=np.array([[1, 1, 2]]))
        frame1 = system1
        frame2 = system2
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert traj.check_vacuum() == False

    def test_box_volumes(self):
        traj = Trajectory(self.frames)
        assert traj.box_volumes[0] > 10**10
        assert traj.box_volumes[1] > 10**10
        assert traj.box_volumes[2] > 10**10

        system1 = AtomicSystem(
            atoms=self.atoms1,
            pos=np.array([[0, 1, 2]]),
            cell=Cell(10, 10, 10)
        )
        system2 = AtomicSystem(
            atoms=self.atoms2,
            pos=np.array([[1, 1, 2]]),
            cell=Cell(11, 11, 11)
        )
        frame1 = system1
        frame2 = system2
        frames = [frame1, frame2]

        traj = Trajectory(frames)
        assert np.allclose(traj.box_volumes, np.array([1000, 11 * 11 * 11]))

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

    def test_window(self, caplog):
        traj = Trajectory(self.frames)

        test_frames = [traj.frames for traj in traj.window(1, 2)]
        assert test_frames == [[self.frame1], [self.frame3]]

        test_frames = [traj.frames for traj in traj.window(2, 1)]
        assert test_frames == [
            [self.frame1, self.frame2], [self.frame2, self.frame3]
        ]

        test_frames = [traj.frames for traj in traj.window(2)]
        assert test_frames == [
            [self.frame1, self.frame2], [self.frame2, self.frame3]
        ]

        test_frames = [traj.frames for traj in traj.window(1)]
        assert test_frames == [[self.frame1], [self.frame2], [self.frame3]]

        test_frames = [traj.frames for traj in traj.window(2, 2)]
        assert test_frames == [[self.frame1, self.frame2]]

        assert_logging(
            caplog,
            Trajectory.__qualname__,
            "WARNING",
            "Not all frames are included in the windows. Check the window size and gap.",
            traj.window(2, 2).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size can not be less than 1 or greater than the length of the trajectory"
            ),
            function=traj.window(0).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size can not be less than 1 or greater than the length of the trajectory"
            ),
            function=traj.window(4).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window gap can not be less than 1 or greater than the length of the trajectory"
            ),
            function=traj.window(1, 0).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window gap can not be less than 1 or greater than the length of the trajectory"
            ),
            function=traj.window(1, 4).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is less than 0 or greater than the length of the trajectory"
            ),
            function=traj.window(1, 1, trajectory_start=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "stop index is less than 0 or greater than the length of the trajectory"
            ),
            function=traj.window(1, 1, trajectory_stop=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"
            ),
            function=traj.window(1, 1, trajectory_start=2,
                                 trajectory_stop=1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            Trajectory.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size is greater than the trajectory_stop - trajectory_start"
            ),
            function=traj.window(3, 1, trajectory_start=1,
                                 trajectory_stop=3).__next__,
        )

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
        assert traj1 == traj3
        assert traj1 != traj2
        assert traj2 != traj3

        assert Trajectory() == Trajectory()
        assert traj1 != Trajectory()
        assert Trajectory() != traj1

        assert Trajectory() != 1

    def test_isclose(self):
        frame1 = AtomicSystem(
            atoms=self.atoms1,
            pos=np.array([[1.0001, 1.0002, 2.0003]]),
            cell=Cell(100.1, 100.1, 100.1)
        )
        frame2 = AtomicSystem(
            atoms=self.atoms1,
            pos=np.array([[1.0002, 1.0003, 2.0004]]),
            cell=Cell(100.1001, 100.1001, 100.1001)
        )
        frame3 = AtomicSystem(
            atoms=self.atoms1,
            pos=np.array([[1.0001, 1.0002, 2.0003]]),
            cell=Cell(100.2, 100.2, 100.2)
        )

        traj1 = Trajectory([frame1, frame1])
        traj2 = Trajectory([frame2, frame2])
        traj3 = Trajectory([frame3, frame3])

        assert traj1.isclose(traj2, atol=1.0001e-4)
        assert not traj1.isclose(traj2, atol=1e-5)

        assert traj1.isclose(traj3, rtol=1e-3)
        assert not traj1.isclose(traj3, rtol=1e-4)

    def test_append(self):
        traj = Trajectory()
        print(len(traj))
        traj.append(self.frame1)
        print(len(traj))

        assert traj.frames == [self.frame1]

    def test_pop(self):
        traj = Trajectory(self.frames)
        assert traj.pop() == self.frame3
        assert traj.frames == [self.frame1, self.frame2]

        assert traj.pop(0) == self.frame1
        assert traj.frames == [self.frame2]

        with pytest.raises(IndexError) as exception:
            traj.pop(1)
        assert str(exception.value) == "pop index out of range"

    def test_property_topology(self):
        frame1 = AtomicSystem()
        frame2 = AtomicSystem()

        traj = Trajectory([frame1, frame2])
        assert traj.topology == Topology()

        frame2 = AtomicSystem(atoms=[Atom("C")])
        traj = Trajectory([frame2, frame2])
        assert traj.topology == Topology(atoms=[Atom("C")])

        traj = Trajectory()
        assert traj.topology == Topology()

    def test_property_box_lengths(self):
        frame1 = AtomicSystem()
        frame2 = AtomicSystem()

        max_float = sys.float_info.max

        traj = Trajectory([frame1, frame2])
        assert np.allclose(
            traj.box_lengths,
            np.array(
                [
                    [max_float, max_float, max_float],
                    [max_float, max_float, max_float]
                ]
            ),
        )

        frame1 = AtomicSystem(cell=Cell(10, 10, 10))
        frame2 = AtomicSystem(cell=Cell(11, 11, 11))

        traj = Trajectory([frame1, frame2])
        assert np.allclose(
            traj.box_lengths, np.array([[10, 10, 10], [11, 11, 11]])
        )

    def test_property_cells(self):
        frame1 = AtomicSystem()
        frame2 = AtomicSystem()

        traj = Trajectory([frame1, frame2])

        assert traj.cells == [Cell(), Cell()]

        frame1 = AtomicSystem(cell=Cell(10, 10, 10))
        frame2 = AtomicSystem(cell=Cell(11, 11, 11))

        traj = Trajectory([frame1, frame2])

        assert traj.cells == [Cell(10, 10, 10), Cell(11, 11, 11)]

    def test_str_repr(self):
        traj = Trajectory(self.frames)
        assert str(traj) == "Trajectory with 3 frames"
        assert repr(traj) == "Trajectory with 3 frames"
        assert str(Trajectory()) == "Trajectory with 0 frames"
        assert repr(Trajectory()) == "Trajectory with 0 frames"

    def test_com_residue_traj(self):
        traj = Trajectory(self.frames)
        assert traj.com_residue_traj == Trajectory(
            [
                self.system1.center_of_mass_residues,
                self.system2.center_of_mass_residues,
                self.system3.center_of_mass_residues
            ]
        )

        traj = Trajectory()
        assert traj.com_residue_traj == Trajectory()
