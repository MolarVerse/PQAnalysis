import pytest
import numpy as np

from ..conftest import assert_logging, assert_logging_with_exception

from PQAnalysis.traj import Trajectory
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem


class TestTrajectoryReader:
    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):
        with pytest.raises(FileNotFoundError) as exception:
            TrajectoryReader("tmp.xyz")
        assert str(exception.value) == "File tmp.xyz not found."

        open("tmp.xyz", "w")
        reader = TrajectoryReader("tmp.xyz")
        assert reader.filename == "tmp.xyz"
        assert reader.frames == []

    # -------------------------------------------------------------------------------- #

    @pytest.mark.usefixtures("tmpdir")
    def test_read(self):

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz")

        traj = reader.read()

        cell = Cell(1.0, 1.0, 1.0)
        atoms = [Atom(atom) for atom in ["h", "o"]]

        print(traj[0].topology.atoms)
        print(traj[1].topology.atoms)
        print(traj[0].cell)
        print(traj[1].cell)

        frame1 = AtomicSystem(
            atoms=atoms, pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), cell=cell
        )
        frame2 = AtomicSystem(
            atoms=atoms, pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]), cell=cell
        )

        assert traj[0] == frame1
        # NOTE: here cell is not none because of the consecutive reading of frames
        # Cell will be taken from the previous frame
        assert traj[1] == frame2

        reader = TrajectoryReader("tmp.xyz", md_format="qmcfc")

        with pytest.raises(FrameReaderError) as exception:
            reader.read()
        assert (
            str(exception.value)
            == "The first atom in one of the frames is not X. Please use PQ (default) md engine instead"
        )

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("X 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("X 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz", md_format="qmcfc")

        traj = reader.read()

        print(traj[0].atoms[0].name.lower())
        print(traj[1].atoms[0].name)
        print(frame1.atoms[0].name)
        print(frame2.atoms[0].name)

        assert traj[0] == frame1[1]
        assert traj[1] == frame2[1]

        traj = reader.read()

        filenames = ["tmp.xyz", "tmp.xyz"]
        reader = TrajectoryReader(filenames, md_format="qmcfc")

        ref_traj = traj + traj
        traj = reader.read()

        assert traj == ref_traj

    # -------------------------------------------------------------------------------- #

    @pytest.mark.usefixtures("tmpdir")
    def test_frame_generator(self, caplog):
        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz")

        cell = Cell(1.0, 1.0, 1.0)
        atoms = [Atom(atom) for atom in ["h", "o"]]

        frame1 = AtomicSystem(
            atoms=atoms, pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), cell=cell
        )
        frame2 = AtomicSystem(
            atoms=atoms, pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]), cell=cell
        )

        test_frames = [frame for frame in reader.frame_generator()]
        assert test_frames == [frame1, frame2]

        reader = TrajectoryReader(["tmp.xyz", "tmp.xyz"])

        test_frames = [frame for frame in reader.frame_generator()]
        assert test_frames == [frame1, frame2, frame1, frame2]

        # vacuum cell
        file = open("tmp2.xyz", "w")
        print("2", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader(["tmp.xyz", "tmp2.xyz"])

        test_frames = [frame for frame in reader.frame_generator()]
        assert test_frames == [frame1, frame2, frame1, frame2]

        test_frames = [
            frame for frame in reader.frame_generator(trajectory_start=1)]
        assert test_frames == [frame2, frame1, frame2]

        test_frames = [
            frame for frame in reader.frame_generator(trajectory_stop=3)]
        assert test_frames == [frame1, frame2, frame1]

        # TODO: test topology set when const_topology and topology is None (Pylint)

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.frame_generator(trajectory_start=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "stop index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.frame_generator(trajectory_stop=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"),
            function=reader.frame_generator(
                trajectory_start=1, trajectory_stop=1
            ).__next__,
        )

    # -------------------------------------------------------------------------------- #

    @pytest.mark.usefixtures("tmpdir")
    def test_window_generator(self, caplog):
        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        filenames = ["tmp.xyz", "tmp.xyz"]
        reader = TrajectoryReader(filenames)

        cell = Cell(1.0, 1.0, 1.0)
        atoms = [Atom(atom) for atom in ["h", "o"]]

        frame1 = AtomicSystem(
            atoms=atoms, pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), cell=cell
        )
        frame2 = AtomicSystem(
            atoms=atoms, pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]), cell=cell
        )

        # TODO: Check why the following test is not working

        test_frames = list(reader.window_generator(window_size=1))
        assert test_frames == [
            Trajectory([frame1]),
            Trajectory([frame2]),
            Trajectory([frame1]),
            Trajectory([frame2]),
        ]

        test_frames = list(reader.window_generator(window_size=2))
        assert test_frames == [
            Trajectory([frame1, frame2]),
            Trajectory([frame2, frame1]),
            Trajectory([frame1, frame2]),
        ]

        test_frames = list(reader.window_generator(window_size=4))
        assert test_frames == [Trajectory([frame1, frame2, frame1, frame2])]

        test_frames = list(reader.window_generator(
            window_size=2, window_gap=2))
        assert test_frames == [
            Trajectory([frame1, frame2]),
            Trajectory([frame1, frame2]),
        ]

        test_frames = list(
            reader.window_generator(
                window_size=2, window_gap=1, trajectory_start=1, trajectory_stop=3
            )
        )
        assert test_frames == [Trajectory([frame2, frame1])]

        assert_logging(
            caplog,
            TrajectoryReader.__qualname__,
            "WARNING",
            "Not all frames are included in the windows. Check the window size and gap.",
            reader.window_generator(1, 2).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(0).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(5).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "window gap can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(3, 0).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "window gap can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(3, 5).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(
                1, 1, trajectory_start=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(
                1, 1, trajectory_start=5, trajectory_stop=6
            ).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "stop index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(
                1, 1, trajectory_stop=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "stop index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(
                1, 1, trajectory_start=0, trajectory_stop=6
            ).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"),
            function=reader.window_generator(
                1, 1, trajectory_start=1, trajectory_stop=0
            ).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"),
            function=reader.window_generator(
                1, 1, trajectory_start=1, trajectory_stop=1
            ).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=IndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size is greater than the trajectory_stop - trajectory_start"
            ),
            function=reader.window_generator(
                4, 1, trajectory_start=0, trajectory_stop=2
            ).__next__,
        )

    # -------------------------------------------------------------------------------- #
