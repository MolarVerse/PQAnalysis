import pytest
import numpy as np

from PQAnalysis.traj import MDEngineFormat, Trajectory, TrajectoryFormat
from PQAnalysis.io import TrajectoryReader
import PQAnalysis.io.traj_file.frame_reader as frame_reader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError, TrajectoryReaderError
from PQAnalysis.io.traj_file.api import calculate_frames_of_trajectory_file
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.exceptions import PQIndexError, PQFileNotFoundError

from ..conftest import assert_logging, assert_logging_with_exception
from . import pytestmark



@pytest.mark.usefixtures("tmpdir")
def test_calculate_number_of_frames():
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

    assert calculate_frames_of_trajectory_file("tmp.xyz") == 2


@pytest.mark.usefixtures("tmpdir")
def test_read_extxyz_trajectory_with_auto_format():
    file = open("tmp.xyz", "w")
    print("2", file=file)
    print(
        'Lattice="1 0 0 0 2 0 0 0 3" '
        'Properties=species:S:1:pos:R:3:forces:R:3',
        file=file,
    )
    print("H 0.0 0.0 0.0 0.1 0.2 0.3", file=file)
    print("O 0.0 1.0 0.0 0.4 0.5 0.6", file=file)
    print("2", file=file)
    print(
        'lattice="1 0 0 0 2 0 0 0 3" '
        'properties=species:S:1:pos:R:3:charge:R:1',
        file=file,
    )
    print("H 1.0 0.0 0.0 -0.1", file=file)
    print("O 0.0 1.0 1.0 -0.2", file=file)
    file.close()

    reader = TrajectoryReader("tmp.xyz")
    traj = reader.read()

    assert reader.traj_format is TrajectoryFormat.EXTXYZ
    assert isinstance(reader.frame_reader, frame_reader.ExtXYZFrameReader)
    assert calculate_frames_of_trajectory_file("tmp.xyz") == 2
    assert len(traj) == 2
    assert traj[0].atoms == [Atom("h"), Atom("o")]
    assert np.allclose(traj[0].pos, [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert np.allclose(traj[0].forces, [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    assert traj[0].cell == Cell(1.0, 2.0, 3.0)
    assert np.allclose(traj[1].pos, [[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]])
    assert np.allclose(traj[1].charges, [-0.1, -0.2])
    assert traj[1].cell == Cell(1.0, 2.0, 3.0)


@pytest.mark.usefixtures("tmpdir")
def test_read_variable_size_extxyz_trajectory_with_auto_format():
    file = open("tmp.extended.xyz", "w")
    print("2", file=file)
    print("Properties=species:S:1:pos:R:3", file=file)
    print("H 0.0 0.0 0.0", file=file)
    print("O 0.0 1.0 0.0", file=file)
    print("1", file=file)
    print("Properties=species:S:1:pos:R:3", file=file)
    print("C 1.0 0.0 0.0", file=file)
    file.close()

    reader = TrajectoryReader("tmp.extended.xyz")
    traj = reader.read()

    assert reader.traj_format is TrajectoryFormat.EXTXYZ
    assert reader.constant_topology is False
    assert calculate_frames_of_trajectory_file("tmp.extended.xyz") == 2
    assert len(traj) == 2
    assert traj.frames[0].atoms == [Atom("h"), Atom("o")]
    assert traj.frames[1].atoms == [Atom("c")]
    assert np.allclose(traj.frames[1].pos, [[1.0, 0.0, 0.0]])


def test_missing_xyz_file_infers_plain_xyz_for_format_only():
    assert (
        TrajectoryFormat.infer_format_from_extension("missing.xyz")
        is TrajectoryFormat.XYZ
    )


@pytest.mark.usefixtures("tmpdir")
def test_calculate_number_of_extxyz_frames_rejects_invalid_atom_count(caplog):
    file = open("tmp.extxyz", "w")
    print("not-an-int", file=file)
    print("Properties=species:S:1:pos:R:3", file=file)
    print("H 0.0 0.0 0.0", file=file)
    file.close()

    reader = TrajectoryReader("tmp.extxyz")

    assert_logging_with_exception(
        caplog,
        TrajectoryReader.__qualname__,
        exception=TrajectoryReaderError,
        logging_level="ERROR",
        message_to_test=(
            "Invalid number of atoms in the first line of file tmp.extxyz."
        ),
        function=reader.calculate_number_of_frames_per_file,
    )


@pytest.mark.usefixtures("tmpdir")
def test_calculate_number_of_extxyz_frames_skips_blank_lines():
    file = open("tmp.extxyz", "w")
    print("", file=file)
    print("1", file=file)
    print("Properties=species:S:1:pos:R:3", file=file)
    print("H 0.0 0.0 0.0", file=file)
    file.close()

    reader = TrajectoryReader("tmp.extxyz")

    assert reader.calculate_number_of_frames_per_file() == [1]


@pytest.mark.usefixtures("tmpdir")
def test_calculate_number_of_extxyz_frames_rejects_missing_comment(caplog):
    file = open("tmp.extxyz", "w")
    print("2", file=file)
    file.close()

    reader = TrajectoryReader("tmp.extxyz")

    assert_logging_with_exception(
        caplog,
        TrajectoryReader.__qualname__,
        exception=TrajectoryReaderError,
        logging_level="ERROR",
        message_to_test=(
            "The number of lines in the file is not divisible "
            "by the number of atoms 2 in one frame."
        ),
        function=reader.calculate_number_of_frames_per_file,
    )


@pytest.mark.usefixtures("tmpdir")
def test_calculate_number_of_extxyz_frames_rejects_truncated_frame(caplog):
    file = open("tmp.extxyz", "w")
    print("2", file=file)
    print("Properties=species:S:1:pos:R:3", file=file)
    print("H 0.0 0.0 0.0", file=file)
    file.close()

    reader = TrajectoryReader("tmp.extxyz")

    assert_logging_with_exception(
        caplog,
        TrajectoryReader.__qualname__,
        exception=TrajectoryReaderError,
        logging_level="ERROR",
        message_to_test=(
            "The number of lines in the file is not divisible "
            "by the number of atoms 2 in one frame."
        ),
        function=reader.calculate_number_of_frames_per_file,
    )



class TestTrajectoryReader:

    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):
        with pytest.raises(PQFileNotFoundError) as exception:
            TrajectoryReader("tmp.xyz")
        assert str(exception.value) == "File tmp.xyz not found."

        open("tmp.xyz", "w")
        reader = TrajectoryReader("tmp.xyz")
        assert reader.filename == "tmp.xyz"
        assert reader.frames == []

        reader = TrajectoryReader("tmp.xyz", md_format="qmcfc")
        assert isinstance(reader.frame_reader, frame_reader.XYZFrameReader)
        assert reader.frame_reader.md_format is MDEngineFormat.QMCFC

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
        # traj[] can also return an AtomicSystem
        print(traj[0].cell)  # pylint: disable=no-member
        print(traj[1].cell)  # pylint: disable=no-member

        frame1 = AtomicSystem(
            atoms=atoms,
            pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            cell=cell
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]),
            cell=cell
        )

        assert traj[0] == frame1
        # NOTE: here cell is not none because of the consecutive reading of frames
        # Cell will be taken from the previous frame
        assert traj[1] == frame2

        reader = TrajectoryReader("tmp.xyz", md_format="qmcfc")

        with pytest.raises(FrameReaderError) as exception:
            reader.read()
        assert (
            str(exception.value) ==
            "The first atom in one of the frames is not X. Please use PQ (default) md engine instead"
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

        print(traj[0].atoms[0].name.lower())  # pylint: disable=no-member
        print(traj[1].atoms[0].name)  # pylint: disable=no-member
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
            atoms=atoms,
            pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            cell=cell
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]),
            cell=cell
        )

        test_frames = [frame for frame in reader.frame_generator()]
        assert test_frames == [frame1, frame2]

        reader = TrajectoryReader(["tmp.xyz", "tmp.xyz"])

        test_frames = [frame for frame in reader.frame_generator()]
        assert test_frames == [frame1, frame2, frame1, frame2]

        # vacuum cell
        file = open("tmp2.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
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
            frame for frame in reader.frame_generator(trajectory_start=1)
        ]
        assert test_frames == [frame2, frame1, frame2]

        test_frames = [
            frame for frame in reader.frame_generator(trajectory_stop=3)
        ]
        assert test_frames == [frame1, frame2, frame1]

        test_frames = [
            frame for frame in reader.frame_generator(trajectory_start=2)
        ]
        assert test_frames == [frame1, frame2]

        test_frames = [
            frame for frame in reader.frame_generator(trajectory_stop=2)
        ]
        assert test_frames == [frame1, frame2]

        # Check file change after setting the reader
        file = open("tmp2.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        file.close()

        test_frames = [frame for frame in reader.frame_generator()]
        assert test_frames == [frame1, frame2, frame1]

        # revert to the original file
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

        # TODO: test topology set when const_topology and topology is None (Pylint)

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.frame_generator(trajectory_start=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "stop index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.frame_generator(trajectory_stop=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"
            ),
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
            atoms=atoms,
            pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            cell=cell
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]),
            cell=cell
        )

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

        test_frames = list(
            reader.window_generator(window_size=2, window_gap=2)
        )
        assert test_frames == [
            Trajectory([frame1, frame2]),
            Trajectory([frame1, frame2]),
        ]

        test_frames = list(
            reader.window_generator(
                window_size=2,
                window_gap=1,
                trajectory_start=1,
                trajectory_stop=3
            )
        )
        assert test_frames == [Trajectory([frame2, frame1])]
        test_frames = list(
            reader.window_generator(
                window_size=2,
                window_gap=1,
                trajectory_start=2,
                trajectory_stop=4
            )
        )
        assert test_frames == [Trajectory([frame1, frame2])]

        filenames = ["tmp.xyz", "tmp.xyz", "tmp.xyz"]
        reader = TrajectoryReader(filenames)

        test_frames = list(
            reader.window_generator(
                window_size=2, window_gap=1, trajectory_stop=2
            )
        )
        assert test_frames == [Trajectory([frame1, frame2])]

        # Test file change after setting the reader
        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        file.close()

        filenames = ["tmp.xyz", "tmp.xyz"]
        reader = TrajectoryReader(filenames)

        test_frames = list(reader.window_generator(1))
        assert test_frames == [
            Trajectory([frame1]),
            Trajectory([frame2]),
            Trajectory([frame1]),
            Trajectory([frame1]),
            Trajectory([frame2]),
            Trajectory([frame1])
        ]

        # revert to the original file
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
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(0).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(5).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window gap can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(3, 0).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window gap can not be less than 1 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(3, 5).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(1, 1,
                                             trajectory_start=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
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
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "stop index is less than 0 or greater than the length of the trajectory"
            ),
            function=reader.window_generator(1, 1,
                                             trajectory_stop=-1).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
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
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"
            ),
            function=reader.window_generator(
                1, 1, trajectory_start=1, trajectory_stop=0
            ).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "start index is greater than or equal to the stop index"
            ),
            function=reader.window_generator(
                1, 1, trajectory_start=1, trajectory_stop=1
            ).__next__,
        )

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=PQIndexError,
            logging_level="ERROR",
            message_to_test=(
                "window size is greater than the trajectory_stop - trajectory_start"
            ),
            function=reader.window_generator(
                4, 1, trajectory_start=0, trajectory_stop=2
            ).__next__,
        )

    # -------------------------------------------------------------------------------- #

    @pytest.mark.usefixtures("tmpdir")
    def test_calculate_number_of_frames_per_file(self, caplog):
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

        file = open("tmp2.xyz", "w")
        file.close()

        filenames = ["tmp.xyz", "tmp2.xyz"]
        reader = TrajectoryReader(filenames)

        assert sum(reader.calculate_number_of_frames_per_file()) == 2

        file = open("tmp2.xyz", "w")
        print("str 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader(filenames)

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=TrajectoryReaderError,
            logging_level="ERROR",
            message_to_test=(
                "Invalid number of atoms in the first line of file tmp2.xyz."
            ),
            function=reader.calculate_number_of_frames_per_file,
        )

        file = open("tmp2.xyz", "w")
        print("", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        file.close()

        reader = TrajectoryReader(filenames)

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=TrajectoryReaderError,
            logging_level="ERROR",
            message_to_test=(
                "Invalid number of atoms in the first line of file tmp2.xyz."
            ),
            function=reader.calculate_number_of_frames_per_file,
        )

        file = open("tmp2.xyz", "w")
        print("2", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        file.close()

        reader = TrajectoryReader(filenames)

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=TrajectoryReaderError,
            logging_level="ERROR",
            message_to_test=(
                "The number of lines in the file is not divisible "
                "by the number of atoms 2 in the first line."
            ),
            function=reader.calculate_number_of_frames_per_file,
        )

    # -------------------------------------------------------------------------------- #
    @pytest.mark.usefixtures("tmpdir")
    def test_cells(self, caplog):
        file = open("tmp.xyz", "w")
        print("2", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz")

        cell = Cell()

        test_cells = reader.cells
        assert test_cells == [cell]

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        file.close()

        cell = Cell(1.0, 1.0, 1.0)

        test_cells = reader.cells
        assert test_cells == [cell]

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 1.0 1.0 0.0", file=file)
        file.close()

        cell = Cell(1.0, 1.0, 1.0)

        test_cells = reader.cells
        assert test_cells == [cell, cell]

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2 2.0 2.0 2.0", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 1.0 1.0 0.0", file=file)
        file.close()

        cell1 = Cell(1.0, 1.0, 1.0)
        cell2 = Cell(2.0, 2.0, 2.0)

        test_cells = reader.cells
        assert test_cells == [cell1, cell2]

        reader = TrajectoryReader(["tmp.xyz", "tmp.xyz"])
        test_cells = reader.cells
        assert test_cells == [cell1, cell2, cell1, cell2]

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0 90.0 90.0 90.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz")

        cell = Cell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
        test_cells = reader.cells
        assert test_cells == [cell]

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0 90.0 90.0 90.0 0.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        file.close()

        assert_logging_with_exception(
            caplog,
            TrajectoryReader.__qualname__,
            exception=TrajectoryReaderError,
            logging_level="ERROR",
            message_to_test=(
                "Invalid number of arguments for box: 8 encountered "
                "in file tmp.xyz:1 = 2 1.0 1.0 1.0 90.0 90.0 90.0 0.0"
            ),
            function=reader._cell_generator().__next__,
        )
