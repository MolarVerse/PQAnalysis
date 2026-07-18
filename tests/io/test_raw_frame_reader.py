import importlib
import logging
import sys

import numpy as np
import pytest

from PQAnalysis.core import Cell
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file import RawTrajectoryReader, _slab_parser_py
from PQAnalysis.io.traj_file import raw_frame_reader
from PQAnalysis.io.traj_file.exceptions import (
    FrameReaderError,
    TrajectoryReaderError,
)
from PQAnalysis.traj import MDEngineFormat, TrajectoryFormat

from . import pytestmark

try:
    from PQAnalysis.io.traj_file import _slab_parser
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _slab_parser = None

#: Both slab parser implementations (the reader is exercised with
#: each one via the parser_module fixture).
PARSER_MODULES = [
    pytest.param(_slab_parser_py, id="python-fallback"),
    pytest.param(
        _slab_parser,
        id="cython",
        marks=pytest.mark.skipif(
            _slab_parser is None,
            reason="compiled slab parser not available",
        ),
    ),
]



@pytest.fixture(params=PARSER_MODULES)
def parser_module(request, monkeypatch):
    """
    Runs the test with each slab parser implementation wired into
    the raw frame reader module.
    """

    module = request.param

    monkeypatch.setattr(raw_frame_reader, "scan_header", module.scan_header)
    monkeypatch.setattr(raw_frame_reader, "parse_body", module.parse_body)

    return module



def assert_raw_stream_matches_frame_generator(
    filenames,
    traj_format=TrajectoryFormat.AUTO,
    md_format=MDEngineFormat.PQ,
):
    """
    Asserts that the raw fast-path stream matches the values and
    cells of TrajectoryReader.frame_generator exactly.
    """

    reader = TrajectoryReader(
        filenames, traj_format=traj_format, md_format=md_format
    )
    raw_reader = RawTrajectoryReader(
        filenames, traj_format=traj_format, md_format=md_format
    )

    frames = list(reader.frame_generator())
    raw_frames = list(raw_reader.raw_frame_generator())

    assert len(raw_frames) == len(frames)

    for frame, (values, cell) in zip(frames, raw_frames):
        if reader.traj_format == TrajectoryFormat.XYZ:
            expected_values = frame.pos
        elif reader.traj_format == TrajectoryFormat.VEL:
            expected_values = frame.vel
        else:
            expected_values = frame.forces

        assert values.dtype == np.float32
        assert values.shape == expected_values.shape
        assert np.array_equal(values, expected_values)

        assert np.array_equal(cell.box_lengths, frame.cell.box_lengths)
        assert np.array_equal(cell.box_angles, frame.cell.box_angles)
        assert np.array_equal(cell.box_matrix, frame.cell.box_matrix)
        assert cell.is_vacuum == frame.cell.is_vacuum

    return raw_frames



@pytest.mark.usefixtures("tmpdir", "parser_module")
class TestRawTrajectoryReaderEquivalence:

    def test_pq_xyz(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("o -0.1234567 0.98765432 12.3456789", file=file)
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)
            print("o -1.1234567 1.98765432 13.3456789", file=file)

        raw_frames = assert_raw_stream_matches_frame_generator(["tmp.xyz"])

        # unchanged header box text -> the same Cell object is reused
        assert raw_frames[0][1] is raw_frames[1][1]

    def test_pq_vel(self):
        with open("tmp.vel", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 0.0012345 -0.0023456 0.0034567", file=file)
            print("o -0.0001234 0.0009876 0.0123456", file=file)
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 0.0022345 -0.0033456 0.0044567", file=file)
            print("o -0.0011234 0.0019876 0.0133456", file=file)

        raw_reader = RawTrajectoryReader("tmp.vel")
        assert raw_reader.traj_format == TrajectoryFormat.VEL

        assert_raw_stream_matches_frame_generator(["tmp.vel"])

    def test_qmcfc_xyz(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("3 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("X 0.0 0.0 0.0", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("o -0.1234567 0.98765432 12.3456789", file=file)
            print("3 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("x 0.0 0.0 0.0", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)
            print("o -1.1234567 1.98765432 13.3456789", file=file)

        raw_frames = assert_raw_stream_matches_frame_generator(
            ["tmp.xyz"], md_format=MDEngineFormat.QMCFC
        )

        # the dummy atom row must be stripped
        assert raw_frames[0][0].shape == (2, 3)

    def test_qmcfc_vel(self):
        with open("tmp.vel", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("X 0.0 0.0 0.0", file=file)
            print("h 0.0012345 -0.0023456 0.0034567", file=file)

        assert_raw_stream_matches_frame_generator(
            ["tmp.vel"], md_format=MDEngineFormat.QMCFC
        )

    def test_qmcfc_first_atom_not_x_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("Xe 0.0 0.0 0.0", file=file)
            print("h 1.0 2.0 3.0", file=file)

        raw_reader = RawTrajectoryReader(
            "tmp.xyz", md_format=MDEngineFormat.QMCFC
        )

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "is not X" in str(exception.value)

    def test_multiple_files(self):
        with open("tmp1.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("o -0.1234567 0.98765432 12.3456789", file=file)

        # second file with a different number of atoms and a vacuum
        # header in its first frame (must inherit the cell of the
        # last frame of the first file)
        with open("tmp2.xyz", "w", encoding="utf-8") as file:
            print("3", file=file)
            print("", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)
            print("o -1.1234567 1.98765432 13.3456789", file=file)
            print("n 5.5555555 6.6666666 7.7777777", file=file)
            print("3 14.4 15.5 16.6", file=file)
            print("", file=file)
            print("h 3.2345678 -4.3456789 5.4567891", file=file)
            print("o -2.1234567 2.98765432 14.3456789", file=file)
            print("n 6.5555555 7.6666666 8.7777777", file=file)

        raw_frames = assert_raw_stream_matches_frame_generator(
            ["tmp1.xyz", "tmp2.xyz"]
        )

        assert len(raw_frames) == 3
        # cross-file vacuum propagation reuses the last cell object
        assert raw_frames[1][1] is raw_frames[0][1]

    def test_npt_box_change(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("1 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("1 11.2 12.3 13.4", file=file)
            print("", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)
            print("1 11.3 12.4 13.5 88.8 89.9 90.1", file=file)
            print("", file=file)
            print("h 3.2345678 -4.3456789 5.4567891", file=file)

        raw_frames = assert_raw_stream_matches_frame_generator(["tmp.xyz"])

        cells = [cell for _, cell in raw_frames]
        assert cells[0] is not cells[1]
        assert cells[1] is not cells[2]
        assert not np.array_equal(cells[0].box_lengths, cells[1].box_lengths)
        assert not np.array_equal(cells[1].box_angles, cells[2].box_angles)

    def test_vacuum_frame_inherits_last_cell(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("1 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("1", file=file)
            print("", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)
            print("1 14.4 15.5 16.6", file=file)
            print("", file=file)
            print("h 3.2345678 -4.3456789 5.4567891", file=file)

        raw_frames = assert_raw_stream_matches_frame_generator(["tmp.xyz"])

        assert raw_frames[1][1] is raw_frames[0][1]
        assert not raw_frames[2][1].is_vacuum

    def test_leading_vacuum_frames(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("1", file=file)
            print("", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("1", file=file)
            print("", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)

        raw_frames = assert_raw_stream_matches_frame_generator(["tmp.xyz"])

        assert raw_frames[0][1].is_vacuum
        assert np.array_equal(raw_frames[0][1].box_matrix, Cell().box_matrix)

    def test_exact_values(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.25 -2.5 3.75", file=file)
            print("o -0.125 0.5 12.25", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")
        values, cell = next(raw_reader.raw_frame_generator())

        assert values.dtype == np.float32
        assert np.array_equal(
            values,
            np.array(
                [[1.25, -2.5, 3.75], [-0.125, 0.5, 12.25]],
                dtype=np.float32,
            ),
        )
        assert np.array_equal(cell.box_lengths, [11.1, 12.2, 13.3])



@pytest.mark.usefixtures("tmpdir")
class TestRawTrajectoryReaderFirstFrame:

    def test_read_first_frame(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("o -0.1234567 0.98765432 12.3456789", file=file)
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 2.2345678 -3.3456789 4.4567891", file=file)
            print("o -1.1234567 1.98765432 13.3456789", file=file)

        reader = TrajectoryReader("tmp.xyz")
        expected_first_frame = next(reader.frame_generator())

        raw_reader = RawTrajectoryReader("tmp.xyz")
        first_frame = raw_reader.read_first_frame()

        assert first_frame.atoms == expected_first_frame.atoms
        assert np.array_equal(first_frame.pos, expected_first_frame.pos)
        assert first_frame.cell == expected_first_frame.cell

        # reading the first frame does not consume any raw frames
        raw_frames = list(raw_reader.raw_frame_generator())
        assert len(raw_frames) == 2
        assert np.array_equal(raw_frames[0][0], expected_first_frame.pos)

    def test_read_first_frame_qmcfc(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("3 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("X 0.0 0.0 0.0", file=file)
            print("h 1.2345678 -2.3456789 3.4567891", file=file)
            print("o -0.1234567 0.98765432 12.3456789", file=file)

        reader = TrajectoryReader("tmp.xyz", md_format=MDEngineFormat.QMCFC)
        expected_first_frame = next(reader.frame_generator())

        raw_reader = RawTrajectoryReader(
            "tmp.xyz", md_format=MDEngineFormat.QMCFC
        )
        first_frame = raw_reader.read_first_frame()

        assert first_frame.n_atoms == 2
        assert first_frame.atoms == expected_first_frame.atoms
        assert np.array_equal(first_frame.pos, expected_first_frame.pos)

    def test_read_first_frame_vel(self):
        with open("tmp.vel", "w", encoding="utf-8") as file:
            print("1 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 0.0012345 -0.0023456 0.0034567", file=file)

        raw_reader = RawTrajectoryReader("tmp.vel")
        first_frame = raw_reader.read_first_frame()

        assert np.array_equal(
            first_frame.vel,
            np.array(
                [[0.0012345, -0.0023456, 0.0034567]],
                dtype=np.float32,
            ),
        )

    def test_read_first_frame_empty_trajectory_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8"):
            pass

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(TrajectoryReaderError) as exception:
            raw_reader.read_first_frame()

        assert "does not contain any frames" in str(exception.value)



@pytest.mark.usefixtures("tmpdir")
class TestRawTrajectoryReaderCountFrames:

    def test_count_frames_matches_trajectory_reader(self):
        with open("tmp1.xyz", "w", encoding="utf-8") as file:
            for _ in range(3):
                print("2 11.1 12.2 13.3", file=file)
                print("", file=file)
                print("h 1.0 2.0 3.0", file=file)
                print("o 4.0 5.0 6.0", file=file)

        with open("tmp2.xyz", "w", encoding="utf-8") as file:
            for _ in range(2):
                print("1 11.1 12.2 13.3", file=file)
                print("", file=file)
                print("h 1.0 2.0 3.0", file=file)

        filenames = ["tmp1.xyz", "tmp2.xyz"]

        reader = TrajectoryReader(filenames)
        raw_reader = RawTrajectoryReader(filenames)

        assert raw_reader.count_frames() == sum(
            reader.calculate_number_of_frames_per_file()
        )
        assert raw_reader.count_frames() == 5

    def test_count_frames_without_trailing_newline(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            file.write("1 11.1 12.2 13.3\n")
            file.write("\n")
            file.write("h 1.0 2.0 3.0\n")
            file.write("1 11.1 12.2 13.3\n")
            file.write("\n")
            file.write("h 4.0 5.0 6.0")  # no trailing newline

        raw_reader = RawTrajectoryReader("tmp.xyz")
        assert raw_reader.count_frames() == 2

    def test_count_frames_empty_file(self):
        with open("tmp1.xyz", "w", encoding="utf-8") as file:
            print("1 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)

        with open("tmp2.xyz", "w", encoding="utf-8"):
            pass

        raw_reader = RawTrajectoryReader(["tmp1.xyz", "tmp2.xyz"])
        assert raw_reader.count_frames() == 1

    def test_count_frames_not_divisible_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)
            print("o 4.0 5.0 6.0", file=file)
            print("2 11.1 12.2 13.3", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(TrajectoryReaderError) as exception:
            raw_reader.count_frames()

        assert "not divisible" in str(exception.value)

    def test_count_frames_invalid_first_line_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("invalid", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(TrajectoryReaderError) as exception:
            raw_reader.count_frames()

        assert "Invalid number of atoms" in str(exception.value)



@pytest.mark.usefixtures("tmpdir", "parser_module")
class TestRawTrajectoryReaderErrors:

    def test_unsupported_traj_format_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("1 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)

        with pytest.raises(TrajectoryReaderError) as exception:
            RawTrajectoryReader("tmp.xyz", traj_format=TrajectoryFormat.CHARGE)

        assert "supports only" in str(exception.value)

    def test_invalid_header_line_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("1 11.1 12.2", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "header line" in str(exception.value)

    def test_incomplete_frame_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "incomplete frame" in str(exception.value)

    def test_invalid_body_line_raises(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)
            print("o 4.0 bad 6.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "Invalid file format" in str(exception.value)

    def test_negative_atom_count_raises_value_error(self):
        # a header count token that is a valid but negative integer
        # literal (with a valid box, so the box validation passes):
        # scan_header reports a bad header, the reader converts the
        # token to an int and rejects the negative count, replicating
        # the islice() error of the line based reader.
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("-3 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(ValueError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "islice()" in str(exception.value)

    def test_invalid_box_count_raises_frame_reader_error_directly(self):
        # when ERROR logging is disabled, self.logger.error does not
        # raise, so the explicit raise after it is reached (the header
        # box has only two values, which is not 0, 3 or 6).
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("1 11.1 12.2", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        original_level = raw_reader.logger.level
        raw_reader.logger.setLevel(logging.CRITICAL + 1)

        try:
            with pytest.raises(FrameReaderError) as exception:
                list(raw_reader.raw_frame_generator())
        finally:
            raw_reader.logger.setLevel(original_level)

        assert "header line" in str(exception.value)

    def test_qmcfc_frame_without_atom_rows_raises_index_error(self):
        # a QMCFC frame whose header declares zero atoms: the frame
        # has no atom row, so no first name is available and the dummy
        # atom stripping fails with the same IndexError as the line
        # based reader.
        with open("tmp.vel", "w", encoding="utf-8") as file:
            print("0 11.1 12.2 13.3", file=file)
            print("", file=file)

        raw_reader = RawTrajectoryReader(
            "tmp.vel", md_format=MDEngineFormat.QMCFC
        )

        with pytest.raises(IndexError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "list index out of range" in str(exception.value)



class TestRawFrameReaderSlabParserFallbackImport:

    """
    Covers the ``except ModuleNotFoundError`` fallback import of
    ``scan_header``/``parse_body`` from the pure Python slab parser
    module at the top of the raw frame reader module.
    """

    def test_slab_parser_fallback_import(self):
        # remember the currently wired parser functions so the module
        # can be restored to exactly this state afterwards.
        original_module = raw_frame_reader.scan_header.__module__

        target = "PQAnalysis.io.traj_file._slab_parser"
        saved = sys.modules.pop(target, None)

        class _Block:

            def find_spec(self, name, path=None, target=None):  # pylint: disable=unused-argument
                if name == "PQAnalysis.io.traj_file._slab_parser":
                    raise ModuleNotFoundError(name)

                return None

        blocker = _Block()
        sys.meta_path.insert(0, blocker)

        try:
            importlib.reload(raw_frame_reader)

            # with the compiled slab parser blocked, the pure Python
            # implementation is imported instead.
            assert raw_frame_reader.scan_header.__module__.endswith(
                "_slab_parser_py"
            )
            assert raw_frame_reader.parse_body.__module__.endswith(
                "_slab_parser_py"
            )
        finally:
            sys.meta_path.remove(blocker)

            if saved is not None:
                sys.modules[target] = saved

            importlib.reload(raw_frame_reader)

        # the module is restored to exactly the state the rest of the
        # test suite relies on.
        assert raw_frame_reader.scan_header.__module__ == original_module
