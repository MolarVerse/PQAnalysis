"""
Tests of the byte-slab frame parser used by the raw fast-path
readers.

The slab parser (Cython and pure Python fallback) must produce
bitwise identical values compared to the line based
:py:meth:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.frame_generator`
path - including adversarial decimal strings near float32 rounding
boundaries, which would differ if the parser rounded twice (parsing
a float64 first and casting to float32 afterwards). The chunked
buffer handling is additionally exercised with tiny chunk sizes, so
that every frame spans multiple chunk boundaries.
"""

import importlib
import os
import sys

import numpy as np
import pytest

from PQAnalysis.analysis.vacf._raw_charge_reader import (
    RawChargeTrajectoryReader,
)
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file import RawTrajectoryReader, _slab_parser_py
from PQAnalysis.io.traj_file import raw_frame_reader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError
from PQAnalysis.traj import MDEngineFormat, TrajectoryFormat

from . import pytestmark  # pylint: disable=unused-import

try:
    from PQAnalysis.io.traj_file import _slab_parser
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _slab_parser = None

try:
    from PQAnalysis.io.traj_file import process_lines as _process_lines_ext
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _process_lines_ext = None

#: Both slab parser implementations.
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

#: Adversarial decimal strings: values near float32 rounding
#: boundaries (where single and double rounding differ), exponents,
#: leading signs, missing digits and many-digit mantissas.
ADVERSARIAL_VALUES = [
    "16777217.0",  # float32 tie at 2**24
    "16777217.000000001",  # above the tie only without double rounding
    "-16777217.000000001",
    "1.000000059604644775390625",  # midpoint of 1.0 and 1 + 2**-23
    "1.00000005960464477539062501",  # just above that midpoint
    "-1.00000005960464477539062501",
    "0.1",
    "0.2",
    "0.30000000000000004",
    "3.4028234663852886e+38",  # float32 max
    "3.4028236e38",  # overflows float32 to inf
    "-3.4028236e38",
    "1.1754943508222875e-38",  # smallest normal float32
    "1e-45",  # subnormal
    "1.4012984643248171e-45",  # smallest subnormal float32
    "7.006492321624085e-46",  # half the smallest subnormal (tie)
    "7.0064923216240854e-46",  # just above that tie
    "1e-46",  # underflows to zero
    "+4.75",
    "-.25",
    "+.5",
    "5.",
    ".5e1",
    "1E5",
    "2e-3",
    "1e+10",
    "123456789.123456789123456789123456789",
    "9.999999999999999999999999e-8",
    "0.00000000000000000000000000000000000000000000001",
    "-0.0",
    "0",
    "42",
    "6.103515624999999e-05",  # near a float32 tie with many digits
    "6.103515625e-05",
    "6.10351562500000001e-05",
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



def write_xyz(filename, frames, header="{n} 11.1 12.2 13.3"):
    """
    Writes frames (lists of atom line strings) to an xyz-family file.
    """

    with open(filename, "w", encoding="utf-8") as file:
        for frame in frames:
            print(header.format(n=len(frame)), file=file)
            print("", file=file)

            for line in frame:
                print(line, file=file)


def assert_matches_frame_generator(filenames, **kwargs):
    """
    Asserts that the raw fast-path stream matches the values and
    cells of TrajectoryReader.frame_generator bitwise.
    """

    reader = TrajectoryReader(filenames, **kwargs)
    raw_reader = RawTrajectoryReader(filenames, **kwargs)

    frames = list(reader.frame_generator())
    raw_frames = list(raw_reader.raw_frame_generator())

    assert len(raw_frames) == len(frames)

    for frame, (values, cell) in zip(frames, raw_frames):
        if reader.traj_format == TrajectoryFormat.VEL:
            expected = frame.vel
        elif reader.traj_format == TrajectoryFormat.FORCE:
            expected = frame.forces
        else:
            expected = frame.pos

        assert values.dtype == np.float32
        assert np.array_equal(values, expected)
        assert np.array_equal(cell.box_matrix, frame.cell.box_matrix)
        assert cell.is_vacuum == frame.cell.is_vacuum

    return raw_frames



@pytest.mark.usefixtures("tmpdir", "parser_module")
class TestSlabParserAdversarialFloats:

    def test_xyz_values_bitwise_identical(self):
        # cycle the adversarial strings through all three columns
        values = ADVERSARIAL_VALUES
        lines = [
            f"a{i} {values[i % len(values)]} "
            f"{values[(i + 1) % len(values)]} "
            f"{values[(i + 2) % len(values)]}"
            for i in range(len(values))
        ]

        # one file with all lines in one frame and one file with one
        # line per frame (the line based reference reader requires a
        # constant frame size per file)
        write_xyz("tmp1.xyz", [lines, lines])
        write_xyz("tmp2.xyz", [[line] for line in lines])

        assert_matches_frame_generator(["tmp1.xyz", "tmp2.xyz"])

    def test_charge_values_bitwise_identical(self):
        lines = [
            f"a{i} {value}" for i, value in enumerate(ADVERSARIAL_VALUES)
        ]

        write_xyz("tmp1.chrg", [lines, lines])
        write_xyz("tmp2.chrg", [[line] for line in lines])

        filenames = ["tmp1.chrg", "tmp2.chrg"]

        reader = TrajectoryReader(
            filenames, traj_format=TrajectoryFormat.CHARGE
        )
        raw_reader = RawChargeTrajectoryReader(filenames)

        frames = list(reader.frame_generator())
        raw_frames = list(raw_reader.raw_frame_generator())

        assert len(raw_frames) == len(frames)

        for frame, (values, _) in zip(frames, raw_frames):
            assert values.dtype == np.float64
            assert np.array_equal(values, frame.charges)

    def test_whitespace_and_extra_tokens(self):
        lines = [
            "h\t1.25\t-2.5\t3.75",
            "  o   -0.125\t 0.5  12.25   ",
            "n 1.0 2.0 3.0 extra tokens 4.0",
        ]

        if _process_lines_ext is not None:
            # sscanf and strtof both stop after the longest valid
            # float prefix of the last token (the pure Python
            # process_lines fallback rejects such tokens instead)
            lines.append("c 1.0 2.0 3.0extra")

        write_xyz("tmp.xyz", [lines])

        assert_matches_frame_generator(["tmp.xyz"])



@pytest.mark.usefixtures("tmpdir", "parser_module")
class TestSlabParserChunkBoundaries:

    @pytest.mark.parametrize("chunk_size", [1, 3, 7, 17, 64, 256])
    def test_frames_spanning_chunk_boundaries(
        self, chunk_size, monkeypatch
    ):
        frames = [
            ["h 1.2345678 -2.3456789 3.4567891", "o 4.0 5.0 6.0"],
            ["h 2.2345678 -3.3456789 4.4567891", "o 5.0 6.0 7.0"],
            ["h 3.2345678 -4.3456789 5.4567891", "o 6.0 7.0 8.0"],
        ]

        write_xyz("tmp.xyz", frames)

        # a second file with a vacuum frame (cell inheritance across
        # files must survive the chunked reading)
        write_xyz("tmp2.xyz", [["n 7.0 8.0 9.0"]], header="{n}")

        monkeypatch.setattr(raw_frame_reader, "_CHUNK_SIZE", chunk_size)

        raw_frames = assert_matches_frame_generator(["tmp.xyz", "tmp2.xyz"])

        assert len(raw_frames) == 4
        assert not raw_frames[3][1].is_vacuum

    @pytest.mark.parametrize("chunk_size", [1, 7, 64])
    def test_qmcfc_and_charges_with_chunk_boundaries(
        self, chunk_size, monkeypatch
    ):
        monkeypatch.setattr(raw_frame_reader, "_CHUNK_SIZE", chunk_size)

        write_xyz(
            "tmp.xyz",
            [
                ["X 0.0 0.0 0.0", "h 1.0 2.0 3.0"],
                ["x 0.0 0.0 0.0", "h 4.0 5.0 6.0"],
            ],
        )

        raw_frames = assert_matches_frame_generator(
            ["tmp.xyz"], md_format=MDEngineFormat.QMCFC
        )

        assert raw_frames[0][0].shape == (1, 3)

        write_xyz("tmp.chrg", [["o -0.89076318", "h 0.44538159"]] * 3)

        raw_reader = RawChargeTrajectoryReader("tmp.chrg")
        raw_frames = list(raw_reader.raw_frame_generator())

        assert len(raw_frames) == 3

        for values, _ in raw_frames:
            assert np.array_equal(
                values, np.array([-0.89076318, 0.44538159])
            )

    def test_file_without_trailing_newline(self, monkeypatch):
        monkeypatch.setattr(raw_frame_reader, "_CHUNK_SIZE", 5)

        with open("tmp.xyz", "w", encoding="utf-8") as file:
            file.write("1 11.1 12.2 13.3\n\nh 1.25 -2.5 3.75")

        raw_reader = RawTrajectoryReader("tmp.xyz")
        values, _ = next(raw_reader.raw_frame_generator())

        assert np.array_equal(
            values, np.array([[1.25, -2.5, 3.75]], dtype=np.float32)
        )



@pytest.mark.usefixtures("tmpdir", "parser_module")
class TestSlabParserErrors:

    def test_trailing_garbage_raises_value_error(self):
        write_xyz("tmp.xyz", [["h 1.0 2.0 3.0"]])

        with open("tmp.xyz", "a", encoding="utf-8") as file:
            print("banana", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")
        generator = raw_reader.raw_frame_generator()

        # the valid frame is still yielded before the error is raised
        values, _ = next(generator)
        assert np.array_equal(
            values, np.array([[1.0, 2.0, 3.0]], dtype=np.float32)
        )

        with pytest.raises(ValueError):
            next(generator)

    def test_bad_count_with_invalid_box_raises_frame_reader_error(self):
        # the box substring is validated before the atom count, so
        # the FrameReaderError of the header line wins
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("banana 1.0 2.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "header line" in str(exception.value)

    def test_bad_count_with_valid_box_raises_value_error(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("banana 1.0 2.0 3.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        with pytest.raises(ValueError):
            list(raw_reader.raw_frame_generator())

    def test_underscore_count_is_accepted_like_int(self):
        # int("1_1") == 11: the slab path must replicate the int()
        # semantics of the line based header parsing
        write_xyz(
            "tmp.xyz",
            [[f"a{i} 1.0 2.0 3.0" for i in range(11)]],
            header="1_1 11.1 12.2 13.3",
        )

        raw_reader = RawTrajectoryReader("tmp.xyz")
        values, _ = next(raw_reader.raw_frame_generator())

        assert values.shape == (11, 3)

    def test_short_body_raises_incomplete_frame(self):
        with open("tmp.xyz", "w", encoding="utf-8") as file:
            print("3 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("h 1.0 2.0 3.0", file=file)
            print("o 4.0 bad 6.0", file=file)

        raw_reader = RawTrajectoryReader("tmp.xyz")

        # the truncation wins over the (also) malformed atom line
        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "incomplete frame" in str(exception.value)

    def test_charge_bad_line_raises_scalar_values_error(self):
        write_xyz("tmp.chrg", [["o -0.5", "h 0.25 0.25"]])

        raw_reader = RawChargeTrajectoryReader("tmp.chrg")

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "scalar values" in str(exception.value)

    def test_charge_qmcfc_first_atom_not_x_raises(self):
        write_xyz("tmp.chrg", [["Xe 0.0", "h 0.25"]])

        raw_reader = RawChargeTrajectoryReader(
            "tmp.chrg", md_format=MDEngineFormat.QMCFC
        )

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "is not X" in str(exception.value)

    def test_charge_incomplete_frame_raises(self):
        with open("tmp.chrg", "w", encoding="utf-8") as file:
            print("2 11.1 12.2 13.3", file=file)
            print("", file=file)
            print("o -0.5", file=file)

        raw_reader = RawChargeTrajectoryReader("tmp.chrg")

        with pytest.raises(FrameReaderError) as exception:
            list(raw_reader.raw_frame_generator())

        assert "incomplete frame" in str(exception.value)



class TestSlabParserContract:

    def test_status_and_mode_constants_are_shared(self):
        if _slab_parser is None:
            pytest.skip("compiled slab parser not available")

        # the compiled module re-imports the constants of the pure
        # Python module, so they can never drift apart
        assert _slab_parser.STATUS_FRAME == _slab_parser_py.STATUS_FRAME
        assert _slab_parser.STATUS_EOF == _slab_parser_py.STATUS_EOF
        assert (
            _slab_parser.STATUS_NEED_MORE == _slab_parser_py.STATUS_NEED_MORE
        )
        assert (
            _slab_parser.STATUS_BAD_HEADER
            == _slab_parser_py.STATUS_BAD_HEADER
        )
        assert _slab_parser.MODE_XYZ == _slab_parser_py.MODE_XYZ
        assert _slab_parser.MODE_CHARGE == _slab_parser_py.MODE_CHARGE



@pytest.mark.skipif(
    not os.environ.get("PQANALYSIS_SLAB_BENCH_DIR"),
    reason="PQANALYSIS_SLAB_BENCH_DIR not set",
)
@pytest.mark.usefixtures("parser_module")
class TestSlabParserBenchFileParity:

    """
    Full bit-parity over all frames of the local benchmark
    trajectories (only run when PQANALYSIS_SLAB_BENCH_DIR is set).
    """

    def test_bench_xyz_and_vel_parity(self):
        bench_dir = os.environ["PQANALYSIS_SLAB_BENCH_DIR"]

        for name in ("traj.xyz", "traj.vel"):
            assert_matches_frame_generator(
                [os.path.join(bench_dir, name)]
            )

    def test_bench_charge_parity(self):
        bench_dir = os.environ["PQANALYSIS_SLAB_BENCH_DIR"]
        filename = os.path.join(bench_dir, "traj.chrg")

        reader = TrajectoryReader(
            filename, traj_format=TrajectoryFormat.CHARGE
        )
        raw_reader = RawChargeTrajectoryReader(filename)

        generator = reader.frame_generator()

        for values, _ in raw_reader.raw_frame_generator():
            assert np.array_equal(values, next(generator).charges)



class TestSlabParserPyLowLevel:

    """
    Directly exercises the pure Python slab parser helpers on hand
    crafted byte buffers to reach the buffer/EOF edge branches that
    the end-to-end reader tests do not hit.
    """

    def test_scan_header_body_line_without_trailing_newline_at_eof(self):
        # a header line that is the last line of the file and has no
        # trailing newline: the parser must treat the end of the
        # buffer as the end of the line (line_end = n_data).
        status, n_atoms, box_bytes, header_token, next_offset = (
            _slab_parser_py.scan_header(b"2 11.1 12.2 13.3", 0, True)
        )

        assert status == _slab_parser_py.STATUS_FRAME
        assert n_atoms == 2
        assert box_bytes == b"11.1 12.2 13.3"
        assert header_token is None
        # the (virtual) body offset is one past the end of the buffer
        assert next_offset == len(b"2 11.1 12.2 13.3") + 1

    def test_scan_header_skips_leading_blank_lines(self):
        # two whitespace-only lines precede the header line and must
        # be skipped before the header is parsed.
        status, n_atoms, box_bytes, _, _ = _slab_parser_py.scan_header(
            b"\n   \n2 11.1 12.2 13.3\n", 0, True
        )

        assert status == _slab_parser_py.STATUS_FRAME
        assert n_atoms == 2
        assert box_bytes == b"11.1 12.2 13.3"

    def test_parse_body_incomplete_comment_line_at_eof_raises(self):
        # the comment line of the frame has no trailing newline and
        # the buffer is at EOF: the frame is incomplete.
        with pytest.raises(EOFError) as exception:
            _slab_parser_py.parse_body(
                b"comment-without-newline",
                0,
                1,
                True,
                False,
                _slab_parser_py.MODE_XYZ,
            )

        assert "incomplete frame" in str(exception.value)



class TestSlabParserPyProcessLinesFallbackImport:

    """
    Covers the ``except ModuleNotFoundError`` fallback import of
    ``process_lines`` at the top of the pure Python slab parser
    module.
    """

    def test_process_lines_fallback_import(self):
        # remember the currently wired process_lines implementation so
        # the module can be restored to exactly this state afterwards.
        original_module = _slab_parser_py.process_lines.__module__

        target = "PQAnalysis.io.traj_file.process_lines"
        saved = sys.modules.pop(target, None)

        class _Block:

            def find_spec(self, name, path=None, target=None):  # pylint: disable=unused-argument
                if name == "PQAnalysis.io.traj_file.process_lines":
                    raise ModuleNotFoundError(name)

                return None

        blocker = _Block()
        sys.meta_path.insert(0, blocker)

        try:
            importlib.reload(_slab_parser_py)

            # with the compiled process_lines blocked, the pure Python
            # implementation is imported instead.
            assert _slab_parser_py.process_lines.__module__.endswith(
                "_process_lines_py"
            )
        finally:
            sys.meta_path.remove(blocker)

            if saved is not None:
                sys.modules[target] = saved

            importlib.reload(_slab_parser_py)

        # the module is restored to exactly the state the rest of the
        # test suite relies on.
        assert _slab_parser_py.process_lines.__module__ == original_module
