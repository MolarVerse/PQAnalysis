"""
A module containing a raw fast-path reader for xyz-family trajectory
files (xyz, vel, force).

The :py:class:`RawTrajectoryReader` streams the numeric per-frame data
of a trajectory as plain numpy arrays together with the corresponding
:py:class:`~PQAnalysis.core.cell.cell.Cell` objects, without building
:py:class:`~PQAnalysis.atomic_system.atomic_system.AtomicSystem` or
:py:class:`~PQAnalysis.core.atom.atom.Atom` objects for every frame.
It is an additive fast path intended for analyses that only need the
raw coordinates or velocities per frame (e.g. MSD and VACF). The
frames are parsed from large byte chunks by the slab parser
(:py:mod:`~PQAnalysis.io.traj_file._slab_parser`). Its default
``float32`` mode is bit-identical to the line-based trajectory reader;
the optional ``float64`` mode instead preserves the precision of the
text with direct ``strtod`` conversions. When the compiled slab parser
is not available, the pure Python implementation
(:py:mod:`~PQAnalysis.io.traj_file._slab_parser_py`), which reuses
the current per-line machinery, is used instead.
"""

import logging

from itertools import islice
from os.path import getsize

import numpy as np

from beartype.typing import Generator, List, Tuple

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Cell
from PQAnalysis.traj import TrajectoryFormat, MDEngineFormat
from PQAnalysis.types import (
    Np1DNumberArray,
    Np2DNumberArray,
    NpnDNumberArray,
)
from PQAnalysis.io.base import BaseReader
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import FrameReaderError, TrajectoryReaderError
from .frame_reader import XYZFrameReader

# the status/mode constants are shared by both slab parser
# implementations and defined once in the pure Python module
from ._slab_parser_py import (
    MODE_XYZ,
    MODE_XYZ64,
    STATUS_BAD_HEADER,
    STATUS_EOF,
    STATUS_NEED_MORE,
)

try:
    from ._slab_parser import (  # pylint: disable=import-error
        parse_body,
        parse_xyz_frames,
        scan_header,
    )
except ModuleNotFoundError:
    from ._slab_parser_py import parse_body, parse_xyz_frames, scan_header

#: The trajectory formats supported by the raw fast-path reader.
RAW_READER_TRAJ_FORMATS = (
    TrajectoryFormat.XYZ,
    TrajectoryFormat.VEL,
    TrajectoryFormat.FORCE,
)

#: The chunk size (in bytes) of the buffered slab reads.
_CHUNK_SIZE = 8 * 1024 * 1024



class RawTrajectoryReader(BaseReader):

    """
    A fast-path reader that streams raw per-frame data of xyz-family
    trajectory files.

    In contrast to
    :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`,
    this reader does not construct AtomicSystem/Atom objects per frame.
    Instead, :py:meth:`raw_frame_generator` yields
    ``(values, cell)`` tuples, where ``values`` is the ``(n_atoms, 3)``
    array of the frame body (positions, velocities or forces,
    depending on the trajectory format) in the configured precision
    and ``cell`` is the unit cell of the frame.

    The reader follows the exact same semantics as
    :py:meth:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.frame_generator`:

    - Multiple files are read one after another.
    - For the QMCFC MD engine format the leading dummy atom row is
      stripped from every frame (and it is checked to be an ``X``
      atom).
    - Frames without box information in the header (vacuum frames)
      inherit the cell of the last frame that had one - also across
      file boundaries.

    As a performance optimization, the reader caches Cell objects by
    the (textual) box information of the header line. Consecutive
    frames with an identical header box string share the same Cell
    object (NPT trajectories with changing boxes still get a new Cell
    per unique box string). The yielded Cell objects must therefore be
    treated as immutable by consumers.

    For topology-dependent setup (e.g. selections),
    :py:meth:`read_first_frame` reads only the first frame of the
    trajectory the normal way and returns it as an AtomicSystem. This
    does not consume any frames of :py:meth:`raw_frame_generator`:
    every call to :py:meth:`raw_frame_generator` always streams the
    trajectory from the very first frame, so analyses can bootstrap
    their topology from :py:meth:`read_first_frame` and afterwards
    still consume every frame of the trajectory exactly once and in
    order.
    """

    # Set up the logger
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    #: The default slab parser body mode of this reader (a name token
    #: plus three float32 values per atom line).
    _SLAB_MODE = MODE_XYZ

    #: The error message used when a frame body line cannot be parsed.
    _BODY_ERROR_MESSAGE = 'Invalid file format in xyz coordinates of Frame.'

    @runtime_type_checking
    def __init__(
        self,
        filename: str | List[str],
        traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
        md_format: MDEngineFormat | str = MDEngineFormat.PQ,
        dtype: str = "float32",
    ) -> None:
        """
        Parameters
        ----------
        filename : str or list of str
            The name of the file to read from or a list of filenames
            to read from.
        traj_format : TrajectoryFormat | str, optional
            The format of the trajectory. Default is
            TrajectoryFormat.AUTO. The format is inferred from the
            file extension. Only the xyz-family formats XYZ, VEL and
            FORCE are supported by this reader.
        md_format : MDEngineFormat | str, optional
            The format of the MD engine. Default is MDEngineFormat.PQ.
        dtype : {"float32", "float64"}, optional
            Numeric precision of xyz-family frame values. ``float32``
            preserves the established trajectory-reader behavior;
            ``float64`` parses directly with ``strtod``. Default is
            ``float32``.

        Raises
        ------
        TrajectoryReaderError
            If the trajectory format is not an xyz-family format.
        ValueError
            If dtype is not ``float32`` or ``float64``.
        """
        if dtype not in {"float32", "float64"}:
            raise ValueError(
                "dtype must be either 'float32' or 'float64', "
                f"got {dtype!r}."
            )

        super().__init__(filename)

        if not self.multiple_files:
            self.filenames = [self.filename]

        self.traj_format = TrajectoryFormat((traj_format, self.filenames[0]))

        if self.traj_format not in RAW_READER_TRAJ_FORMATS:
            self.logger.error(
                (
                    "The raw trajectory reader supports only the "
                    f"{[f.value for f in RAW_READER_TRAJ_FORMATS]} "
                    f"trajectory formats, got {self.traj_format}."
                ),
                exception=TrajectoryReaderError,
            )

        self.md_format = MDEngineFormat(md_format)

        self.dtype = dtype
        self._slab_mode = (
            MODE_XYZ64 if self._SLAB_MODE == MODE_XYZ and dtype == "float64"
            else self._SLAB_MODE
        )

        # Cache of Cell objects keyed by the box substring of the
        # header line, so that unchanged boxes reuse the same Cell
        # object instead of rebuilding it for every frame.
        self._cell_cache = {}

    def read_first_frame(self) -> AtomicSystem:
        """
        Reads only the first frame of the trajectory the normal way.

        This is meant as a topology bootstrap for analyses that use
        :py:meth:`raw_frame_generator`: the first frame is read as a
        full AtomicSystem (including Atom objects), so that
        selections/topologies can be built from it. The raw frame
        stream is not affected by this method - it always starts at
        the first frame.

        Returns
        -------
        AtomicSystem
            The first frame of the trajectory.

        Raises
        ------
        TrajectoryReaderError
            If the trajectory contains no frames.
        """

        frame_reader = XYZFrameReader(md_format=self.md_format)

        for filename in self.filenames:
            with open(filename, "r", encoding="utf-8") as file:
                header_line = self._next_header_line(file)

                if header_line is None:
                    continue

                n_atoms, _, _ = self._parse_header_line(header_line)

                frame_lines = [header_line]
                frame_lines.extend(islice(file, n_atoms + 1))

                return frame_reader.read(
                    "".join(frame_lines),
                    traj_format=self.traj_format,
                )

        self.logger.error(
            "The trajectory does not contain any frames.",
            exception=TrajectoryReaderError,
        )

        return None  # pragma: no cover - logger.error raises

    @runtime_type_checking
    def raw_frame_generator(self) -> Generator[Tuple[Np2DNumberArray, Cell]]:
        """
        A generator that yields the raw data of the trajectory frames.

        For every frame a tuple ``(values, cell)`` is yielded, where
        ``values`` is the ``(n_atoms, 3)`` array parsed from the frame
        body (positions, velocities or forces, depending on the
        trajectory format) and ``cell`` is the unit cell of the frame.
        The array dtype is selected at construction. In ``float32``
        mode the values and cells are bit-identical to those produced
        by the line-based trajectory reader; ``float64`` mode retains
        additional digits present in the source text.

        The generator always starts at the first frame of the
        trajectory, so it can be restarted by simply calling this
        method again.

        Yields
        ------
        Generator[Tuple[Np2DNumberArray, Cell]]
            The raw values and the cell of the frames of the
            trajectory.

        Raises
        ------
        FrameReaderError
            If a frame of the trajectory is incomplete or its body
            cannot be parsed.
        ValueError
            If the atom count of a frame header cannot be parsed as
            an integer.
        """

        strip_dummy_atom = self.md_format == MDEngineFormat.QMCFC
        last_cell = None

        for filename in self.filenames:
            with open(filename, "rb") as file:
                buffer = b""
                offset = 0
                at_eof = False
                forced_n_atoms = -1

                while True:
                    (
                        status,
                        n_atoms,
                        box_bytes,
                        header_token,
                        body_offset,
                    ) = scan_header(buffer, offset, at_eof, forced_n_atoms)

                    if status == STATUS_EOF:
                        break

                    if status == STATUS_NEED_MORE:
                        buffer, offset, at_eof = self._refill(
                            file, buffer, offset
                        )
                        continue

                    if status == STATUS_BAD_HEADER:
                        # replicate the error order of the line based
                        # header parsing: the box substring is
                        # validated (and cached) before the atom
                        # count is converted
                        self._cell_from_box_bytes(box_bytes)

                        forced_n_atoms = int(header_token.decode("utf-8"))

                        if forced_n_atoms < 0:
                            raise ValueError(
                                "Indices for islice() must be None or "
                                "an integer: 0 <= x <= sys.maxsize."
                            )

                        continue

                    cell, cell_is_vacuum = self._cell_from_box_bytes(box_bytes)

                    try:
                        (
                            body_status,
                            values,
                            first_name,
                            next_offset,
                        ) = parse_body(
                            buffer,
                            body_offset,
                            n_atoms,
                            at_eof,
                            strip_dummy_atom,
                            self._slab_mode,
                        )
                    except EOFError:
                        self.logger.error(
                            (
                                f"Unexpected end of file {filename}: "
                                "incomplete frame."
                            ),
                            exception=FrameReaderError,
                        )
                    except ValueError:
                        self.logger.error(
                            self._BODY_ERROR_MESSAGE,
                            exception=FrameReaderError,
                        )

                    if body_status == STATUS_NEED_MORE:
                        buffer, offset, at_eof = self._refill(
                            file, buffer, offset
                        )
                        continue

                    if strip_dummy_atom:
                        values = self._strip_dummy_values(first_name, values)

                    if cell_is_vacuum and last_cell is not None:
                        cell = last_cell

                    last_cell = cell
                    offset = next_offset
                    forced_n_atoms = -1

                    yield values, cell

    # The checks stay together so every compatibility failure has one fallback.
    # pylint: disable-next=too-complex,too-many-locals,too-many-return-statements,too-many-branches,too-many-statements
    def try_read_all_frames(
        self,
        expected_n_atoms: int,
        max_bytes: int,
        expected_n_frames: int | None = None,
        include_cells: bool = True,
        *,
        include_box_lengths: bool = False,
    ) -> Tuple[NpnDNumberArray, List[Cell] | Np2DNumberArray] | None:
        """
        Read a fixed-topology xyz-family trajectory as one batch.

        The peak input-buffer and numeric-array working set must fit below
        ``max_bytes``. If it does not, or if the file requires one of the
        permissive line-reader compatibility paths, ``None`` is returned and
        the caller can use :meth:`raw_frame_generator` instead.

        Numeric tokens are converted by the same slab-parser routines as
        the streaming path. The returned values therefore have identical
        float32 or float64 bit patterns in the same frame order.

        Parameters
        ----------
        expected_n_atoms : int
            Number of physical atoms expected after QMCFC dummy removal.
        max_bytes : int
            Maximum combined size of the input byte buffers and output
            value arrays.
        expected_n_frames : int | None, optional
            Expected total frame count, used as an additional consistency
            check when the caller already knows it.
        include_cells : bool, optional
            Build the per-frame cell list. When false, unique box headers are
            still validated but the returned cell list is empty.
        include_box_lengths : bool, optional
            Return per-frame orthorhombic box lengths as a float64 array
            instead of Cell objects. Explicit triclinic boxes and a vacuum
            cell without a preceding periodic cell are not compatible with
            this mode and cause the bounded batch path to return ``None``.
            This option requires ``include_cells=False``.

        Returns
        -------
        tuple or None
            ``(values, cell_data)`` with values shaped
            ``(n_frames, expected_n_atoms, 3)``. ``cell_data`` is either the
            Cell list or, when requested, an ``(n_frames, 3)`` float64 box
            length array. Returns ``None`` when the bounded batch path cannot
            be used.
        """
        if include_cells and include_box_lengths:
            raise ValueError(
                "include_cells and include_box_lengths are mutually exclusive"
            )

        if self._slab_mode not in {MODE_XYZ, MODE_XYZ64}:
            return None

        if expected_n_atoms < 0 or max_bytes <= 0:
            return None

        input_bytes = sum(getsize(filename) for filename in self.filenames)
        if input_bytes > max_bytes:
            return None

        strip_dummy_atom = self.md_format == MDEngineFormat.QMCFC
        raw_n_atoms = expected_n_atoms + int(strip_dummy_atom)
        frame_size = raw_n_atoms + 2
        buffers = []
        frame_counts = []
        output_bytes = 0
        itemsize = 4 if self.dtype == "float32" else 8

        single_file_count = (
            expected_n_frames if expected_n_frames is not None and
            len(self.filenames) == 1 else None
        )

        for filename in self.filenames:
            with open(filename, "rb") as file:
                data = file.read()

            if data and not data.endswith(b"\n"):
                data += b"\n"
                input_bytes += 1

            if single_file_count is None:
                n_lines = data.count(b"\n")
                n_frames, remainder = divmod(n_lines, frame_size)
                if remainder != 0:
                    return None
            else:
                n_frames = single_file_count

            output_bytes += n_frames * expected_n_atoms * 3 * itemsize

            if include_box_lengths:
                output_bytes += n_frames * 3 * np.dtype(np.float64).itemsize

            peak_output_bytes = output_bytes * (
                2 if len(self.filenames) > 1 else 1
            )
            if input_bytes + peak_output_bytes > max_bytes:
                return None

            buffers.append(data)
            frame_counts.append(n_frames)

        if (
            expected_n_frames is not None and
            sum(frame_counts) != expected_n_frames
        ):
            return None

        value_batches = []
        box_length_batches = []
        box_length_cache = {}
        cells = []
        last_cell = None
        last_box_lengths = None

        try:
            for data, n_frames in zip(buffers, frame_counts):
                values, box_headers, first_names = parse_xyz_frames(
                    data,
                    n_frames,
                    raw_n_atoms,
                    strip_dummy_atom,
                    self._slab_mode,
                )

                if strip_dummy_atom:
                    for first_name in first_names:
                        self._validate_dummy_name(first_name)

                if include_cells:
                    for box_bytes in box_headers:
                        cell, cell_is_vacuum = self._cell_from_box_bytes(
                            box_bytes
                        )

                        if cell_is_vacuum and last_cell is not None:
                            cell = last_cell

                        last_cell = cell
                        cells.append(cell)
                elif include_box_lengths:
                    parsed_box_lengths = self._parse_orthorhombic_box_lengths(
                        box_headers,
                        box_length_cache,
                        last_box_lengths,
                    )

                    if parsed_box_lengths is None:
                        return None

                    box_lengths, last_box_lengths = parsed_box_lengths
                    box_length_batches.append(box_lengths)
                else:
                    for box_bytes in dict.fromkeys(box_headers):
                        self._cell_from_box_bytes(box_bytes)

                value_batches.append(values)
        except (EOFError, ValueError):
            return None

        if len(value_batches) == 1:
            all_values = value_batches[0]
        else:
            all_values = np.concatenate(value_batches, axis=0)

        if not include_box_lengths:
            return all_values, cells

        if len(box_length_batches) == 1:
            all_box_lengths = box_length_batches[0]
        else:
            all_box_lengths = np.concatenate(box_length_batches, axis=0)

        return all_values, all_box_lengths

    @staticmethod
    def _parse_orthorhombic_box_lengths(
        box_headers: List[bytes],
        cache: dict[bytes, Tuple[float, float, float]],
        inherited: Tuple[float, float, float] | None,
    ) -> Tuple[
        Np2DNumberArray,
        Tuple[float, float, float] | None,
    ] | None:
        """Parse exact box lengths without constructing Cell matrices."""
        frame_box_lengths = []
        last_box_lengths = inherited

        for box_bytes in box_headers:
            box_lengths = cache.get(box_bytes)

            if box_lengths is None:
                box_values = tuple(
                    float(value) for value in box_bytes.split()
                )

                if len(box_values) == 0:
                    if last_box_lengths is None:
                        return None

                    frame_box_lengths.append(last_box_lengths)
                    continue

                if len(box_values) not in {3, 6}:
                    return None

                if len(box_values) == 6 and box_values[3:] != (
                    90.0,
                    90.0,
                    90.0,
                ):
                    return None

                box_lengths = box_values[:3]

                if (
                    0.0 in box_lengths or
                    (box_lengths[0] * box_lengths[1]) * box_lengths[2] > 1e100
                ):
                    return None

                cache[box_bytes] = box_lengths

            last_box_lengths = box_lengths
            frame_box_lengths.append(box_lengths)

        box_lengths = np.asarray(
            frame_box_lengths,
            dtype=np.float64,
        ).reshape((-1, 3))

        return box_lengths, last_box_lengths

    def count_frames(self) -> int:
        """
        Counts the number of frames of the trajectory.

        The count is done with a cheap single-pass block scan of the
        files, without materializing the lines of the files. The
        number of atoms is taken from the first line of every file,
        exactly as in the frame counting of
        :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`.

        Returns
        -------
        int
            The total number of frames of the trajectory.

        Raises
        ------
        TrajectoryReaderError
            If the number of lines of a file is not divisible by its
            frame size or the number of atoms in the first line of a
            file is invalid.
        """

        return sum(
            self._count_frames_in_file(filename) for filename in self.filenames
        )

    def _count_frames_in_file(self, filename: str) -> int:
        """
        Counts the number of frames in a single trajectory file.

        Parameters
        ----------
        filename : str
            The name of the file to count the frames of.

        Returns
        -------
        int
            The number of frames in the file.

        Raises
        ------
        TrajectoryReaderError
            If the number of lines in the file is not divisible by
            the frame size or the number of atoms in the first line
            is invalid.
        """

        n_lines = self._count_lines(filename)

        if n_lines == 0:
            return 0

        with open(filename, "r", encoding="utf-8") as file:
            try:
                n_atoms = int(file.readline().split()[0])
            except (ValueError, IndexError):
                self.logger.error(
                    (
                        "Invalid number of atoms in the first line "
                        f"of file {filename}."
                    ),
                    exception=TrajectoryReaderError,
                )

        # +2 for the cell/atom_count + comment lines
        frame_size = n_atoms + 2

        n_frames, remainder = divmod(n_lines, frame_size)

        if remainder != 0:
            self.logger.error(
                (
                    "The number of lines in the file is not divisible "
                    f"by the number of atoms {n_atoms} "
                    "in the first line."
                ),
                exception=TrajectoryReaderError,
            )

        return n_frames

    @staticmethod
    def _count_lines(filename: str) -> int:
        """
        Counts the lines of a file with a block scan.

        A trailing line without a final newline character is counted
        as a line, matching the semantics of ``readlines()``.

        Parameters
        ----------
        filename : str
            The name of the file to count the lines of.

        Returns
        -------
        int
            The number of lines in the file.
        """

        block_size = 1 << 20
        n_lines = 0
        last_block = b"\n"

        with open(filename, "rb") as file:
            while True:
                block = file.read(block_size)

                if not block:
                    break

                n_lines += block.count(b"\n")
                last_block = block

        if not last_block.endswith(b"\n"):
            n_lines += 1

        return n_lines

    @staticmethod
    def _next_header_line(file) -> str | None:
        """
        Reads the next non-blank line of a file.

        Parameters
        ----------
        file : io.TextIOBase
            The file object to read from.

        Returns
        -------
        str | None
            The next non-blank line or None if the end of the file
            was reached.
        """

        while True:
            line = file.readline()

            if line == "":
                return None

            if line.strip() != "":
                return line

    def _parse_header_line(
        self,
        header_line: str,
    ) -> Tuple[int, Cell, bool]:
        """
        Parses the header line of a frame.

        The Cell object is cached by the box substring of the header
        line, so that frames with textually identical box information
        share the same Cell object.

        Parameters
        ----------
        header_line : str
            The header line to parse.

        Returns
        -------
        n_atoms : int
            The number of atoms in the frame.
        cell : Cell
            The cell of the frame. A vacuum cell if the header line
            contains no box information.
        cell_is_vacuum : bool
            Whether the cell of the frame is a vacuum cell.

        Raises
        ------
        FrameReaderError
            If the header line is not valid. Either it contains too
            many or too few values.
        """

        split_header = header_line.split(None, 1)
        box_text = split_header[1] if len(split_header) == 2 else ""

        cached_cell = self._cell_cache.get(box_text)

        if cached_cell is None:
            cached_cell = self._build_cell(box_text)
            self._cell_cache[box_text] = cached_cell

        return (int(split_header[0]), *cached_cell)

    def _build_cell(self, box_text: str) -> Tuple[Cell, bool]:
        """
        Builds a Cell object from the box substring of a header line.

        Parameters
        ----------
        box_text : str
            The header line substring after the atom count.

        Returns
        -------
        cell : Cell
            The cell described by the box substring. A vacuum cell if
            the substring is empty.
        cell_is_vacuum : bool
            Whether the cell is a vacuum cell.

        Raises
        ------
        FrameReaderError
            If the box substring does not contain 0, 3 or 6 values.
        """

        box_values = box_text.split()

        if len(box_values) == 0:
            cell = Cell()
        elif len(box_values) in {3, 6}:
            cell = Cell(*(float(value) for value in box_values))
        else:
            self.logger.error(
                'Invalid file format in header line of Frame.',
                exception=FrameReaderError,
            )
            raise FrameReaderError(
                'Invalid file format in header line of Frame.'
            )

        return cell, cell.is_vacuum

    @staticmethod
    def _refill(file, buffer: bytes, offset: int) -> Tuple[bytes, int, bool]:
        """
        Reads the next chunk of a file into the parse buffer.

        The already consumed part of the buffer (everything before
        ``offset``) is dropped and the next chunk is appended, so
        that a frame spanning a chunk boundary can be re-parsed from
        its start. When the end of the file is reached, the buffer
        is terminated with a newline (if it does not already end
        with one), so that the slab parsers only ever see complete
        lines.

        Parameters
        ----------
        file : io.BufferedReader
            The (binary mode) file object to read from.
        buffer : bytes
            The current parse buffer.
        offset : int
            The offset of the first unconsumed byte of the buffer.

        Returns
        -------
        buffer : bytes
            The refilled parse buffer.
        offset : int
            The new parse offset (always 0).
        at_eof : bool
            Whether the end of the file was reached.
        """

        chunk = file.read(_CHUNK_SIZE)
        buffer = buffer[offset:] + chunk
        at_eof = chunk == b""

        if at_eof and buffer != b"" and not buffer.endswith(b"\n"):
            buffer += b"\n"

        return buffer, 0, at_eof

    def _cell_from_box_bytes(self, box_bytes: bytes) -> Tuple[Cell, bool]:
        """
        Returns the (cached) Cell of a raw header box substring.

        Parameters
        ----------
        box_bytes : bytes
            The raw box substring of the header line.

        Returns
        -------
        cell : Cell
            The cell described by the box substring. A vacuum cell
            if the substring is empty.
        cell_is_vacuum : bool
            Whether the cell is a vacuum cell.

        Raises
        ------
        FrameReaderError
            If the box substring does not contain 0, 3 or 6 values.
        """

        cached_cell = self._cell_cache.get(box_bytes)

        if cached_cell is None:
            cached_cell = self._build_cell(box_bytes.decode("utf-8"))
            self._cell_cache[box_bytes] = cached_cell

        return cached_cell

    def _strip_dummy_values(
        self,
        first_name: bytes | None,
        values: Np2DNumberArray | Np1DNumberArray,
    ) -> Np2DNumberArray | Np1DNumberArray:
        """
        Strips the leading QMCFC dummy atom row from the values.

        Parameters
        ----------
        first_name : bytes | None
            The raw name token of the first atom line of the frame,
            or None if the frame has no atom lines.
        values : numpy.ndarray
            The parsed values of the frame body.

        Returns
        -------
        numpy.ndarray
            The values without the leading dummy atom row.

        Raises
        ------
        FrameReaderError
            If the first atom of the frame is not X.
        """

        self._validate_dummy_name(first_name)

        return values[1:]

    def _validate_dummy_name(self, first_name: bytes | None) -> None:
        """Validate the leading QMCFC dummy atom name."""
        if first_name is None:
            # A QMCFC frame without any atom row matches the IndexError
            # raised by the line-based dummy atom handling.
            raise IndexError('list index out of range')

        if first_name.decode("utf-8").upper() != 'X':
            self.logger.error(
                (
                    'The first atom in one of the frames is not X. '
                    'Please use PQ (default) md engine instead'
                ),
                exception=FrameReaderError,
            )
