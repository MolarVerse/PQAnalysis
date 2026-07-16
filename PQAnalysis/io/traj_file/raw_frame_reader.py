"""
A module containing a raw fast-path reader for xyz-family trajectory
files (xyz, vel, force).

The :py:class:`RawTrajectoryReader` streams the numeric per-frame data
of a trajectory as plain numpy arrays together with the corresponding
:py:class:`~PQAnalysis.core.cell.cell.Cell` objects, without building
:py:class:`~PQAnalysis.atomic_system.atomic_system.AtomicSystem` or
:py:class:`~PQAnalysis.core.atom.atom.Atom` objects for every frame.
It is an additive fast path intended for analyses that only need the
raw coordinates/velocities per frame (e.g. MSD and VACF) and produces
bit-identical values compared to
:py:meth:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.frame_generator`,
because it reuses the exact same line parsing routine
(:py:func:`~PQAnalysis.io.traj_file.process_lines.process_lines`).
"""

import logging

from itertools import islice

from beartype.typing import Generator, List, Tuple

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Cell
from PQAnalysis.traj import TrajectoryFormat, MDEngineFormat
from PQAnalysis.types import Np2DNumberArray
from PQAnalysis.io.base import BaseReader
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import FrameReaderError, TrajectoryReaderError
from .frame_reader import XYZFrameReader

try:
    from .process_lines import process_lines  # pylint: disable=import-error
except ModuleNotFoundError:
    from ._process_lines_py import process_lines

#: The trajectory formats supported by the raw fast-path reader.
RAW_READER_TRAJ_FORMATS = (
    TrajectoryFormat.XYZ,
    TrajectoryFormat.VEL,
    TrajectoryFormat.FORCE,
)



class RawTrajectoryReader(BaseReader):

    """
    A fast-path reader that streams raw per-frame data of xyz-family
    trajectory files.

    In contrast to
    :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`,
    this reader does not construct AtomicSystem/Atom objects per frame.
    Instead, :py:meth:`raw_frame_generator` yields
    ``(values, cell)`` tuples, where ``values`` is the ``(n_atoms, 3)``
    float32 array of the frame body (positions, velocities or forces,
    depending on the trajectory format) and ``cell`` is the unit cell
    of the frame.

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

    @runtime_type_checking
    def __init__(
        self,
        filename: str | List[str],
        traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
        md_format: MDEngineFormat | str = MDEngineFormat.PQ,
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

        Raises
        ------
        TrajectoryReaderError
            If the trajectory format is not an xyz-family format.
        """
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
    def raw_frame_generator(
        self
    ) -> Generator[Tuple[Np2DNumberArray, Cell]]:
        """
        A generator that yields the raw data of the trajectory frames.

        For every frame a tuple ``(values, cell)`` is yielded, where
        ``values`` is the ``(n_atoms, 3)`` float32 array parsed from
        the frame body (positions, velocities or forces, depending on
        the trajectory format) and ``cell`` is the unit cell of the
        frame. The values and cells are bit-identical to the ones
        produced by
        :py:meth:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.frame_generator`.

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
        """

        strip_dummy_atom = self.md_format == MDEngineFormat.QMCFC
        last_cell = None

        for filename in self.filenames:
            with open(filename, "r", encoding="utf-8") as file:

                while True:
                    header_line = self._next_header_line(file)

                    if header_line is None:
                        break

                    n_atoms, cell, cell_is_vacuum = self._parse_header_line(
                        header_line
                    )

                    # comment line + n_atoms atom lines
                    body_lines = list(islice(file, n_atoms + 1))

                    if len(body_lines) != n_atoms + 1:
                        self.logger.error(
                            (
                                f"Unexpected end of file {filename}: "
                                "incomplete frame."
                            ),
                            exception=FrameReaderError,
                        )

                    values = self._parse_body_lines(body_lines, n_atoms)

                    if strip_dummy_atom:
                        values = self._strip_dummy_atom(body_lines, values)

                    if cell_is_vacuum and last_cell is not None:
                        cell = last_cell

                    last_cell = cell

                    yield values, cell

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
            self._count_frames_in_file(filename)
            for filename in self.filenames
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

    def _parse_body_lines(
        self,
        body_lines: List[str],
        n_atoms: int,
    ) -> Np2DNumberArray:
        """
        Parses the numeric values of the frame body.

        Parameters
        ----------
        body_lines : List[str]
            The lines of the frame body (comment line + atom lines).
        n_atoms : int
            The number of atoms in the frame.

        Returns
        -------
        Np2DNumberArray
            The (n_atoms, 3) float32 array of values parsed from the
            atom lines.

        Raises
        ------
        FrameReaderError
            If an atom line cannot be parsed.
        """

        try:
            return process_lines(body_lines[1:], n_atoms)
        except ValueError:
            self.logger.error(
                'Invalid file format in xyz coordinates of Frame.',
                exception=FrameReaderError,
            )

        return None  # pragma: no cover - logger.error raises

    def _strip_dummy_atom(
        self,
        body_lines: List[str],
        values: Np2DNumberArray,
    ) -> Np2DNumberArray:
        """
        Strips the leading QMCFC dummy atom row from the values.

        Parameters
        ----------
        body_lines : List[str]
            The lines of the frame body (comment line + atom lines).
        values : Np2DNumberArray
            The parsed values of the frame body.

        Returns
        -------
        Np2DNumberArray
            The values without the leading dummy atom row.

        Raises
        ------
        FrameReaderError
            If the first atom of the frame is not X.
        """

        first_atom = body_lines[1].split(None, 1)[0]

        if first_atom.upper() != 'X':
            self.logger.error(
                (
                    'The first atom in one of the frames is not X. '
                    'Please use PQ (default) md engine instead'
                ),
                exception=FrameReaderError,
            )

        return values[1:]
