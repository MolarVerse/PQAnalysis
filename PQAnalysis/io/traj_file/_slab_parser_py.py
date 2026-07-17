"""
Pure Python fallback of the byte-slab frame parser.

This module implements the exact same function contract as the Cython
module :py:mod:`PQAnalysis.io.traj_file._slab_parser` on top of the
current per-line machinery
(:py:func:`~PQAnalysis.io.traj_file.process_lines.process_lines` for
the xyz-family vector lines and Python's ``float`` for the scalar
charge lines). It is used by
:py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
when the compiled extension is not available.

Both implementations operate on byte buffers (chunks of a trajectory
file) instead of decoded text lines. A frame is parsed from a given
offset in two steps:

1. :py:func:`scan_header` skips blank lines, locates the header line
   of the next frame and extracts the atom count and the raw box
   substring of the header.
2. :py:func:`parse_body` parses the frame body (comment line plus
   atom lines) into a numpy array.

Both functions never consume data on failure: when a frame is not
fully contained in the buffer, they return
:py:data:`STATUS_NEED_MORE` and the caller re-invokes them with a
refilled buffer. The caller guarantees that the buffer ends with a
newline once the end of the file is reached (``at_eof`` is True).
"""

import re

import numpy as np

from beartype.typing import Tuple

from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray

try:
    from .process_lines import process_lines  # pylint: disable=import-error
except ModuleNotFoundError:
    from ._process_lines_py import process_lines

#: Status code: a frame (header or body) was parsed successfully.
STATUS_FRAME = 0
#: Status code: the end of the trajectory file was reached cleanly.
STATUS_EOF = 1
#: Status code: the buffer ends within the frame - more data needed.
STATUS_NEED_MORE = 2
#: Status code: the atom count of the header line is not a valid
#: non-negative integer literal.
STATUS_BAD_HEADER = 3

#: Body mode: xyz-family lines with a name token and three float32
#: values per line.
MODE_XYZ = 0
#: Body mode: charge lines with a name token and exactly one float64
#: value per line.
MODE_CHARGE = 1

#: The integer literals accepted by the fast atom count parsing
#: (mirroring C ``strtol`` with base 10 on a fully consumed token).
_INT_TOKEN_RE = re.compile(rb"[+-]?[0-9]+\Z")



def _next_line_end(data: bytes, pos: int) -> int:
    """
    Finds the end of the line starting at ``pos``.

    Parameters
    ----------
    data : bytes
        The buffer to scan.
    pos : int
        The offset of the line start.

    Returns
    -------
    int
        The offset of the newline character terminating the line, or
        ``-1`` if the buffer ends before the next newline.
    """

    return data.find(b"\n", pos)



def scan_header(
    data: bytes,
    offset: int,
    at_eof: bool,
    forced_n_atoms: int = -1,
) -> Tuple[int, int, bytes | None, bytes | None, int]:
    """
    Scans for the header line of the next frame.

    Blank (whitespace-only) lines before the header line are skipped,
    exactly as the line based reader skips blank lines between
    frames.

    Parameters
    ----------
    data : bytes
        The buffer to scan.
    offset : int
        The offset to start scanning at.
    at_eof : bool
        Whether the buffer contains the complete rest of the file.
    forced_n_atoms : int, optional
        The atom count to use when the count token of the header line
        is not a plain integer literal but was successfully parsed by
        the caller (e.g. an integer literal with underscores).
        Default is -1 (disabled).

    Returns
    -------
    status : int
        One of :py:data:`STATUS_FRAME`, :py:data:`STATUS_EOF`,
        :py:data:`STATUS_NEED_MORE` and
        :py:data:`STATUS_BAD_HEADER`.
    n_atoms : int
        The atom count of the frame (-1 unless the status is
        :py:data:`STATUS_FRAME`).
    box_bytes : bytes | None
        The raw box substring of the header line (the part after the
        atom count token, without leading whitespace). None unless
        the status is :py:data:`STATUS_FRAME` or
        :py:data:`STATUS_BAD_HEADER`.
    header_token : bytes | None
        The raw atom count token of the header line. Only set for
        :py:data:`STATUS_BAD_HEADER`.
    next_offset : int
        The offset of the frame body (the comment line) for
        :py:data:`STATUS_FRAME` and :py:data:`STATUS_BAD_HEADER`,
        otherwise the (unconsumed) input offset.
    """

    n_data = len(data)
    pos = offset

    while True:
        if pos >= n_data:
            if at_eof:
                return (STATUS_EOF, -1, None, None, pos)

            return (STATUS_NEED_MORE, -1, None, None, offset)

        line_end = _next_line_end(data, pos)

        if line_end < 0:
            if not at_eof:
                return (STATUS_NEED_MORE, -1, None, None, offset)

            line_end = n_data

        line = data[pos:line_end]

        if line.strip() != b"":
            break

        pos = line_end + 1

    fields = line.split(None, 1)
    token = fields[0]
    box_bytes = fields[1] if len(fields) == 2 else b""
    body_offset = line_end + 1

    if _INT_TOKEN_RE.match(token):
        n_atoms = int(token)

        if n_atoms >= 0:
            return (STATUS_FRAME, n_atoms, box_bytes, None, body_offset)

    if forced_n_atoms >= 0:
        return (STATUS_FRAME, forced_n_atoms, box_bytes, None, body_offset)

    return (STATUS_BAD_HEADER, -1, box_bytes, token, body_offset)



def parse_body(
    data: bytes,
    offset: int,
    n_atoms: int,
    at_eof: bool,
    want_first_name: bool,
    mode: int,
) -> Tuple[
    int,
    Np2DNumberArray | Np1DNumberArray | None,
    bytes | None,
    int,
]:
    """
    Parses the body (comment line plus atom lines) of a frame.

    The completeness of the frame is checked before any line is
    parsed, so that a truncated frame at the end of the file is
    always reported as incomplete - even if it also contains
    malformed lines.

    Parameters
    ----------
    data : bytes
        The buffer to parse from.
    offset : int
        The offset of the comment line of the frame.
    n_atoms : int
        The number of atom lines of the frame.
    at_eof : bool
        Whether the buffer contains the complete rest of the file.
    want_first_name : bool
        Whether to extract the name token of the first atom line
        (needed for the QMCFC dummy atom check).
    mode : int
        The body mode, either :py:data:`MODE_XYZ` or
        :py:data:`MODE_CHARGE`.

    Returns
    -------
    status : int
        :py:data:`STATUS_FRAME` or :py:data:`STATUS_NEED_MORE`.
    values : numpy.ndarray | None
        The parsed values of the frame body: a ``(n_atoms, 3)``
        float32 array for :py:data:`MODE_XYZ` or a ``(n_atoms,)``
        float64 array for :py:data:`MODE_CHARGE`. None unless the
        status is :py:data:`STATUS_FRAME`.
    first_name : bytes | None
        The name token of the first atom line if requested and the
        frame has at least one atom line, otherwise None.
    next_offset : int
        The offset directly after the frame for
        :py:data:`STATUS_FRAME`, otherwise the (unconsumed) input
        offset.

    Raises
    ------
    EOFError
        If the buffer contains the complete rest of the file and the
        frame is incomplete.
    ValueError
        If an atom line cannot be parsed.
    """

    line_end = _next_line_end(data, offset)

    if line_end < 0:
        if at_eof:
            raise EOFError("incomplete frame")

        return (STATUS_NEED_MORE, None, None, offset)

    pos = line_end + 1

    line_bounds = []
    scan = pos

    for _ in range(n_atoms):
        line_end = _next_line_end(data, scan)

        if line_end < 0:
            if at_eof:
                raise EOFError("incomplete frame")

            return (STATUS_NEED_MORE, None, None, offset)

        line_bounds.append((scan, line_end))
        scan = line_end + 1

    first_name = None

    if want_first_name and n_atoms > 0:
        start, end = line_bounds[0]
        first_fields = data[start:end].split(None, 1)
        first_name = first_fields[0] if first_fields else b""

    if mode == MODE_XYZ:
        lines = [
            data[start:end].decode("utf-8") for start, end in line_bounds
        ]
        values = process_lines(lines, n_atoms)
    else:
        values = np.empty(n_atoms, dtype=np.float64)

        for i, (start, end) in enumerate(line_bounds):
            fields = data[start:end].split()

            if len(fields) != 2:
                raise ValueError("Could not parse line")

            values[i] = float(fields[1])

    return (STATUS_FRAME, values, first_name, scan)
