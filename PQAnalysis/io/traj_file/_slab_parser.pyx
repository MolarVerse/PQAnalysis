# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
"""
A Cython byte-slab frame parser for xyz-family trajectory files.

This module implements the same function contract as the pure Python
fallback :py:mod:`PQAnalysis.io.traj_file._slab_parser_py` (see its
module docstring for the contract), but parses the frames directly
from the byte buffer: newlines are located with ``memchr`` and the
numeric tokens are converted in place with ``strtol``/``strtof``/
``strtod`` without decoding or splitting the lines.

The default float32 xyz values are parsed with ``strtof`` (single
rounding), which is bitwise identical to the ``sscanf("%f")``
conversion of the line-based readers. The float64 xyz mode and charge
values are parsed directly with ``strtod`` so no precision is discarded.
"""

import numpy as np

from libc.string cimport memchr

from PQAnalysis.io.traj_file._slab_parser_py import (
    MODE_CHARGE,
    MODE_XYZ,
    MODE_XYZ64,
    STATUS_BAD_HEADER,
    STATUS_EOF,
    STATUS_FRAME,
    STATUS_NEED_MORE,
)

cdef extern from "<stdlib.h>":
    float strtof(const char *nptr, char **endptr) nogil
    double strtod(const char *nptr, char **endptr) nogil
    long strtol(const char *nptr, char **endptr, int base) nogil

# C-level copies of the shared status/mode constants of
# _slab_parser_py (the single source of truth for their values)
cdef int _MODE_XYZ = MODE_XYZ
cdef int _MODE_CHARGE = MODE_CHARGE
cdef int _MODE_XYZ64 = MODE_XYZ64


cdef inline bint _is_space(char c) nogil:
    """C locale whitespace as used by sscanf/strtof token skipping."""
    return c == c' ' or (c >= c'\t' and c <= c'\r')


def scan_header(
    bytes data,
    Py_ssize_t offset,
    bint at_eof,
    Py_ssize_t forced_n_atoms=-1,
):
    """
    Scans for the header line of the next frame.

    Same contract as
    :py:func:`PQAnalysis.io.traj_file._slab_parser_py.scan_header`.

    Parameters
    ----------
    data : bytes
        The buffer to scan.
    offset : int
        The offset to start scanning at.
    at_eof : bool
        Whether the buffer contains the complete rest of the file.
    forced_n_atoms : int, optional
        The atom count to use when the count token of the header
        line is not a plain integer literal but was successfully
        parsed by the caller. Default is -1 (disabled).

    Returns
    -------
    tuple
        ``(status, n_atoms, box_bytes, header_token, next_offset)``,
        see the fallback implementation for the exact semantics.
    """

    cdef const char* base = data
    cdef Py_ssize_t n_data = len(data)
    cdef Py_ssize_t pos = offset
    cdef Py_ssize_t line_end, i, tok_start, tok_end
    cdef const char* found
    cdef char* endptr
    cdef long count

    while True:
        if pos >= n_data:
            if at_eof:
                return (STATUS_EOF, -1, None, None, pos)

            return (STATUS_NEED_MORE, -1, None, None, offset)

        found = <const char*> memchr(base + pos, c'\n', n_data - pos)

        if found == NULL:
            if not at_eof:
                return (STATUS_NEED_MORE, -1, None, None, offset)

            line_end = n_data
        else:
            line_end = found - base

        i = pos

        while i < line_end and _is_space(base[i]):
            i += 1

        if i < line_end:
            break

        pos = line_end + 1

    tok_start = i

    while i < line_end and not _is_space(base[i]):
        i += 1

    tok_end = i

    while i < line_end and _is_space(base[i]):
        i += 1

    box_bytes = data[i:line_end]
    body_offset = line_end + 1

    count = strtol(base + tok_start, &endptr, 10)

    if <const char*> endptr == base + tok_end and count >= 0:
        return (STATUS_FRAME, count, box_bytes, None, body_offset)

    if forced_n_atoms >= 0:
        return (STATUS_FRAME, forced_n_atoms, box_bytes, None, body_offset)

    return (
        STATUS_BAD_HEADER,
        -1,
        box_bytes,
        data[tok_start:tok_end],
        body_offset,
    )


def parse_body(
    bytes data,
    Py_ssize_t offset,
    Py_ssize_t n_atoms,
    bint at_eof,
    bint want_first_name,
    int mode,
):
    """
    Parses the body (comment line plus atom lines) of a frame.

    Same contract as
    :py:func:`PQAnalysis.io.traj_file._slab_parser_py.parse_body`.

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
        Whether to extract the name token of the first atom line.
    mode : int
        The body mode: ``MODE_XYZ``, ``MODE_XYZ64`` or
        ``MODE_CHARGE``.

    Returns
    -------
    tuple
        ``(status, values, first_name, next_offset)``, see the
        fallback implementation for the exact semantics.

    Raises
    ------
    EOFError
        If the buffer contains the complete rest of the file and the
        frame is incomplete.
    ValueError
        If an atom line cannot be parsed.
    """

    cdef const char* base = data
    cdef Py_ssize_t n_data = len(data)
    cdef const char* found
    cdef Py_ssize_t pos, scan, i

    # comment line
    found = NULL

    if offset < n_data:
        found = <const char*> memchr(base + offset, c'\n', n_data - offset)

    if found == NULL:
        if at_eof:
            raise EOFError("incomplete frame")

        return (STATUS_NEED_MORE, None, None, offset)

    pos = (found - base) + 1

    # pre-scan: all atom lines must be complete before any of them is
    # parsed, so that a truncated frame at the end of the file is
    # always reported as incomplete (EOFError) even if it also
    # contains malformed lines
    scan = pos

    for i in range(n_atoms):
        found = NULL

        if scan < n_data:
            found = <const char*> memchr(base + scan, c'\n', n_data - scan)

        if found == NULL:
            if at_eof:
                raise EOFError("incomplete frame")

            return (STATUS_NEED_MORE, None, None, offset)

        scan = (found - base) + 1

    if mode == _MODE_XYZ:
        values, first_name = _parse_xyz_lines(
            data, pos, n_atoms, want_first_name
        )
    elif mode == _MODE_XYZ64:
        values, first_name = _parse_xyz64_lines(
            data, pos, n_atoms, want_first_name
        )
    else:
        values, first_name = _parse_charge_lines(
            data, pos, n_atoms, want_first_name
        )

    return (STATUS_FRAME, values, first_name, scan)


def parse_xyz_frames(
    bytes data,
    Py_ssize_t n_frames,
    Py_ssize_t n_atoms,
    bint strip_first,
    int mode,
):
    """Parse a complete fixed-topology xyz-family file in one call.

    The numeric token conversions and frame order are identical to
    repeated :func:`parse_body` calls, but all frame values are written
    directly into one three-dimensional NumPy array.
    """
    if mode != _MODE_XYZ and mode != _MODE_XYZ64:
        raise ValueError("Batch parsing supports only xyz-family data.")

    if n_frames < 0 or n_atoms < 0 or (strip_first and n_atoms == 0):
        raise ValueError("Invalid batch shape.")

    cdef Py_ssize_t n_output_atoms = n_atoms - <Py_ssize_t> strip_first

    if mode == _MODE_XYZ:
        values = np.empty(
            (n_frames, n_output_atoms, 3),
            dtype=np.float32,
        )
    else:
        values = np.empty(
            (n_frames, n_output_atoms, 3),
            dtype=np.float64,
        )

    cdef float[:, :, ::1] out32
    cdef double[:, :, ::1] out64

    if mode == _MODE_XYZ:
        out32 = values
    else:
        out64 = values

    box_headers = []
    first_names = []

    cdef const char* base = data
    cdef Py_ssize_t n_data = len(data)
    cdef Py_ssize_t pos = 0
    cdef Py_ssize_t line_end
    cdef Py_ssize_t scan
    cdef Py_ssize_t token_start
    cdef Py_ssize_t token_end
    cdef Py_ssize_t box_start
    cdef Py_ssize_t name_start
    cdef Py_ssize_t name_end
    cdef Py_ssize_t frame
    cdef Py_ssize_t atom
    cdef Py_ssize_t output_atom
    cdef Py_ssize_t component
    cdef const char* found
    cdef char* endptr
    cdef long frame_atoms
    cdef float value32
    cdef double value64

    for frame in range(n_frames):
        # Skip blank lines before the frame header.
        while True:
            if pos >= n_data:
                raise EOFError("incomplete frame")

            found = <const char*> memchr(base + pos, c'\n', n_data - pos)
            if found == NULL:
                raise EOFError("incomplete frame")

            line_end = found - base
            scan = pos

            while scan < line_end and _is_space(base[scan]):
                scan += 1

            if scan < line_end:
                break

            pos = line_end + 1

        token_start = scan
        while scan < line_end and not _is_space(base[scan]):
            scan += 1
        token_end = scan

        while scan < line_end and _is_space(base[scan]):
            scan += 1
        box_start = scan

        frame_atoms = strtol(base + token_start, &endptr, 10)
        if (
            <const char*> endptr != base + token_end
            or frame_atoms < 0
        ):
            raise ValueError("Invalid frame atom count.")

        if frame_atoms != n_atoms:
            raise ValueError("Frame atom count does not match the batch.")

        box_headers.append(data[box_start:line_end])
        pos = line_end + 1

        # Comment line.
        found = <const char*> memchr(base + pos, c'\n', n_data - pos)
        if found == NULL:
            raise EOFError("incomplete frame")
        pos = (found - base) + 1

        # Check completeness before parsing any atom row, matching
        # parse_body's error ordering for truncated frames.
        scan = pos
        for atom in range(n_atoms):
            found = <const char*> memchr(base + scan, c'\n', n_data - scan)
            if found == NULL:
                raise EOFError("incomplete frame")
            scan = (found - base) + 1

        for atom in range(n_atoms):
            found = <const char*> memchr(base + pos, c'\n', n_data - pos)
            line_end = found - base

            while pos < line_end and _is_space(base[pos]):
                pos += 1

            if pos == line_end:
                raise ValueError("Could not parse line")

            name_start = pos
            while pos < line_end and not _is_space(base[pos]):
                pos += 1
            name_end = pos

            if strip_first and atom == 0:
                first_names.append(data[name_start:name_end])

            output_atom = atom - <Py_ssize_t> strip_first

            for component in range(3):
                while pos < line_end and _is_space(base[pos]):
                    pos += 1

                if pos == line_end:
                    raise ValueError("Could not parse line")

                if mode == _MODE_XYZ:
                    value32 = strtof(base + pos, &endptr)
                else:
                    value64 = strtod(base + pos, &endptr)

                if <const char*> endptr == base + pos:
                    raise ValueError("Could not parse line")

                if not strip_first or atom > 0:
                    if mode == _MODE_XYZ:
                        out32[frame, output_atom, component] = value32
                    else:
                        out64[frame, output_atom, component] = value64

                pos = <const char*> endptr - base

            pos = line_end + 1

        pos = scan

    while pos < n_data and _is_space(base[pos]):
        pos += 1

    if pos != n_data:
        raise ValueError("File contains more frames than expected.")

    return values, box_headers, first_names


cdef _parse_xyz_lines(
    bytes data,
    Py_ssize_t pos,
    Py_ssize_t n_atoms,
    bint want_first_name,
):
    """
    Parses ``n_atoms`` xyz-family atom lines into a float32 array.

    Every line must consist of a name token followed by at least
    three float values; the floats are parsed with ``strtof``
    (bitwise identical to the ``sscanf("%f")`` conversions of
    ``process_lines``).
    """

    cdef const char* base = data
    cdef Py_ssize_t n_data = len(data)

    values = np.empty((n_atoms, 3), dtype=np.float32)

    cdef float[:, ::1] out = values
    cdef Py_ssize_t p = pos
    cdef Py_ssize_t line_end, tok_start, i, j
    cdef const char* found
    cdef char* endptr
    cdef float value

    first_name = None

    for i in range(n_atoms):
        # guaranteed by the pre-scan of parse_body
        found = <const char*> memchr(base + p, c'\n', n_data - p)
        line_end = found - base

        while p < line_end and _is_space(base[p]):
            p += 1

        if p == line_end:
            raise ValueError("Could not parse line")

        tok_start = p

        while p < line_end and not _is_space(base[p]):
            p += 1

        if want_first_name and i == 0:
            first_name = data[tok_start:p]

        for j in range(3):
            while p < line_end and _is_space(base[p]):
                p += 1

            if p == line_end:
                raise ValueError("Could not parse line")

            value = strtof(base + p, &endptr)

            if <const char*> endptr == base + p:
                raise ValueError("Could not parse line")

            out[i, j] = value
            p = <const char*> endptr - base

        p = line_end + 1

    return values, first_name


cdef _parse_xyz64_lines(
    bytes data,
    Py_ssize_t pos,
    Py_ssize_t n_atoms,
    bint want_first_name,
):
    """Parse xyz-family atom lines directly into a float64 array."""
    cdef const char* base = data
    cdef Py_ssize_t n_data = len(data)

    values = np.empty((n_atoms, 3), dtype=np.float64)

    cdef double[:, ::1] out = values
    cdef Py_ssize_t p = pos
    cdef Py_ssize_t line_end, tok_start, i, j
    cdef const char* found
    cdef char* endptr
    cdef double value

    first_name = None

    for i in range(n_atoms):
        # guaranteed by the pre-scan of parse_body
        found = <const char*> memchr(base + p, c'\n', n_data - p)
        line_end = found - base

        while p < line_end and _is_space(base[p]):
            p += 1

        if p == line_end:
            raise ValueError("Could not parse line")

        tok_start = p

        while p < line_end and not _is_space(base[p]):
            p += 1

        if want_first_name and i == 0:
            first_name = data[tok_start:p]

        for j in range(3):
            while p < line_end and _is_space(base[p]):
                p += 1

            if p == line_end:
                raise ValueError("Could not parse line")

            value = strtod(base + p, &endptr)

            if <const char*> endptr == base + p:
                raise ValueError("Could not parse line")

            out[i, j] = value
            p = <const char*> endptr - base

        p = line_end + 1

    return values, first_name


cdef _parse_charge_lines(
    bytes data,
    Py_ssize_t pos,
    Py_ssize_t n_atoms,
    bint want_first_name,
):
    """
    Parses ``n_atoms`` charge atom lines into a float64 array.

    Every line must consist of a name token followed by exactly one
    float value; the floats are parsed with ``strtod`` (bitwise
    identical to Python's correctly rounded ``float``).
    """

    cdef const char* base = data
    cdef Py_ssize_t n_data = len(data)

    values = np.empty(n_atoms, dtype=np.float64)

    cdef double[::1] out = values
    cdef Py_ssize_t p = pos
    cdef Py_ssize_t line_end, tok_start, i
    cdef const char* found
    cdef char* endptr
    cdef double value

    first_name = None

    for i in range(n_atoms):
        # guaranteed by the pre-scan of parse_body
        found = <const char*> memchr(base + p, c'\n', n_data - p)
        line_end = found - base

        while p < line_end and _is_space(base[p]):
            p += 1

        if p == line_end:
            raise ValueError("Could not parse line")

        tok_start = p

        while p < line_end and not _is_space(base[p]):
            p += 1

        if want_first_name and i == 0:
            first_name = data[tok_start:p]

        while p < line_end and _is_space(base[p]):
            p += 1

        if p == line_end:
            raise ValueError("Could not parse line")

        value = strtod(base + p, &endptr)

        if <const char*> endptr == base + p:
            raise ValueError("Could not parse line")

        p = <const char*> endptr - base

        # exactly-two-token semantics: the rest of the line has to be
        # whitespace only
        while p < line_end and _is_space(base[p]):
            p += 1

        if p != line_end:
            raise ValueError("Could not parse line")

        out[i] = value
        p = line_end + 1

    return values, first_name
