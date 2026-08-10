# cython: boundscheck=False, wraparound=False, cdivision=True
"""Cython kernel for legacy-compatible total momentum norms.

The operation order matches ``thh_tools/thh_momentum/equipartition.jl``:
each float64 velocity component is multiplied by its atom mass and added
to the corresponding momentum component in atom order. The three-vector
norm follows Julia's ``generic_norm2`` order before applying ``scale``.
"""

import numpy as np
cimport numpy as np

from libc.math cimport fabs, isfinite, isinf, isnan, sqrt
from libc.string cimport memchr

cdef extern from "<stdlib.h>":
    double strtod(const char *nptr, char **endptr)
    long strtol(const char *nptr, char **endptr, int base)


cdef inline bint _is_space(char value) noexcept nogil:
    """C-locale whitespace accepted by the trajectory parser."""
    return value == c' ' or (value >= c'\t' and value <= c'\r')


cdef inline double _legacy_norm3(
    double x,
    double y,
    double z,
) noexcept nogil:
    """Replicates Julia's ``generic_norm2`` for three float64 values."""
    cdef double max_abs = fabs(x)
    cdef double value = fabs(y)
    cdef double square
    cdef double total

    if value > max_abs:
        max_abs = value

    value = fabs(z)
    if value > max_abs:
        max_abs = value

    if max_abs == 0.0 or isinf(max_abs) or isnan(max_abs):
        return max_abs

    square = max_abs * max_abs

    if isfinite((3.0 * max_abs) * max_abs) and square != 0.0:
        total = x * x
        total = total + y * y
        total = total + z * z
        return sqrt(total)

    value = fabs(x) / max_abs
    total = value * value
    value = fabs(y) / max_abs
    total = total + value * value
    value = fabs(z) / max_abs
    total = total + value * value

    return max_abs * sqrt(total)


def legacy_momentum_norm(
    const np.float64_t[:, ::1] values,
    const np.int64_t[::1] indices,
    const np.float64_t[::1] masses,
    double scale,
):
    """Calculates one scaled momentum norm in legacy operation order."""
    cdef Py_ssize_t n_selected = indices.shape[0]
    cdef Py_ssize_t selected_index
    cdef Py_ssize_t row
    cdef double mass
    cdef double momentum_x = 0.0
    cdef double momentum_y = 0.0
    cdef double momentum_z = 0.0

    if masses.shape[0] != n_selected:
        raise ValueError("indices and masses must have the same length")

    with nogil:
        for selected_index in range(n_selected):
            row = <Py_ssize_t> indices[selected_index]
            mass = masses[selected_index]
            momentum_x = momentum_x + mass * values[row, 0]
            momentum_y = momentum_y + mass * values[row, 1]
            momentum_z = momentum_z + mass * values[row, 2]

    return _legacy_norm3(momentum_x, momentum_y, momentum_z) * scale


def legacy_momentum_file(
    bytes data,
    Py_ssize_t n_atoms,
    bint strip_first,
    const np.int64_t[::1] indices,
    const np.float64_t[::1] masses,
    double scale,
):
    """Parse and reduce a fixed-topology velocity file in one pass."""
    cdef Py_ssize_t n_selected = indices.shape[0]
    if masses.shape[0] != n_selected:
        raise ValueError("indices and masses must have the same length")

    if n_atoms < 0:
        raise ValueError("Invalid atom count.")

    cdef Py_ssize_t raw_n_atoms = n_atoms + <Py_ssize_t> strip_first
    values_array = np.empty((n_atoms, 3), dtype=np.float64)
    cdef double[:, ::1] values = values_array

    norms = []
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
    cdef Py_ssize_t atom
    cdef Py_ssize_t output_atom
    cdef Py_ssize_t component
    cdef Py_ssize_t selected_index
    cdef Py_ssize_t row
    cdef const char* found
    cdef char* endptr
    cdef long frame_atoms
    cdef double mass
    cdef double value
    cdef double momentum_x
    cdef double momentum_y
    cdef double momentum_z

    while True:
        # Locate the next non-blank header line.
        while True:
            if pos >= n_data:
                return (
                    np.asarray(norms, dtype=np.float64),
                    box_headers,
                    first_names,
                )

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

        if frame_atoms != raw_n_atoms:
            raise ValueError("Frame atom count does not match the topology.")

        box_headers.append(data[box_start:line_end])
        pos = line_end + 1

        # Comment line.
        found = <const char*> memchr(base + pos, c'\n', n_data - pos)
        if found == NULL:
            raise EOFError("incomplete frame")
        pos = (found - base) + 1

        # Preserve the streaming parser's incomplete-frame error order.
        scan = pos
        for atom in range(raw_n_atoms):
            found = <const char*> memchr(base + scan, c'\n', n_data - scan)
            if found == NULL:
                raise EOFError("incomplete frame")
            scan = (found - base) + 1

        for atom in range(raw_n_atoms):
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

                value = strtod(base + pos, &endptr)
                if <const char*> endptr == base + pos:
                    raise ValueError("Could not parse line")

                if not strip_first or atom > 0:
                    values[output_atom, component] = value

                pos = <const char*> endptr - base

            pos = line_end + 1

        pos = scan
        momentum_x = 0.0
        momentum_y = 0.0
        momentum_z = 0.0

        for selected_index in range(n_selected):
            row = <Py_ssize_t> indices[selected_index]
            mass = masses[selected_index]
            momentum_x = momentum_x + mass * values[row, 0]
            momentum_y = momentum_y + mass * values[row, 1]
            momentum_z = momentum_z + mass * values[row, 2]

        norms.append(
            _legacy_norm3(momentum_x, momentum_y, momentum_z) * scale
        )
