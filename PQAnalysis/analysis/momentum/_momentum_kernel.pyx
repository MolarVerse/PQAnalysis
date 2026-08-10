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
