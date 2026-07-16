# cython: boundscheck=False, wraparound=False, cdivision=True
"""
Cython kernels for the hot loops of the VACF analysis.

The kernels keep the exact numeric semantics of the analysis: the
velocities are parsed as float32 by the trajectory readers and all
accumulation is performed in float64. A numpy fallback with identical
signatures lives in
:py:mod:`PQAnalysis.analysis.vacf._vacf_kernel_py` and is used when
this extension is not available.

The per-origin dot products of :py:func:`accumulate_frame` are
accumulated in float64 with a fixed 4-way unrolled summation order.
The numpy fallback uses ``np.einsum``/``np.sum`` (SIMD/pairwise
summation) instead, so kernel and fallback can differ by floating
point rounding on the order of the machine epsilon - far below the
parity tolerances of the analysis. All other kernels are bitwise
identical to the fallback.
"""

import numpy as np

from libc.stdio cimport sscanf
from libc.string cimport memcpy, memmove



cdef inline double _dot(
    const double* a,
    const double* b,
    Py_ssize_t n,
) noexcept nogil:
    """
    Dot product of two float64 buffers of length ``n``.

    The summation is 4-way unrolled with a fixed combination order
    ``(s0 + s1) + (s2 + s3)`` so that the result is deterministic on
    every platform.
    """
    cdef double s0 = 0.0
    cdef double s1 = 0.0
    cdef double s2 = 0.0
    cdef double s3 = 0.0
    cdef Py_ssize_t k = 0

    while k + 4 <= n:
        s0 = s0 + a[k] * b[k]
        s1 = s1 + a[k + 1] * b[k + 1]
        s2 = s2 + a[k + 2] * b[k + 2]
        s3 = s3 + a[k + 3] * b[k + 3]
        k = k + 4

    while k < n:
        s0 = s0 + a[k] * b[k]
        k = k + 1

    return (s0 + s1) + (s2 + s3)


def accumulate_frame(
    double[::1] corr,
    double[:, :, ::1] origin_vel,
    double[::1] origin_norm,
    long long[::1] origin_frame,
    Py_ssize_t n_active,
    const double[:, ::1] vel,
    long long frame_number,
    bint spawn,
    long long window_size,
):
    """
    Performs the per-frame update of the sliding-origin estimator.

    If ``spawn`` is set, the frame is registered as a new time origin
    (velocities, aggregate squared norm and frame number are stored in
    the origin bookkeeping arrays). Afterwards every active origin
    ``i`` contributes ``sum_j v_j(t) . v_j(t0_i) / sum_j |v_j(t0_i)|^2``
    to the lag ``t - t0_i`` of ``corr``. The oldest origin is retired
    (shifted out of the bookkeeping arrays) after it has contributed
    to the lag ``window_size``.

    Parameters
    ----------
    corr : np.ndarray of float64, shape (window_size + 1,)
        The lag accumulator, updated in place.
    origin_vel : np.ndarray of float64, shape (n_slots, n_target, 3)
        The velocities of the active origins, updated in place.
    origin_norm : np.ndarray of float64, shape (n_slots,)
        The aggregate squared velocity norms of the active origins,
        updated in place.
    origin_frame : np.ndarray of int64, shape (n_slots,)
        The 1-based frame numbers of the active origins, updated in
        place.
    n_active : int
        The number of active origins before this frame.
    vel : np.ndarray of float64, shape (n_target, 3)
        The (charge weighted) selected velocities of the frame.
    frame_number : int
        The 1-based number of the frame.
    spawn : bool
        Whether a new time origin is spawned at this frame.
    window_size : int
        The correlation window length in frames.

    Returns
    -------
    int
        The number of active origins after this frame, or ``-1`` if a
        new origin was to be spawned but its aggregate squared
        velocity norm is zero (in which case no state was modified).
    """
    cdef Py_ssize_t m = vel.shape[0]
    cdef Py_ssize_t n = m * 3
    cdef const double* v = &vel[0, 0]
    cdef double* buffer = &origin_vel[0, 0, 0]
    cdef double norm
    cdef double scalar
    cdef Py_ssize_t o
    cdef long long lag

    if spawn:
        norm = _dot(v, v, n)

        if norm == 0.0:
            return -1

        memcpy(buffer + n_active * n, v, n * sizeof(double))
        origin_norm[n_active] = norm
        origin_frame[n_active] = frame_number
        n_active = n_active + 1

    if n_active == 0:
        return 0

    for o in range(n_active):
        scalar = _dot(buffer + o * n, v, n)
        lag = frame_number - origin_frame[o]
        corr[lag] += scalar / origin_norm[o]

    if frame_number - origin_frame[0] == window_size:
        # retire the oldest origin
        memmove(buffer, buffer + n, (n_active - 1) * n * sizeof(double))

        for o in range(n_active - 1):
            origin_norm[o] = origin_norm[o + 1]
            origin_frame[o] = origin_frame[o + 1]

        n_active = n_active - 1

    return n_active


def weight_frame(
    const float[:, ::1] values,
    const Py_ssize_t[::1] indices,
    charges,
):
    """
    Builds the (charge weighted) selected float64 velocities of a
    frame from the raw float32 values of the trajectory reader.

    Bitwise identical to
    ``np.asarray(values, dtype=np.float64)[indices] * charges[:, None]``
    (the float32 -> float64 cast is exact and the weighting is a
    plain elementwise IEEE product).

    Parameters
    ----------
    values : np.ndarray of float32, shape (n_atoms, 3)
        The raw values of the frame.
    indices : np.ndarray of intp, shape (n_target,)
        The indices of the selected atoms.
    charges : np.ndarray of float64, shape (n_target,), or None
        The charges of the selected atoms, or None for no weighting.

    Returns
    -------
    np.ndarray of float64, shape (n_target, 3)
        The (charge weighted) selected velocities of the frame.
    """
    cdef Py_ssize_t m = indices.shape[0]

    out_array = np.empty((m, 3), dtype=np.float64)

    cdef double[:, ::1] out = out_array
    cdef const double[::1] q
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef double q_i

    if charges is None:
        for i in range(m):
            j = indices[i]
            out[i, 0] = <double> values[j, 0]
            out[i, 1] = <double> values[j, 1]
            out[i, 2] = <double> values[j, 2]
    else:
        q = charges

        for i in range(m):
            j = indices[i]
            q_i = q[i]
            out[i, 0] = <double> values[j, 0] * q_i
            out[i, 1] = <double> values[j, 1] * q_i
            out[i, 2] = <double> values[j, 2] * q_i

    return out_array


def parse_charge_lines(lines, Py_ssize_t n_atoms):
    """
    Parses the float64 charge values of 'name charge' body lines.

    Every line must consist of exactly two whitespace separated
    tokens; the second token is parsed as a float64 (via the
    correctly rounded C ``strtod``, bitwise identical to Python's
    ``float``).

    Parameters
    ----------
    lines : list of str
        The atom lines of the charge frame body.
    n_atoms : int
        The number of atoms in the frame.

    Returns
    -------
    np.ndarray of float64, shape (n_atoms,)
        The charge values parsed from the lines.

    Raises
    ------
    ValueError
        If a line does not consist of a name token followed by
        exactly one parsable charge value.
    """
    charges_array = np.empty(n_atoms, dtype=np.float64)

    cdef double[::1] charges = charges_array
    cdef double value
    cdef int n_parsed
    cdef int end
    cdef Py_ssize_t i
    cdef const char* c_line
    cdef char c

    for i in range(n_atoms):
        line_bytes = lines[i].encode('utf-8')
        c_line = line_bytes

        end = -1
        n_parsed = sscanf(c_line, "%*s %lf%n", &value, &end)

        if n_parsed != 1 or end < 0:
            raise ValueError("Could not parse line")

        # exactly-two-token semantics: the rest of the line has to be
        # whitespace only
        c = c_line[end]
        while c != 0:
            if not (c == 32 or (c >= 9 and c <= 13)):
                raise ValueError("Could not parse line")
            end = end + 1
            c = c_line[end]

        charges[i] = value

    return charges_array
