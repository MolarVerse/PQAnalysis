# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: cdivision=True
"""
Cython kernel for the MSD accumulation hot loop.

The kernel advances the running Diffcalc-style MSD state by one
trajectory frame: it gathers the selected atom positions of the frame
(float32 -> float64, exact), unwraps them with a running minimum image
convention applied to the per-frame displacement vectors, spawns,
replaces and drains time origins following the legacy Diffcalc
bookkeeping and accumulates the per-axis squared displacements of all
active time origins into the lag bins of the msd accumulator.

The unwrapping uses the exact op semantics of the pure numpy
implementation: ``fractional = d @ B^-1.T`` and
``shift += -rint(fractional) @ B.T``, where the C ``rint`` rounds
half-even under the default rounding mode, matching ``np.rint``.

A pure Python/numpy fallback with the identical signature lives in
:py:mod:`PQAnalysis.analysis.msd._msd_kernel_py`.
"""

import numpy as np

cimport numpy as np

from libc.math cimport rint
from libc.string cimport memcpy, memmove


def msd_frame_update(
    const np.float32_t[:, ::1] values,
    const np.int64_t[::1] indices,
    const np.float64_t[:, ::1] box,
    const np.float64_t[:, ::1] inv_box,
    long long is_vacuum,
    np.float64_t[:, ::1] pos,
    np.float64_t[:, ::1] prev_pos,
    np.float64_t[:, ::1] shift,
    np.float64_t[:, ::1] unwrapped,
    np.float64_t[:, :, ::1] origins,
    np.float64_t[:, ::1] msd,
    np.int64_t[::1] state,
    long long counter,
    long long gap,
    long long window,
    long long n_start,
    long long stop_frame,
):
    """
    Advances the running MSD state by one trajectory frame.

    The selected rows of ``values`` are gathered into ``pos`` (exact
    float32 -> float64 conversion). For every frame after the first
    one the per-frame displacement ``pos - prev_pos`` is folded back
    into the minimum image convention and the resulting change is
    accumulated into ``shift`` (skipped for vacuum cells). The
    unwrapped coordinates ``pos + shift`` are stored in ``unwrapped``
    and ``prev_pos`` is updated to ``pos``. Frames with
    ``counter < n_start`` only update the unwrapping state. All other
    frames run the legacy Diffcalc time origin bookkeeping (spawn,
    final-lag accumulation, replace or drain every ``gap`` frames) and
    accumulate the per-axis squared displacements of all active time
    origins into the lag bins of ``msd``.

    Parameters
    ----------
    values : np.float32 array of shape (n_atoms, 3), C-contiguous
        The raw frame values (positions) of all atoms of the frame.
    indices : np.int64 array of shape (n_sel,)
        The indices of the selected atoms.
    box : np.float64 array of shape (3, 3), C-contiguous
        The box matrix of the current frame. Ignored for vacuum
        cells.
    inv_box : np.float64 array of shape (3, 3), C-contiguous
        The inverse box matrix of the current frame. Ignored for
        vacuum cells.
    is_vacuum : int
        Whether the cell of the current frame is a vacuum cell (no
        unwrapping is applied then).
    pos : np.float64 array of shape (n_sel, 3), C-contiguous
        Scratch buffer, filled with the selected positions of the
        current frame.
    prev_pos : np.float64 array of shape (n_sel, 3), C-contiguous
        The selected positions of the previous frame; updated in
        place to the current frame. Its input content is ignored for
        ``counter == 1``.
    shift : np.float64 array of shape (n_sel, 3), C-contiguous
        The cumulative unwrapping shift vectors; updated in place.
    unwrapped : np.float64 array of shape (n_sel, 3), C-contiguous
        Scratch buffer, filled with the unwrapped coordinates
        ``pos + shift``.
    origins : np.float64 array of shape (n_origins_max, n_sel, 3), C-contiguous
        The unwrapped coordinates of the active time origins, oldest
        origin at index 0; updated in place.
    msd : np.float64 array of shape (window + 1, 3), C-contiguous
        The raw (unnormalized) per-axis MSD accumulator; updated in
        place.
    state : np.int64 array of shape (2,)
        The origin bookkeeping state ``[n_active, last]``; updated in
        place.
    counter : int
        The 1-based frame counter of the current frame.
    gap : int
        The gap between two time origins in frames.
    window : int
        The correlation window size in frames.
    n_start : int
        The first frame (1-based frame counter) at which processing
        starts. Earlier frames only update the unwrapping state.
    stop_frame : int
        The last frame (1-based frame counter) at which a time origin
        may spawn.
    """

    cdef Py_ssize_t n_sel = indices.shape[0]
    cdef Py_ssize_t n_origins_max = origins.shape[0]
    cdef long long n_active = state[0]
    cdef long long last = state[1]
    cdef bint do_unwrap = counter > 1 and is_vacuum == 0
    cdef Py_ssize_t a, o, row
    cdef long long lag
    cdef double p0, p1, p2
    cdef double d0, d1, d2
    cdef double f0, f1, f2
    cdef double r0, r1, r2
    cdef double sx, sy, sz
    cdef double dx, dy, dz

    # gather the selected positions (exact float32 -> float64) and
    # unwrap them with the running minimum image convention
    for a in range(n_sel):
        row = <Py_ssize_t> indices[a]
        p0 = <double> values[row, 0]
        p1 = <double> values[row, 1]
        p2 = <double> values[row, 2]

        if do_unwrap:
            d0 = p0 - prev_pos[a, 0]
            d1 = p1 - prev_pos[a, 1]
            d2 = p2 - prev_pos[a, 2]

            f0 = d0 * inv_box[0, 0] + d1 * inv_box[0, 1] + d2 * inv_box[0, 2]
            f1 = d0 * inv_box[1, 0] + d1 * inv_box[1, 1] + d2 * inv_box[1, 2]
            f2 = d0 * inv_box[2, 0] + d1 * inv_box[2, 1] + d2 * inv_box[2, 2]

            r0 = rint(f0)
            r1 = rint(f1)
            r2 = rint(f2)

            # adding an all-zero shift change never changes the
            # accumulated shift, so it can be skipped
            if r0 != 0.0 or r1 != 0.0 or r2 != 0.0:
                shift[a, 0] -= r0 * box[0, 0] + r1 * box[0, 1] + r2 * box[0, 2]
                shift[a, 1] -= r0 * box[1, 0] + r1 * box[1, 1] + r2 * box[1, 2]
                shift[a, 2] -= r0 * box[2, 0] + r1 * box[2, 1] + r2 * box[2, 2]

        pos[a, 0] = p0
        pos[a, 1] = p1
        pos[a, 2] = p2

        prev_pos[a, 0] = p0
        prev_pos[a, 1] = p1
        prev_pos[a, 2] = p2

        unwrapped[a, 0] = p0 + shift[a, 0]
        unwrapped[a, 1] = p1 + shift[a, 1]
        unwrapped[a, 2] = p2 + shift[a, 2]

    if counter < n_start:
        return

    if counter % gap == 0:
        if n_active != <long long> n_origins_max and counter <= stop_frame:
            # starting stage - add new origin
            if n_active == 0:
                last = counter

            memcpy(
                &origins[n_active, 0, 0],
                &unwrapped[0, 0],
                n_sel * 3 * sizeof(double),
            )
            n_active += 1

        elif n_active > 0 and last + window == counter:
            # oldest origin reached the full window:
            # accumulate its final lag term
            sx = 0.0
            sy = 0.0
            sz = 0.0

            for a in range(n_sel):
                dx = unwrapped[a, 0] - origins[0, a, 0]
                dy = unwrapped[a, 1] - origins[0, a, 1]
                dz = unwrapped[a, 2] - origins[0, a, 2]
                sx += dx * dx
                sy += dy * dy
                sz += dz * dz

            msd[window, 0] += sx
            msd[window, 1] += sy
            msd[window, 2] += sz

            if n_active > 1:
                memmove(
                    &origins[0, 0, 0],
                    &origins[1, 0, 0],
                    (n_active - 1) * n_sel * 3 * sizeof(double),
                )

            if counter > stop_frame:
                # stopping stage - drop without replacement
                n_active -= 1
            else:
                # running stage - replace by a new origin
                memcpy(
                    &origins[n_active - 1, 0, 0],
                    &unwrapped[0, 0],
                    n_sel * 3 * sizeof(double),
                )

            last += gap

        state[0] = n_active
        state[1] = last

    # accumulate the squared displacements of all active origins
    for o in range(<Py_ssize_t> n_active):
        lag = counter - last - gap * <long long> o

        sx = 0.0
        sy = 0.0
        sz = 0.0

        for a in range(n_sel):
            dx = unwrapped[a, 0] - origins[o, a, 0]
            dy = unwrapped[a, 1] - origins[o, a, 1]
            dz = unwrapped[a, 2] - origins[o, a, 2]
            sx += dx * dx
            sy += dy * dy
            sz += dz * dz

        msd[lag, 0] += sx
        msd[lag, 1] += sy
        msd[lag, 2] += sz
