# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: cdivision=True
"""
Cython kernel for the RDF distance-histogram hot loop.

The kernel accumulates the minimum image distance histogram of one
trajectory frame: for every reference atom the distances to all
target atoms (excluding the reference atom itself) are computed and
scattered into the bins ``floor((distance - r_min) / delta_r)``.

The numeric semantics replicate the original numpy hot loop bit for
bit: the pair displacements are computed in float32 (the dtype of
the raw frame values) and widened to float64, the minimum image
convention uses the exact op order of
:py:meth:`PQAnalysis.core.cell.cell.Cell.image` (an orthorhombic
box-length branch with round-half-even ``rint`` matching
``np.round``, and a fractional-coordinate branch for triclinic
cells) and the binning replicates ``np.floor_divide`` exactly,
including its ``fmod`` based edge-case handling. Two shortcuts skip
work only where the result is provably unchanged: the imaging
division is skipped when ``|d| <= L/2`` (the rounded quotient is
guaranteed to be 0 then) and pairs are rejected before the square
root when the squared distance lies safely beyond the last bin edge.

The extension must be compiled without floating point contraction
(``-ffp-contract=off``) so that no multiply-add sequence is fused
into an FMA, which would round differently than the separate numpy
operations.

A pure Python/numpy fallback with the identical signature lives in
:py:mod:`PQAnalysis.analysis.rdf._rdf_kernel_py`.
"""

import numpy as np

cimport numpy as np

from libc.float cimport DBL_MAX
from libc.math cimport fabs, floor, fmod, rint, sqrt


cdef inline double _floor_divide(double a, double b) noexcept nogil:
    """
    Replicates ``np.floor_divide`` for float64 scalars (b != 0).

    This is the exact algorithm of numpy's ``npy_floor_divide``: the
    quotient is derived from the exact ``fmod`` remainder, adjusted
    to floor semantics for opposite signs and snapped to the nearest
    integer.
    """
    cdef double mod = fmod(a, b)
    cdef double div = (a - mod) / b
    cdef double floordiv

    if mod != 0.0:
        if (b < 0.0) != (mod < 0.0):
            mod += b
            div -= 1.0

    floordiv = floor(div)

    if div - floordiv > 0.5:
        floordiv += 1.0

    return floordiv


def rdf_frame_histogram(
    const np.float32_t[:, ::1] values,
    const np.int64_t[::1] reference_indices,
    const np.int64_t[::1] target_indices,
    const np.float64_t[::1] box_lengths,
    const np.float64_t[:, ::1] box,
    const np.float64_t[:, ::1] inv_box,
    long long is_orthorhombic,
    double r_min,
    double delta_r,
    long long n_bins,
    np.int64_t[::1] hist,
):
    """
    Accumulates the distance histogram of one trajectory frame.

    For every reference atom the minimum image distances to all
    target atoms (excluding the reference atom itself) are computed
    and scattered into the histogram bins
    ``floor((distance - r_min) / delta_r)``. Distances outside of
    ``[0, n_bins)`` bins are discarded. The histogram accumulator is
    updated in place.

    Parameters
    ----------
    values : np.float32 array of shape (n_atoms, 3), C-contiguous
        The raw frame values (positions) of all atoms of the frame.
    reference_indices : np.int64 array of shape (n_ref,)
        The indices of the reference atoms.
    target_indices : np.int64 array of shape (n_tgt,)
        The indices of the target atoms.
    box_lengths : np.float64 array of shape (3,)
        The box lengths of the current frame. Only used for
        orthorhombic (and vacuum) cells.
    box : np.float64 array of shape (3, 3), C-contiguous
        The box matrix of the current frame. Only used for
        non-orthorhombic cells.
    inv_box : np.float64 array of shape (3, 3), C-contiguous
        The inverse box matrix of the current frame. Only used for
        non-orthorhombic cells.
    is_orthorhombic : int
        Whether all box angles of the cell are exactly 90 degrees
        (vacuum cells included), selecting the box-length imaging
        branch of :py:meth:`PQAnalysis.core.cell.cell.Cell.image`.
    r_min : float
        The minimum (starting) radius of the RDF analysis.
    delta_r : float
        The spacing between the histogram bins.
    n_bins : int
        The number of histogram bins.
    hist : np.int64 array of shape (n_bins,)
        The histogram accumulator; updated in place.
    """

    cdef Py_ssize_t n_ref = reference_indices.shape[0]
    cdef Py_ssize_t n_tgt = target_indices.shape[0]
    cdef Py_ssize_t i, j, row
    cdef np.int64_t ref_index
    cdef float ref_x, ref_y, ref_z
    cdef float delta32_x, delta32_y, delta32_z
    cdef double dx, dy, dz
    cdef double f0, f1, f2
    cdef double dist_sq, dist, shifted, quotient, bin_index, fraction

    cdef double length_x = box_lengths[0]
    cdef double length_y = box_lengths[1]
    cdef double length_z = box_lengths[2]
    cdef double half_x = 0.5 * length_x
    cdef double half_y = 0.5 * length_y
    cdef double half_z = 0.5 * length_z

    cdef double b00 = box[0, 0], b01 = box[0, 1], b02 = box[0, 2]
    cdef double b10 = box[1, 0], b11 = box[1, 1], b12 = box[1, 2]
    cdef double b20 = box[2, 0], b21 = box[2, 1], b22 = box[2, 2]
    cdef double i00 = inv_box[0, 0], i01 = inv_box[0, 1], i02 = inv_box[0, 2]
    cdef double i10 = inv_box[1, 0], i11 = inv_box[1, 1], i12 = inv_box[1, 2]
    cdef double i20 = inv_box[2, 0], i21 = inv_box[2, 1], i22 = inv_box[2, 2]

    cdef double n_bins_d = <double> n_bins

    # the |d| <= L/2 imaging shortcut requires finite box lengths
    # (rint of the quotient is then guaranteed to be exactly 0)
    cdef bint use_half_box_shortcut = (
        is_orthorhombic != 0
        and length_x <= DBL_MAX and length_x == length_x
        and length_y <= DBL_MAX and length_y == length_y
        and length_z <= DBL_MAX and length_z == length_z
    )

    # conservative squared cutoff: any pair with a squared distance
    # above it certainly falls beyond the last bin edge (the slack
    # generously covers all rounding of the exact path); pairs below
    # the cutoff always take the exact path
    cdef double reject = (r_min + delta_r * n_bins_d) * (1.0 + 1e-9) + 1e-9
    cdef double reject_sq = reject * reject

    for i in range(n_ref):
        ref_index = reference_indices[i]
        row = <Py_ssize_t> ref_index
        ref_x = values[row, 0]
        ref_y = values[row, 1]
        ref_z = values[row, 2]

        for j in range(n_tgt):
            if target_indices[j] == ref_index:
                continue

            row = <Py_ssize_t> target_indices[j]

            # float32 pair displacement (exact float32 arithmetic,
            # as in the original loop), widened to float64
            delta32_x = values[row, 0] - ref_x
            delta32_y = values[row, 1] - ref_y
            delta32_z = values[row, 2] - ref_z

            dx = <double> delta32_x
            dy = <double> delta32_y
            dz = <double> delta32_z

            if is_orthorhombic:
                # d - L * rint(d / L); for |d| <= L/2 the rounded
                # quotient is exactly 0 and d is unchanged
                if use_half_box_shortcut:
                    if fabs(dx) > half_x:
                        dx = dx - length_x * rint(dx / length_x)
                    if fabs(dy) > half_y:
                        dy = dy - length_y * rint(dy / length_y)
                    if fabs(dz) > half_z:
                        dz = dz - length_z * rint(dz / length_z)
                else:
                    dx = dx - length_x * rint(dx / length_x)
                    dy = dy - length_y * rint(dy / length_y)
                    dz = dz - length_z * rint(dz / length_z)
            else:
                # fractional = d @ inv_box.T
                f0 = dx * i00 + dy * i01 + dz * i02
                f1 = dx * i10 + dy * i11 + dz * i12
                f2 = dx * i20 + dy * i21 + dz * i22

                f0 -= rint(f0)
                f1 -= rint(f1)
                f2 -= rint(f2)

                # d = fractional @ box.T
                dx = f0 * b00 + f1 * b01 + f2 * b02
                dy = f0 * b10 + f1 * b11 + f2 * b12
                dz = f0 * b20 + f1 * b21 + f2 * b22

            dist_sq = (dx * dx + dy * dy) + dz * dz

            if dist_sq > reject_sq:
                continue

            dist = sqrt(dist_sq)

            shifted = dist - r_min
            quotient = shifted / delta_r
            bin_index = floor(quotient)
            fraction = quotient - bin_index

            # the plain floor of the quotient equals np.floor_divide
            # whenever the quotient is safely away from an integer
            # (|quotient| < 2^33 bounds the rounding error of the
            # division by 2^-20); otherwise take the exact replica
            if not (
                fraction > 1e-6
                and fraction < 0.999999
                and fabs(quotient) < 8.589934592e9
            ):
                bin_index = _floor_divide(shifted, delta_r)

            if bin_index >= 0.0 and bin_index < n_bins_d:
                hist[<Py_ssize_t> bin_index] += 1
