# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: cdivision=True
"""
Cython kernel for the ADF angle-histogram hot loop.

The kernel accumulates the minimum image angle histogram of one
trajectory frame: for every reference (center) atom ``i`` the
``j-i-k`` angle at the apex ``i`` is computed for every ligand-1 atom
``j`` and ligand-2 atom ``k`` (excluding ``j == i``, ``k == i`` and
``k == j``) and scattered into the bins ``floor(angle / delta_angle)``.
Ordered ligand pairs are counted, exactly as the legacy tool did.

The numeric semantics replicate the numpy fallback bit for bit: the
pair displacements are computed in float32 (the dtype of the raw frame
values) and widened to float64, the minimum image convention uses the
exact op order of :py:meth:`PQAnalysis.core.cell.cell.Cell.image` (an
orthorhombic box-length branch with round-half-even ``rint`` matching
``np.round``, and a fractional-coordinate branch for triclinic cells),
the ``i-j`` and ``i-k`` norms and the ``v_ij . v_ik`` dot product are
reduced in the fixed ``x + y + z`` order and the cosine is formed as
``dot / r_ij / r_ik`` (the legacy division order). The cosine is
clamped to ``[-1, 1]`` before ``acos`` to avoid NaNs from rounding of
nearly (anti-)parallel vectors (the legacy tool omitted this clamp, a
latent NaN bug that is intentionally not reproduced). The angle is
converted to degrees with the exact ``np.degrees`` factor and binned
with a plain ``floor`` (as the legacy tool did, unlike the RDF
distance binning which uses ``np.floor_divide``).

The extension must be compiled without floating point contraction
(``-ffp-contract=off``) so that no multiply-add sequence is fused into
an FMA, which would round differently than the separate numpy
operations.

A pure Python/numpy fallback with the identical signature lives in
:py:mod:`PQAnalysis.analysis.adf._adf_kernel_py`.
"""

import numpy as np

cimport numpy as np

from libc.float cimport DBL_MAX
from libc.math cimport acos, fabs, floor, rint, sqrt

#: The exact float64 radian-to-degree factor of ``np.degrees`` (equal
#: to ``np.degrees(1.0)``), shared with the numpy fallback.
cdef double RAD2DEG = 57.29577951308232


def adf_frame_histogram(
    const np.float32_t[:, ::1] values,
    const np.int64_t[::1] reference_indices,
    const np.int64_t[::1] target_indices,
    const np.int64_t[::1] target_indices_2,
    const np.float64_t[::1] box_lengths,
    const np.float64_t[:, ::1] box,
    const np.float64_t[:, ::1] inv_box,
    long long is_orthorhombic,
    long long gate_1,
    double r_min_1,
    double r_max_1,
    long long gate_2,
    double r_min_2,
    double r_max_2,
    double delta_angle,
    long long n_bins,
    np.int64_t[::1] hist,
):
    """
    Accumulates the angle histogram of one trajectory frame.

    For every reference (center) atom the ``j-i-k`` angles spanned by
    the ligand-1 atoms ``j`` and the ligand-2 atoms ``k`` are computed
    (excluding ``j == i``, ``k == i`` and ``k == j``) and scattered
    into the angle bins ``floor(angle / delta_angle)``. Ordered ligand
    pairs are counted. Angles outside of ``[0, n_bins)`` bins (an exact
    ``180`` degree collinear triplet lands in bin ``n_bins``) are
    discarded. The histogram accumulator is updated in place.

    Parameters
    ----------
    values : np.float32 array of shape (n_atoms, 3), C-contiguous
        The raw frame values (positions) of all atoms of the frame.
    reference_indices : np.int64 array of shape (n_ref,)
        The indices of the reference (center ``i``) atoms.
    target_indices : np.int64 array of shape (n_tgt,)
        The indices of the ligand-1 (``j``) atoms.
    target_indices_2 : np.int64 array of shape (n_tgt2,)
        The indices of the ligand-2 (``k``) atoms.
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
    gate_1 : int
        Whether the radial gate on the ``i-j`` distance is active.
    r_min_1 : float
        The lower ``i-j`` gate radius (inclusive). Only used if
        ``gate_1``.
    r_max_1 : float
        The upper ``i-j`` gate radius (exclusive). Only used if
        ``gate_1``.
    gate_2 : int
        Whether the radial gate on the ``i-k`` distance is active.
    r_min_2 : float
        The lower ``i-k`` gate radius (inclusive). Only used if
        ``gate_2``.
    r_max_2 : float
        The upper ``i-k`` gate radius (exclusive). Only used if
        ``gate_2``.
    delta_angle : float
        The width (in degrees) of the histogram angle bins.
    n_bins : int
        The number of histogram angle bins.
    hist : np.int64 array of shape (n_bins,)
        The histogram accumulator; updated in place.
    """

    cdef Py_ssize_t n_ref = reference_indices.shape[0]
    cdef Py_ssize_t n_tgt = target_indices.shape[0]
    cdef Py_ssize_t n_tgt2 = target_indices_2.shape[0]
    cdef Py_ssize_t ii, jj, kk, row
    cdef np.int64_t ref_index, j_index, k_index
    cdef float ref_x, ref_y, ref_z
    cdef float delta32_x, delta32_y, delta32_z
    cdef double ij_x, ij_y, ij_z, ik_x, ik_y, ik_z
    cdef double f0, f1, f2
    cdef double r_ij, r_ik, dot, cosine, angle, bin_index

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

    for ii in range(n_ref):
        ref_index = reference_indices[ii]
        row = <Py_ssize_t> ref_index
        ref_x = values[row, 0]
        ref_y = values[row, 1]
        ref_z = values[row, 2]

        for jj in range(n_tgt):
            j_index = target_indices[jj]
            if j_index == ref_index:
                continue

            row = <Py_ssize_t> j_index

            # float32 pair displacement widened to float64
            delta32_x = values[row, 0] - ref_x
            delta32_y = values[row, 1] - ref_y
            delta32_z = values[row, 2] - ref_z

            ij_x = <double> delta32_x
            ij_y = <double> delta32_y
            ij_z = <double> delta32_z

            if is_orthorhombic:
                if use_half_box_shortcut:
                    if fabs(ij_x) > half_x:
                        ij_x = ij_x - length_x * rint(ij_x / length_x)
                    if fabs(ij_y) > half_y:
                        ij_y = ij_y - length_y * rint(ij_y / length_y)
                    if fabs(ij_z) > half_z:
                        ij_z = ij_z - length_z * rint(ij_z / length_z)
                else:
                    ij_x = ij_x - length_x * rint(ij_x / length_x)
                    ij_y = ij_y - length_y * rint(ij_y / length_y)
                    ij_z = ij_z - length_z * rint(ij_z / length_z)
            else:
                f0 = ij_x * i00 + ij_y * i01 + ij_z * i02
                f1 = ij_x * i10 + ij_y * i11 + ij_z * i12
                f2 = ij_x * i20 + ij_y * i21 + ij_z * i22

                f0 -= rint(f0)
                f1 -= rint(f1)
                f2 -= rint(f2)

                ij_x = f0 * b00 + f1 * b01 + f2 * b02
                ij_y = f0 * b10 + f1 * b11 + f2 * b12
                ij_z = f0 * b20 + f1 * b21 + f2 * b22

            r_ij = sqrt((ij_x * ij_x + ij_y * ij_y) + ij_z * ij_z)

            if gate_1 and (r_ij < r_min_1 or r_ij >= r_max_1):
                continue

            for kk in range(n_tgt2):
                k_index = target_indices_2[kk]
                if k_index == ref_index or k_index == j_index:
                    continue

                row = <Py_ssize_t> k_index

                delta32_x = values[row, 0] - ref_x
                delta32_y = values[row, 1] - ref_y
                delta32_z = values[row, 2] - ref_z

                ik_x = <double> delta32_x
                ik_y = <double> delta32_y
                ik_z = <double> delta32_z

                if is_orthorhombic:
                    if use_half_box_shortcut:
                        if fabs(ik_x) > half_x:
                            ik_x = ik_x - length_x * rint(ik_x / length_x)
                        if fabs(ik_y) > half_y:
                            ik_y = ik_y - length_y * rint(ik_y / length_y)
                        if fabs(ik_z) > half_z:
                            ik_z = ik_z - length_z * rint(ik_z / length_z)
                    else:
                        ik_x = ik_x - length_x * rint(ik_x / length_x)
                        ik_y = ik_y - length_y * rint(ik_y / length_y)
                        ik_z = ik_z - length_z * rint(ik_z / length_z)
                else:
                    f0 = ik_x * i00 + ik_y * i01 + ik_z * i02
                    f1 = ik_x * i10 + ik_y * i11 + ik_z * i12
                    f2 = ik_x * i20 + ik_y * i21 + ik_z * i22

                    f0 -= rint(f0)
                    f1 -= rint(f1)
                    f2 -= rint(f2)

                    ik_x = f0 * b00 + f1 * b01 + f2 * b02
                    ik_y = f0 * b10 + f1 * b11 + f2 * b12
                    ik_z = f0 * b20 + f1 * b21 + f2 * b22

                r_ik = sqrt((ik_x * ik_x + ik_y * ik_y) + ik_z * ik_z)

                if gate_2 and (r_ik < r_min_2 or r_ik >= r_max_2):
                    continue

                dot = ij_x * ik_x + ij_y * ik_y + ij_z * ik_z

                cosine = dot / r_ij / r_ik

                if cosine > 1.0:
                    cosine = 1.0
                elif cosine < -1.0:
                    cosine = -1.0

                angle = acos(cosine) * RAD2DEG

                bin_index = floor(angle / delta_angle)

                if bin_index >= 0.0 and bin_index < n_bins_d:
                    hist[<Py_ssize_t> bin_index] += 1
