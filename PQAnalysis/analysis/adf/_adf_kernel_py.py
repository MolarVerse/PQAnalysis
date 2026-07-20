"""
Pure Python/numpy fallback for the ADF angle-histogram kernel.

This module mirrors the API of the Cython extension
:py:mod:`PQAnalysis.analysis.adf._adf_kernel` and is used when the
extension is not available. It also serves as the per-frame worker of
the original (non fast-path) ADF loop, so that the raw-frame fast path
and the original path produce bit-identical results.

For every reference (center) atom ``i`` the minimum image vectors to
all ligand-1 atoms ``j`` (``v_ij``) and all ligand-2 atoms ``k``
(``v_ik``) are computed, the ``j-i-k`` angle at the apex ``i`` is
evaluated as ``degrees(acos(v_ij . v_ik / (r_ij * r_ik)))`` and
scattered into the angle bins ``floor(angle / delta_angle)``. Optional
radial gates on ``r_ij`` and ``r_ik`` restrict the counted triplets.

The pair displacements are computed in the dtype of the raw frame
values (float32 for file based trajectories) and imaged into the unit
cell with the exact operations of
:py:meth:`PQAnalysis.core.cell.cell.Cell.image` (which promote the
displacements to float64), replicated bit for bit from
:py:mod:`PQAnalysis.analysis.rdf._rdf_kernel_py`. The dot products and
norms are reduced in the fixed ``x + y + z`` order (no BLAS, no fused
multiply-add) so that the numeric result is bit-identical to the
scalar Cython kernel. The acos argument is clamped to ``[-1, 1]``
before ``arccos`` to avoid NaNs from rounding of nearly (anti-)parallel
vectors (the legacy tool omitted this clamp, a latent NaN bug that is
intentionally not reproduced).
"""

import numpy as np

#: The exact float64 radian-to-degree factor of ``np.degrees`` (equal
#: to ``np.degrees(1.0)``); used verbatim by the Cython kernel too.
_RAD2DEG = 57.29577951308232


def _minimum_image(delta, box_lengths, box, inv_box, is_orthorhombic):
    """
    Images the pair displacements into the unit cell.

    Replicates the operations of
    :py:meth:`PQAnalysis.core.cell.cell.Cell.image` exactly (the
    float32 displacements are promoted to float64 by the float64 box
    data), matching :py:mod:`PQAnalysis.analysis.rdf._rdf_kernel_py`.

    Parameters
    ----------
    delta : np.ndarray of shape (n, 3)
        The raw pair displacements (target minus reference position).
    box_lengths : np.float64 array of shape (3,)
        The box lengths of the current frame.
    box : np.float64 array of shape (3, 3)
        The box matrix of the current frame.
    inv_box : np.float64 array of shape (3, 3)
        The inverse box matrix of the current frame.
    is_orthorhombic : int
        Whether all box angles of the cell are exactly 90 degrees.

    Returns
    -------
    np.float64 array of shape (n, 3)
        The minimum image displacements.
    """

    if is_orthorhombic:
        return delta - box_lengths * np.round(delta / box_lengths)

    fractional = delta @ inv_box.T
    fractional -= np.round(fractional)

    return fractional @ box.T


def _norm(vectors):
    """
    Euclidean norm reduced in the fixed ``x + y + z`` order.

    Parameters
    ----------
    vectors : np.float64 array of shape (n, 3)
        The vectors to reduce.

    Returns
    -------
    np.float64 array of shape (n,)
        The norms of the vectors.
    """

    return np.sqrt(
        vectors[:, 0] * vectors[:, 0]
        + vectors[:, 1] * vectors[:, 1]
        + vectors[:, 2] * vectors[:, 2]
    )


def adf_frame_histogram(
    values,
    reference_indices,
    target_indices,
    target_indices_2,
    box_lengths,
    box,
    inv_box,
    is_orthorhombic,
    gate_1,
    r_min_1,
    r_max_1,
    gate_2,
    r_min_2,
    r_max_2,
    delta_angle,
    n_bins,
    hist,
):
    """
    Accumulates the angle histogram of one trajectory frame.

    For every reference (center) atom the ``j-i-k`` angles spanned by
    the ligand-1 atoms ``j`` and the ligand-2 atoms ``k`` are computed
    (excluding ``j == i``, ``k == i`` and ``k == j``) and scattered
    into the angle bins ``floor(angle / delta_angle)``. Ordered ligand
    pairs are counted, exactly as the legacy tool did. Angles that fall
    outside of ``[0, n_bins)`` bins (an exact ``180`` degree collinear
    triplet lands in bin ``n_bins``) are discarded. The histogram
    accumulator is updated in place.

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
        (vacuum cells included).
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

    for reference_index in reference_indices:
        reference_position = values[reference_index]

        # ligand-1 atoms j (exclude j == i)
        j_indices = target_indices[target_indices != reference_index]
        v_ij = _minimum_image(
            values[j_indices] - reference_position,
            box_lengths,
            box,
            inv_box,
            is_orthorhombic,
        )
        r_ij = _norm(v_ij)

        if gate_1:
            gate_1_mask = (r_ij >= r_min_1) & (r_ij < r_max_1)
            j_indices = j_indices[gate_1_mask]
            v_ij = v_ij[gate_1_mask]
            r_ij = r_ij[gate_1_mask]

        # ligand-2 atoms k (exclude k == i)
        k_indices = target_indices_2[target_indices_2 != reference_index]
        v_ik = _minimum_image(
            values[k_indices] - reference_position,
            box_lengths,
            box,
            inv_box,
            is_orthorhombic,
        )
        r_ik = _norm(v_ik)

        if gate_2:
            gate_2_mask = (r_ik >= r_min_2) & (r_ik < r_max_2)
            k_indices = k_indices[gate_2_mask]
            v_ik = v_ik[gate_2_mask]
            r_ik = r_ik[gate_2_mask]

        if j_indices.shape[0] == 0 or k_indices.shape[0] == 0:
            continue

        # dot products v_ij . v_ik reduced in the fixed x + y + z order
        dot = (
            v_ij[:, 0][:, None] * v_ik[:, 0][None, :]
            + v_ij[:, 1][:, None] * v_ik[:, 1][None, :]
            + v_ij[:, 2][:, None] * v_ik[:, 2][None, :]
        )

        # cosine of the j-i-k angle: dot / r_ij / r_ik (legacy op order)
        cosine = dot / r_ij[:, None] / r_ik[None, :]

        # clamp to a valid acos domain (documented, avoids NaNs)
        cosine = np.clip(cosine, -1.0, 1.0)

        angle = np.arccos(cosine) * _RAD2DEG

        bin_indices = np.floor(angle / delta_angle).astype(np.int64)

        # ordered ligand pairs, but never pair an atom with itself
        keep = j_indices[:, None] != k_indices[None, :]
        bin_indices = bin_indices[keep]

        bin_indices = bin_indices[(bin_indices >= 0) & (bin_indices < n_bins)]

        hist += np.bincount(bin_indices, minlength=n_bins)
