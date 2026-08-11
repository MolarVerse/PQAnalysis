"""
Pure Python/numpy fallback for the RDF distance-histogram kernel.

This module mirrors the API of the Cython extension
:py:mod:`PQAnalysis.analysis.rdf._rdf_kernel` and is used when the
extension is not available. It implements the per-frame histogram
update as a refactoring of the original RDF hot loop with the exact
numpy operations of that loop: pair displacements retain the dtype of
the raw frame values and are imaged into the unit cell
with the operations of
:py:meth:`PQAnalysis.core.cell.cell.Cell.image` (which promote the
float32 displacements to float64), reduced with ``np.linalg.norm`` and
binned with ``np.floor_divide``/``np.bincount`` exactly as in
:py:meth:`PQAnalysis.analysis.rdf.rdf.RDF._add_to_bins`, so its
results are bit-identical to that implementation.

The separate :func:`legacy_rdf_frame_histogram` and
:func:`legacy_rdf_batch_histogram` functions preserve the operation order of
the legacy ``thh_tools/RDF`` C loop for file-backed orthorhombic trajectories.
"""

import math

import numpy as np


def rdf_frame_histogram(
    values,
    reference_indices,
    target_indices,
    box_lengths,
    box,
    inv_box,
    is_orthorhombic,
    r_min,
    delta_r,
    n_bins,
    hist,
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
    values : np.float32 | np.float64 array of shape (n_atoms, 3), C-contiguous
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

    for reference_index in reference_indices:
        selected_target_indices = target_indices[
            target_indices != reference_index]

        # Pair displacements retain the source-array precision.
        delta = values[selected_target_indices] - values[reference_index]

        # minimum image convention with the exact operations of
        # Cell.image (the float32 displacements are promoted to
        # float64 by the float64 box data)
        if is_orthorhombic:
            delta = delta - box_lengths * np.round(delta / box_lengths)
        else:
            fractional = delta @ inv_box.T
            fractional -= np.round(fractional)
            delta = fractional @ box.T

        distances = np.linalg.norm(delta, axis=-1)

        # binning, exactly as in RDF._add_to_bins
        bin_indices = np.floor_divide(
            distances - r_min, delta_r
        ).astype(int)

        bin_indices = bin_indices[
            (bin_indices < n_bins) & (bin_indices >= 0)]

        hist += np.bincount(bin_indices, minlength=n_bins)


def _round_away_from_zero(value):
    """Replicates C99 ``round`` for a finite float64 value."""
    if value >= 0.0:
        return math.floor(value + 0.5)

    return math.ceil(value - 0.5)


# pylint: disable=too-many-arguments,too-many-positional-arguments
# pylint: disable=too-many-locals
def legacy_rdf_frame_histogram(
    values,
    reference_indices,
    target_indices,
    box_lengths,
    delta_r,
    n_bins,
    hist,
):
    """Accumulates one frame with the legacy RDF operation order."""
    length_x = float(box_lengths[0])
    length_y = float(box_lengths[1])
    length_z = float(box_lengths[2])
    cutoff = min(length_x, length_y, length_z) / 2.0

    for reference_index in reference_indices:
        ref_index = int(reference_index)
        ref_x = float(values[ref_index, 0])
        ref_y = float(values[ref_index, 1])
        ref_z = float(values[ref_index, 2])

        for target_index_value in target_indices:
            target_index = int(target_index_value)

            if target_index == ref_index:
                continue

            target_x = float(values[target_index, 0])
            target_y = float(values[target_index, 1])
            target_z = float(values[target_index, 2])

            image_x = target_x + length_x * _round_away_from_zero(
                (ref_x - target_x) / length_x
            )
            image_y = target_y + length_y * _round_away_from_zero(
                (ref_y - target_y) / length_y
            )
            image_z = target_z + length_z * _round_away_from_zero(
                (ref_z - target_z) / length_z
            )

            dx = ref_x - image_x
            dy = ref_y - image_y
            dz = ref_z - image_z
            distance = math.sqrt((dx * dx + dy * dy) + dz * dz)

            if distance <= cutoff:
                bin_index = math.floor(distance / delta_r)

                # The legacy allocation has one hidden cutoff bin that
                # is never written. Ignore it instead of writing beyond
                # the returned histogram.
                if 0 <= bin_index < n_bins:
                    hist[bin_index] += 1


def legacy_rdf_batch_histogram(
    values,
    reference_indices,
    target_indices,
    box_lengths,
    delta_r,
    n_bins,
    hist,
):
    """Accumulates a frame batch with the legacy RDF operation order."""
    if len(box_lengths) != len(values):
        raise ValueError(
            "The number of RDF boxes must match the number of frames."
        )

    for frame_values, frame_box_lengths in zip(values, box_lengths):
        legacy_rdf_frame_histogram(
            frame_values,
            reference_indices,
            target_indices,
            frame_box_lengths,
            delta_r,
            n_bins,
            hist,
        )
