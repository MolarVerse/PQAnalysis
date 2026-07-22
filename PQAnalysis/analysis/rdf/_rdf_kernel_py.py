"""
Pure Python/numpy fallback for the RDF distance-histogram kernel.

This module mirrors the API of the Cython extension
:py:mod:`PQAnalysis.analysis.rdf._rdf_kernel` and is used when the
extension is not available. It implements the per-frame histogram
update as a refactoring of the original RDF hot loop with the exact
numpy operations of that loop: the pair displacements are computed in
float32 (the dtype of the raw frame values), imaged into the unit cell
with the operations of
:py:meth:`PQAnalysis.core.cell.cell.Cell.image` (which promote the
displacements to float64), reduced with ``np.linalg.norm`` and binned
with ``np.floor_divide``/``np.bincount`` exactly as in
:py:meth:`PQAnalysis.analysis.rdf.rdf.RDF._add_to_bins`, so its
results are bit-identical to that implementation.
"""

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

    for reference_index in reference_indices:
        selected_target_indices = target_indices[
            target_indices != reference_index]

        # float32 pair displacements, exactly as in the original
        # loop (frame.pos is float32 for file based trajectories)
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
