"""
Pure Python/numpy fallback for the MSD accumulation kernel.

This module mirrors the API of the Cython extension
:py:mod:`PQAnalysis.analysis.msd._msd_kernel` and is used when the
extension is not available. It implements the frame update with the
exact numpy operations of the original pure numpy MSD hot loop
(``fractional = d @ B^-1.T`` and ``shift += -np.rint(fractional) @ B.T``
with round-half-even ``np.rint``), so its results are bit-identical to
that implementation.
"""

import numpy as np


def msd_frame_update(
    values,
    indices,
    box,
    inv_box,
    is_vacuum,
    pos,
    prev_pos,
    shift,
    unwrapped,
    origins,
    msd,
    state,
    counter,
    gap,
    window,
    n_start,
    stop_frame,
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

    # gather the selected positions (exact float32 -> float64)
    pos[:] = values[indices]

    if counter > 1 and not is_vacuum:
        # unwrap: fold the per-frame displacements back into the
        # minimum image convention and accumulate the change
        fractional = (pos - prev_pos) @ inv_box.T
        shift += -np.rint(fractional) @ box.T

    prev_pos[:] = pos
    np.add(pos, shift, out=unwrapped)

    if counter < n_start:
        return

    n_active = int(state[0])
    last = int(state[1])

    if counter % gap == 0:
        if n_active != origins.shape[0] and counter <= stop_frame:
            # starting stage - add new origin
            if n_active == 0:
                last = counter

            origins[n_active] = unwrapped
            n_active += 1

        elif n_active > 0 and last + window == counter:
            # oldest origin reached the full window:
            # accumulate its final lag term
            disp = unwrapped - origins[0]
            msd[window] += np.einsum('ax,ax->x', disp, disp)

            origins[:n_active - 1] = origins[1:n_active]

            if counter > stop_frame:
                # stopping stage - drop without replacement
                n_active -= 1
            else:
                # running stage - replace by a new origin
                origins[n_active - 1] = unwrapped

            last += gap

        state[0] = n_active
        state[1] = last

    if n_active > 0:
        # accumulate the squared displacements of all active origins
        lags = counter - last - gap * np.arange(n_active)

        disp = unwrapped - origins[:n_active]

        msd[lags] += np.einsum('oax,oax->ox', disp, disp)
