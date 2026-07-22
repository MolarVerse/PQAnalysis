"""
Pure Python/numpy fallback for the VACF accumulation kernels.

This module mirrors the API of the Cython extension
:py:mod:`PQAnalysis.analysis.vacf._vacf_kernel` and is used when the
extension is not available. It is the numeric reference
implementation: :py:func:`accumulate_frame` performs the per-origin
dot products with ``np.einsum``/``np.sum`` exactly like the original
pure-numpy VACF hot loop, while the Cython kernel uses a fixed 4-way
unrolled float64 summation - the two can differ by floating point
rounding on the order of the machine epsilon. All other functions
are bitwise identical between the two implementations.
"""

import numpy as np



def accumulate_frame(
    corr,
    origin_vel,
    origin_norm,
    origin_frame,
    n_active,
    vel,
    frame_number,
    spawn,
    window_size,
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
    if spawn:
        norm = np.sum(vel * vel)

        if norm == 0.0:
            return -1

        origin_vel[n_active] = vel
        origin_norm[n_active] = norm
        origin_frame[n_active] = frame_number
        n_active += 1

    if n_active == 0:
        return 0

    scalars = np.einsum(
        "omd,md->o",
        origin_vel[:n_active],
        vel,
    )
    lags = frame_number - origin_frame[:n_active]
    corr[lags] += scalars / origin_norm[:n_active]

    if lags[0] == window_size:
        # retire the oldest origin
        origin_vel[:n_active - 1] = origin_vel[1:n_active]
        origin_norm[:n_active - 1] = origin_norm[1:n_active]
        origin_frame[:n_active - 1] = origin_frame[1:n_active]
        n_active -= 1

    return n_active



def weight_frame(values, indices, charges):
    """
    Builds the (charge weighted) selected float64 velocities of a
    frame from the raw float32 values of the trajectory reader.

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
    vel = np.asarray(values, dtype=np.float64)[indices]

    if charges is not None:
        vel = vel * charges[:, None]

    return vel



def parse_charge_lines(lines, n_atoms):
    """
    Parses the float64 charge values of 'name charge' body lines.

    Every line must consist of exactly two whitespace separated
    tokens; the second token is parsed as a float64.

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
    charges = np.empty(n_atoms, dtype=np.float64)

    for i in range(n_atoms):
        fields = lines[i].split()

        if len(fields) != 2:
            raise ValueError("Could not parse line")

        charges[i] = float(fields[1])

    return charges
