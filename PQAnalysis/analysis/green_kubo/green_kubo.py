"""
A module containing the GreenKubo class. The GreenKubo class is used
to calculate transport coefficients from the Green-Kubo relations of a
target selection of a given velocity trajectory. Currently the
self-diffusion coefficient is implemented, which is obtained from the
time integral of the un-normalized velocity auto-correlation function

    D = (1/3) integral_0^inf Cvv(t) dt,

    Cvv(t) = (1/N) sum_{i in selection} < v_i(0) . v_i(t) >,

where Cvv is the velocity auto-correlation function in absolute units
(NOT the C(0) = 1 normalized vacf), averaged over the N selected atoms
and over time origins with the full three-component dot product. The
factor 1/3 turns the isotropic integral into the scalar diffusion
coefficient.

The module is intentionally structured so that further Green-Kubo
transport coefficients (electrical conductivity from the charge-flux
auto-correlation, shear viscosity from the stress auto-correlation)
can be added as additional front-ends on top of the same
auto-correlation and running-integral machinery.
"""

import itertools
import logging

# 3rd party imports
import numpy as np
from beartype.typing import Generator, Tuple
from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis import config
from PQAnalysis.types import (
    Np1DNumberArray,
    Np2DNumberArray,
    NpnDNumberArray,
    PositiveInt,
    PositiveReal,
)
from PQAnalysis.traj import Trajectory, TrajectoryFormat
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import RawTrajectoryReader, TrajectoryReader
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import GreenKuboError

#: float: Conversion factor of a velocity auto-correlation running
#: integral from (Angstrom / s)^2 * ps to m^2 / s. The velocity
#: trajectory (.vel) stores velocities in Angstrom / s (thermal atomic
#: speeds are of the order of 1e12 - 1e13 Angstrom / s), the lag-time
#: axis is built with the time_step in ps, so the running integral
#: G(t) = integral Cvv dt has the units (Angstrom / s)^2 * ps. The
#: derivation of the factor is:
#:
#:   (Angstrom / s)^2 * ps = (1e-10 m / s)^2 * (1e-12 s)
#:                         = 1e-20 m^2 s^-2 * 1e-12 s
#:                         = 1e-32 m^2 / s.
#:
#: With this factor D = (1/3) * ANGSTROM2_PER_S2_PS_TO_M2_PER_S * G is
#: obtained in m^2 / s, matching the m^2 / s convention of the MSD
#: (Einstein) module (ANGSTROM2_PER_PS_TO_M2_PER_S = 1e-8): a velocity
#: series v = dx/dt of an Angstrom/ps position trajectory expressed in
#: Angstrom/s (i.e. multiplied by 1e12) reproduces the Einstein
#: diffusion coefficient of that trajectory exactly.
ANGSTROM2_PER_S2_PS_TO_M2_PER_S = 1.0e-32

#: float: Conversion factor from m^2 / s to cm^2 / s (the common
#: diffusion unit): 1 m^2 / s = 1e4 cm^2 / s.
M2_PER_S_TO_CM2_PER_S = 1.0e4



def velocity_acf_fft(
    velocities: NpnDNumberArray,
    window_size: PositiveInt,
) -> Np1DNumberArray:
    """
    Calculates the un-normalized velocity auto-correlation function
    with the Wiener-Khinchin (FFT) estimator.

    For the ``(n_frames, n_sel, 3)`` velocity array the aggregate
    auto-correlation ``S(lag) = sum_t sum_j v_j(t) . v_j(t + lag)`` is
    computed per atom and component via a real FFT, each lag is divided
    by its number of time origins ``n_frames - lag`` and finally
    divided by the number of selected atoms. Every frame serves as a
    time origin. In contrast to the vacf module the result is NOT
    normalized by its lag-zero value and is kept in the native
    (Angstrom / s)^2 units.

    Parameters
    ----------
    velocities : NpnDNumberArray
        The selected velocities of shape ``(n_frames, n_sel, 3)``.
    window_size : PositiveInt
        The correlation window length in frames. The correlation
        function is returned for the lags ``0..window_size``.

    Returns
    -------
    Np1DNumberArray
        The un-normalized velocity auto-correlation function
        ``Cvv[0..window_size]`` in (Angstrom / s)^2.

    Raises
    ------
    GreenKuboError
        If the aggregate mean squared velocity of the trajectory is
        zero, so that the diffusion coefficient is not defined.
    """
    n_frames = velocities.shape[0]
    n_sel = velocities.shape[1]

    n_fft = 2 * n_frames
    spectrum = np.fft.rfft(velocities, n=n_fft, axis=0)
    autocorr = np.fft.irfft(
        spectrum * np.conj(spectrum),
        n=n_fft,
        axis=0,
    )[:window_size + 1]

    raw = np.sum(autocorr.real, axis=(1, 2))
    counts = n_frames - np.arange(window_size + 1)
    cvv = raw / counts / n_sel

    if cvv[0] == 0.0:
        raise GreenKuboError(
            "The aggregate mean squared velocity of the trajectory is "
            "zero. The Green-Kubo diffusion coefficient is not defined."
        )

    return cvv



def velocity_acf_direct(
    velocities: NpnDNumberArray,
    window_size: PositiveInt,
    gap: PositiveInt,
) -> Np1DNumberArray:
    """
    Calculates the un-normalized velocity auto-correlation function
    with a direct sliding-time-origin estimator.

    A new time origin is spawned every ``gap`` frames. For every origin
    the three-component dot products ``v_j(t0) . v_j(t0 + lag)`` are
    accumulated over all selected atoms ``j`` into the lag bins that
    still fit into the trajectory. Every lag is divided by its actual
    number of contributing origins and finally by the number of
    selected atoms. For ``gap == 1`` this reproduces the FFT estimator
    of :py:func:`velocity_acf_fft` up to floating point noise; for
    larger gaps it is the same estimator on a sparser set of origins.

    Parameters
    ----------
    velocities : NpnDNumberArray
        The selected velocities of shape ``(n_frames, n_sel, 3)``.
    window_size : PositiveInt
        The correlation window length in frames. The correlation
        function is returned for the lags ``0..window_size``.
    gap : PositiveInt
        The spacing between two time origins in frames.

    Returns
    -------
    Np1DNumberArray
        The un-normalized velocity auto-correlation function
        ``Cvv[0..window_size]`` in (Angstrom / s)^2.

    Raises
    ------
    GreenKuboError
        If the aggregate mean squared velocity of the trajectory is
        zero, so that the diffusion coefficient is not defined.
    """
    n_frames = velocities.shape[0]
    n_sel = velocities.shape[1]

    cvv = np.zeros(window_size + 1, dtype=np.float64)
    counts = np.zeros(window_size + 1, dtype=np.float64)

    for t0 in range(0, n_frames, gap):
        max_lag = min(window_size, n_frames - 1 - t0)
        segment = velocities[t0:t0 + max_lag + 1]
        contribution = np.sum(segment * velocities[t0][None], axis=(1, 2))
        cvv[:max_lag + 1] += contribution
        counts[:max_lag + 1] += 1.0

    cvv = cvv / counts / n_sel

    if cvv[0] == 0.0:
        raise GreenKuboError(
            "The aggregate mean squared velocity of the trajectory is "
            "zero. The Green-Kubo diffusion coefficient is not defined."
        )

    return cvv



def cumulative_trapezoid(values: Np1DNumberArray, dx: PositiveReal) -> Np1DNumberArray:
    """
    Calculates the cumulative trapezoidal integral of a series with a
    leading zero, i.e. ``G[0] = 0`` and
    ``G[k] = G[k-1] + 0.5 * (values[k-1] + values[k]) * dx``.

    Parameters
    ----------
    values : Np1DNumberArray
        The integrand sampled on an equidistant grid.
    dx : PositiveReal
        The spacing of the grid.

    Returns
    -------
    Np1DNumberArray
        The cumulative integral, same length as ``values``.
    """
    increments = 0.5 * (values[1:] + values[:-1]) * dx

    return np.concatenate(([0.0], np.cumsum(increments)))



class GreenKubo:  # pylint: disable=too-many-instance-attributes

    """
    A class for calculating transport coefficients from the Green-Kubo
    relations of a target selection of a velocity trajectory.

    Currently the self-diffusion coefficient is implemented. It is
    obtained from the time integral of the un-normalized velocity
    auto-correlation function

        ``Cvv(lag) = (1 / n_sel) * (1 / n_origins(lag))
        * sum_i sum_j v_j(t0_i) . v_j(t0_i + lag)``

    where the sum over ``j`` runs over all selected atoms and the sum
    over ``i`` over all time origins (every frame is an origin for the
    default ``fft`` estimator). The running integral

        ``G(lag) = cumulative_trapezoid(Cvv, dt = time_step)``

    yields the running diffusion coefficient
    ``D_running(lag) = (1/3) * unit_factor * G(lag)`` in m^2 / s, whose
    plateau is reported as the diffusion coefficient. The full running
    curve is returned so that the user can inspect the plateau, and the
    scalar diffusion coefficient is taken as the mean of ``D_running``
    over a trailing fit window ``[fit_start, fit_stop]`` (given as
    fractions of the correlation window).

    The statistical uncertainty of the diffusion coefficient is
    estimated by *block averaging*, the standard method for Green-Kubo
    transport coefficients: the trajectory frames are split into
    ``n_blocks`` contiguous, roughly-equal, non-overlapping blocks, the
    same velocity-ACF -> running-integral -> plateau pipeline is run
    independently on each block to obtain a per-block diffusion
    coefficient ``D_b`` and the reported uncertainty is the standard
    error of the mean over the blocks
    ``std(D_b, ddof=1) / sqrt(n_blocks)``
    (:py:attr:`diffusion_coefficient_stderr`). This is a genuine
    sampling error, in contrast to the (kept, but honestly named)
    :py:attr:`diffusion_coefficient_plateau_spread`, which is the
    standard deviation of the strongly auto-correlated running integral
    over the fit window and only measures the flatness of the plateau,
    NOT the run-to-run statistical error of ``D``.

    Each block must be long enough to cover the correlation window
    (``block_length > window_size``). If the requested ``n_blocks`` is
    too large for the trajectory and window, it is clamped down to the
    largest number of blocks that still satisfy this constraint and a
    warning is logged; if not even two such blocks fit, a
    :py:class:`~PQAnalysis.analysis.green_kubo.exceptions.GreenKuboError`
    is raised.

    Two estimators are available for the auto-correlation function: the
    default ``fft`` (Wiener-Khinchin) estimator and a ``direct``
    sliding-time-origin estimator; both are un-normalized and kept in
    the native (Angstrom / s)^2 velocity units.

    The class can be initialized with either a Trajectory object or a
    TrajectoryReader. For velocity TrajectoryReaders the raw-frame fast
    path
    (:py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`)
    streams the float32 values without building an AtomicSystem per
    frame, exactly like the vacf module.
    """

    _window_size_default = 1000
    _gap_default = 1
    _method_default = "fft"
    _fit_start_default = 0.5
    _fit_stop_default = 1.0
    _coefficient_default = "diffusion"
    _n_blocks_default = 5

    _supported_methods = ("fft", "direct")
    _supported_coefficients = ("diffusion",)

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-statements
        self,
        traj: Trajectory | TrajectoryReader,
        time_step: PositiveReal,
        window_size: PositiveInt | None = None,
        target_species: SelectionCompatible = None,
        gap: PositiveInt | None = None,
        fit_start: PositiveReal | None = None,
        fit_stop: PositiveReal | None = None,
        method: str | None = None,
        coefficient: str | None = None,
        n_blocks: PositiveInt | None = None,
        use_full_atom_info: bool | None = False,
    ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The velocity trajectory to analyze. If a TrajectoryReader
            is provided, the trajectory is read lazily frame by frame.
            All frames must provide velocities (e.g. .vel files) in
            Angstrom / s.
        time_step : PositiveReal
            The time step between two frames in ps. It is used to build
            the lag-time axis and to integrate the auto-correlation
            function.
        window_size : PositiveInt | None, optional
            The correlation window length in frames. The correlation
            function and its running integral are calculated for the
            lags ``0..window_size``, by default None (1000, capped to
            ``n_frames - 1``).
        target_species : SelectionCompatible, optional
            The target species of the Green-Kubo analysis, by default
            None (all atoms).
        gap : PositiveInt | None, optional
            The spacing between two time origins in frames for the
            ``direct`` estimator, by default None (1). It is ignored by
            the ``fft`` estimator, which uses every frame as an origin.
        fit_start : PositiveReal | None, optional
            The start of the plateau fit window as a fraction of the
            correlation window (0..1), by default None (0.5).
        fit_stop : PositiveReal | None, optional
            The end of the plateau fit window as a fraction of the
            correlation window (0..1), by default None (1.0).
        method : str | None, optional
            The auto-correlation estimator, either ``fft`` (default,
            Wiener-Khinchin) or ``direct`` (sliding time origins), by
            default None (``fft``).
        coefficient : str | None, optional
            The transport coefficient to calculate. Currently only
            ``diffusion`` is supported, by default None (``diffusion``).
        n_blocks : PositiveInt | None, optional
            The number of contiguous, non-overlapping blocks the
            trajectory is split into for the block-averaged statistical
            uncertainty of the diffusion coefficient, by default None
            (5). It must be at least 2. If it is too large for the
            trajectory and the correlation window (each block must be
            longer than ``window_size``) it is clamped down at run time
            with a warning; if not even two blocks fit a GreenKuboError
            is raised.
        use_full_atom_info : bool | None, optional
            Whether to use the full atom information for the target
            selection, by default None (False).

        Raises
        ------
        GreenKuboError
            If the time step is not positive.
        GreenKuboError
            If the method is unknown.
        GreenKuboError
            If the transport coefficient is unknown.
        GreenKuboError
            If the fit window fractions are not within ``0 <= fit_start
            < fit_stop <= 1``.
        GreenKuboError
            If the number of blocks is smaller than 2.
        GreenKuboError
            If the trajectory is empty.
        GreenKuboError
            If the trajectory is too short for the requested window
            (``n_frames <= window_size``).
        GreenKuboError
            If the target selection does not select any atoms.
        """

        ##############
        # dummy init #
        ##############

        self.lag_times = np.array([])
        self.cvv = np.array([])
        self.d_running = np.array([])
        self.diffusion_coefficient = 0.0
        self.diffusion_coefficient_stderr = 0.0
        self.diffusion_coefficient_plateau_spread = 0.0
        self.block_diffusion_coefficients = np.array([])
        self.n_origins = 0

        #####################################################
        # Initialize parameters with default values if None #
        #####################################################

        self.time_step = time_step
        self.gap = gap if gap is not None else self._gap_default
        self.method = (
            method if method is not None else self._method_default
        ).lower()
        self.coefficient = (
            coefficient
            if coefficient is not None else self._coefficient_default
        ).lower()
        self.fit_start = (
            fit_start if fit_start is not None else self._fit_start_default
        )
        self.fit_stop = (
            fit_stop if fit_stop is not None else self._fit_stop_default
        )
        self.n_blocks = (
            n_blocks if n_blocks is not None else self._n_blocks_default
        )

        if use_full_atom_info is None:
            use_full_atom_info = False
        self.use_full_atom_info = use_full_atom_info

        self._check_parameters()

        ############################################
        # Initialize trajectory iterator/generator #
        ############################################

        self._raw_reader = None
        self._frame_generator = None

        if isinstance(traj, TrajectoryReader):
            if traj.traj_format == TrajectoryFormat.VEL:
                # additive fast path: stream the raw float32 values of
                # the velocity trajectory without building an
                # AtomicSystem per frame (bit-identical values)
                self._raw_reader = RawTrajectoryReader(
                    traj.filenames,
                    traj_format=traj.traj_format,
                    md_format=traj.md_format,
                )
                self.n_frames = self._raw_reader.count_frames()
            else:
                self.n_frames = sum(
                    traj.calculate_number_of_frames_per_file()
                )
                self._frame_generator = traj.frame_generator()
        elif len(traj) > 0:
            self.n_frames = len(traj)
            self._frame_generator = iter(traj)
        else:
            self.logger.error(
                "Trajectory cannot be of length 0.",
                exception=GreenKuboError,
            )

        self.window_size = (
            window_size
            if window_size is not None
            else min(self._window_size_default, self.n_frames - 1)
        )

        if self.window_size >= self.n_frames:
            self.logger.error(
                (
                    f"The trajectory contains only {self.n_frames} "
                    "frame(s), but at least window_size + 1 = "
                    f"{self.window_size + 1} frames are needed so that "
                    "every lag has at least one time origin."
                ),
                exception=GreenKuboError,
            )

        if self._raw_reader is not None:
            self._first_frame = self._raw_reader.read_first_frame()
        else:
            self._first_frame = next(self._frame_generator)

        if traj.topology is not None:
            self.topology = traj.topology
        else:
            self.topology = self._first_frame.topology

        ################################
        # Initialize Selection objects #
        ################################

        self.target_species = target_species
        self.target_selection = Selection(target_species)
        self.target_indices = self.target_selection.select(
            self.topology,
            self.use_full_atom_info,
        )

        if len(self.target_indices) == 0:
            self.logger.error(
                "The target selection does not select any atoms.",
                exception=GreenKuboError,
            )

        self._target_indices_intp = np.ascontiguousarray(
            self.target_indices,
            dtype=np.intp,
        )

    def _check_parameters(self) -> None:
        """
        Checks the consistency of the setup parameters.

        This method is called by the __init__ method of the GreenKubo
        class after the scalar parameters have been initialized.

        Raises
        ------
        GreenKuboError
            If the time step is not positive.
        GreenKuboError
            If the method is unknown.
        GreenKuboError
            If the transport coefficient is unknown.
        GreenKuboError
            If the fit window fractions are not within ``0 <= fit_start
            < fit_stop <= 1``.
        GreenKuboError
            If the number of blocks is smaller than 2.
        """
        if self.time_step <= 0.0:
            self.logger.error(
                "The time_step must be a positive real number.",
                exception=GreenKuboError,
            )

        if self.method not in self._supported_methods:
            self.logger.error(
                (
                    f"Unknown method '{self.method}'. Possible methods "
                    f"are: {', '.join(self._supported_methods)}."
                ),
                exception=GreenKuboError,
            )

        if self.coefficient not in self._supported_coefficients:
            self.logger.error(
                (
                    f"Unknown transport coefficient '{self.coefficient}'. "
                    "Currently supported coefficients are: "
                    f"{', '.join(self._supported_coefficients)}."
                ),
                exception=GreenKuboError,
            )

        if not 0.0 <= self.fit_start < self.fit_stop <= 1.0:
            self.logger.error(
                (
                    f"The fit window fractions fit_start={self.fit_start} "
                    f"and fit_stop={self.fit_stop} must satisfy "
                    "0 <= fit_start < fit_stop <= 1."
                ),
                exception=GreenKuboError,
            )

        if self.n_blocks < 2:
            self.logger.error(
                (
                    f"The number of blocks n_blocks={self.n_blocks} must "
                    "be at least 2 so that a block-averaged standard "
                    "error of the diffusion coefficient can be computed."
                ),
                exception=GreenKuboError,
            )

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms of the Green-Kubo analysis."""
        return self.topology.n_atoms

    @property
    def fit_slice(self) -> slice:
        """slice: The index slice of the plateau fit window."""
        start = int(round(self.fit_start * self.window_size))
        stop = int(round(self.fit_stop * self.window_size))

        # guarantee at least one point in the fit window
        stop = max(stop, start + 1)

        return slice(start, stop + 1)

    @property
    def diffusion_coefficient_cm2_per_s(self) -> float:
        """float: The plateau diffusion coefficient in cm^2 / s."""
        return self.diffusion_coefficient * M2_PER_S_TO_CM2_PER_S

    @property
    def diffusion_coefficient_stderr_cm2_per_s(self) -> float:
        """
        float: The block-averaged standard error of the diffusion
        coefficient in cm^2 / s.
        """
        return self.diffusion_coefficient_stderr * M2_PER_S_TO_CM2_PER_S

    @timeit_in_class
    def run(
        self
    ) -> Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]:
        """
        Runs the Green-Kubo analysis.

        This method reads the velocities, calculates the un-normalized
        velocity auto-correlation function, integrates it to the
        running diffusion coefficient and extracts the plateau estimate
        of the diffusion coefficient. Its block-averaged standard error
        is obtained by running the same pipeline independently on
        ``n_blocks`` contiguous blocks of the trajectory.

        This method will display a progress bar by default. This can be
        disabled by setting with_progress_bar to False.

        Returns
        -------
        lag_times : Np1DNumberArray
            The lag times ``lag * time_step`` in ps for the lags
            ``0..window_size``.
        cvv : Np1DNumberArray
            The un-normalized velocity auto-correlation function in
            (Angstrom / s)^2.
        d_running : Np1DNumberArray
            The running diffusion coefficient in m^2 / s.
        """
        velocities = self._load_velocities()

        if self.method == "fft":
            self.cvv = velocity_acf_fft(velocities, self.window_size)
            self.n_origins = self.n_frames
        else:
            self.cvv = velocity_acf_direct(
                velocities,
                self.window_size,
                self.gap,
            )
            self.n_origins = len(range(0, self.n_frames, self.gap))

        self.lag_times = np.arange(self.window_size + 1) * self.time_step

        integral = cumulative_trapezoid(self.cvv, self.time_step)
        prefactor = ANGSTROM2_PER_S2_PS_TO_M2_PER_S / 3.0
        self.d_running = prefactor * integral

        plateau = self.d_running[self.fit_slice]
        self.diffusion_coefficient = float(np.mean(plateau))
        # the spread of the (strongly auto-correlated) running integral
        # over the fit window only measures the flatness of the plateau,
        # NOT the run-to-run statistical error of D
        self.diffusion_coefficient_plateau_spread = float(np.std(plateau))

        self._block_average(velocities)

        return self.lag_times, self.cvv, self.d_running

    def _block_average(self, velocities: NpnDNumberArray) -> None:
        """
        Estimates the statistical uncertainty of the diffusion
        coefficient by block averaging.

        The trajectory is split into ``n_blocks`` (clamped to fit the
        correlation window) contiguous, roughly-equal, non-overlapping
        blocks, the same velocity-ACF -> running-integral -> plateau
        pipeline is run independently on each block and the reported
        :py:attr:`diffusion_coefficient_stderr` is the standard error of
        the mean over the per-block diffusion coefficients
        ``std(D_b, ddof=1) / sqrt(n_blocks)``.

        Parameters
        ----------
        velocities : NpnDNumberArray
            The selected velocities of shape ``(n_frames, n_sel, 3)``.
        """
        n_blocks = self._resolve_n_blocks()

        blocks = np.array_split(velocities, n_blocks, axis=0)
        self.block_diffusion_coefficients = np.array(
            [self._plateau_diffusion(block) for block in blocks],
            dtype=np.float64,
        )

        self.diffusion_coefficient_stderr = float(
            np.std(self.block_diffusion_coefficients, ddof=1)
            / np.sqrt(n_blocks)
        )

    def _resolve_n_blocks(self) -> int:
        """
        Returns the number of blocks actually used for the block
        average.

        Each block must be longer than the correlation window so that
        every lag has at least one time origin. The requested
        ``n_blocks`` is clamped down (with a warning) to the largest
        number of blocks that still satisfy ``block_length >
        window_size``; if not even two such blocks fit a GreenKuboError
        is raised.

        Returns
        -------
        int
            The number of blocks used for the block average.

        Raises
        ------
        GreenKuboError
            If the trajectory is too short to form at least two blocks
            that each cover the correlation window.
        """
        max_blocks = self.n_frames // (self.window_size + 1)

        if self.n_blocks <= max_blocks:
            return self.n_blocks

        if max_blocks < 2:
            self.logger.error(
                (
                    f"The trajectory ({self.n_frames} frames, correlation "
                    f"window {self.window_size}) is too short to be split "
                    f"into {self.n_blocks} blocks that each cover the "
                    "correlation window. At least 2 blocks longer than "
                    f"window_size = {self.window_size} frames are needed, "
                    f"but only {max_blocks} fit(s). Reduce n_blocks or the "
                    "window size, or provide a longer trajectory."
                ),
                exception=GreenKuboError,
            )

        self.logger.warning(
            f"The requested number of blocks n_blocks={self.n_blocks} is "
            f"too large for a trajectory of {self.n_frames} frames with "
            f"correlation window {self.window_size}; each block must be "
            f"longer than the window. Clamping to {max_blocks} blocks."
        )

        return max_blocks

    def _plateau_diffusion(self, velocities: NpnDNumberArray) -> float:
        """
        Calculates the plateau diffusion coefficient of a single block
        of velocities with the same pipeline as the full trajectory.

        Parameters
        ----------
        velocities : NpnDNumberArray
            The selected velocities of one block with shape
            ``(block_length, n_sel, 3)``.

        Returns
        -------
        float
            The plateau diffusion coefficient of the block in m^2 / s.
        """
        if self.method == "fft":
            cvv = velocity_acf_fft(velocities, self.window_size)
        else:
            cvv = velocity_acf_direct(velocities, self.window_size, self.gap)

        prefactor = ANGSTROM2_PER_S2_PS_TO_M2_PER_S / 3.0
        d_running = prefactor * cumulative_trapezoid(cvv, self.time_step)

        return float(np.mean(d_running[self.fit_slice]))

    def _load_velocities(self) -> NpnDNumberArray:
        """
        Reads all selected velocities into a dense float64 array.

        Returns
        -------
        NpnDNumberArray
            The selected velocities of shape
            ``(n_frames, n_sel, 3)`` in Angstrom / s.
        """
        return np.stack(list(self._velocities()))

    def _velocities(self) -> Generator[Np2DNumberArray, None, None]:
        """
        Yields the selected velocities of all frames.

        Dispatches to the raw fast-path stream if the analysis was
        constructed from a velocity TrajectoryReader and to the
        AtomicSystem based stream otherwise. Both streams yield
        bit-identical float64 arrays.

        Returns
        -------
        Generator[Np2DNumberArray, None, None]
            The selected velocities of all frames.
        """
        if self._raw_reader is not None:
            return self._raw_velocities()

        return self._frame_velocities()

    def _raw_velocities(self) -> Generator[Np2DNumberArray, None, None]:
        """
        Yields the selected velocities of all frames from the raw
        fast-path reader.

        Yields
        ------
        Np2DNumberArray
            The selected velocities of one frame with shape
            ``(n_sel, 3)``.

        Raises
        ------
        GreenKuboError
            If a frame does not provide velocities for all atoms.
        """
        n_atoms = self.n_atoms
        indices = self._target_indices_intp

        for values, _cell in tqdm(
            self._raw_reader.raw_frame_generator(),
            total=self.n_frames,
            disable=not config.with_progress_bar):

            if values.shape[0] != n_atoms:
                self.logger.error(
                    (
                        "A frame of the velocity trajectory does not "
                        f"provide velocities for all {n_atoms} atoms. "
                        "Please provide a velocity trajectory (e.g. "
                        ".vel files)."
                    ),
                    exception=GreenKuboError,
                )

            yield np.asarray(values, dtype=np.float64)[indices]

    def _frame_velocities(self) -> Generator[Np2DNumberArray, None, None]:
        """
        Yields the selected velocities of all frames from the
        AtomicSystem based stream.

        Yields
        ------
        Np2DNumberArray
            The selected velocities of one frame with shape
            ``(n_sel, 3)``.

        Raises
        ------
        GreenKuboError
            If a frame does not provide velocities for all atoms.
        """
        frames = itertools.chain([self._first_frame], self._frame_generator)

        for frame in tqdm(
            frames,
            total=self.n_frames,
            disable=not config.with_progress_bar):
            vel = np.asarray(frame.vel, dtype=np.float64)

            if vel.ndim != 2 or vel.shape[0] != self.n_atoms:
                self.logger.error(
                    (
                        "A frame of the velocity trajectory does not "
                        f"provide velocities for all {self.n_atoms} "
                        "atoms. Please provide a velocity trajectory "
                        "(e.g. .vel files)."
                    ),
                    exception=GreenKuboError,
                )

            yield vel[self.target_indices]
