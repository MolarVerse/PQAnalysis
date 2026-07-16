"""
A module containing the VACF class. The VACF class is used to
calculate the normalized velocity auto-correlation function
of a target selection of a given velocity trajectory. In the
charge-flux mode the atomic velocities are weighted with
(possibly time-dependent) atomic partial charges, which turns
the auto-correlation function into a charge-flux (current)
auto-correlation function.

The implementation is a port of the legacy ``FreqCalc`` and
``Fluxfreqcalc`` tools of the ``thh_tools`` collection. Both
legacy tools share the identical sliding-time-origin estimator -
they only differ in the charge weighting of the velocities at
read time.
"""

import itertools
import logging

# 3rd party imports
import numpy as np
from beartype.typing import Generator, Tuple

# local absolute imports
from PQAnalysis.types import (
    Np1DNumberArray,
    Np2DNumberArray,
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
from .exceptions import VACFError
from ._raw_charge_reader import RawChargeTrajectoryReader

try:
    from ._vacf_kernel import (  # pylint: disable=import-error
        accumulate_frame,
        weight_frame,
    )
except ModuleNotFoundError:
    from ._vacf_kernel_py import accumulate_frame, weight_frame



class VACF:

    """
    A class for calculating the normalized velocity (or charge-flux)
    auto-correlation function of a target selection of a velocity
    trajectory.

    The auto-correlation function is estimated with sliding time
    origins. A new origin is spawned every ``gap`` frames as long as a
    full correlation window of ``window_size`` frames can still be
    accommodated by the trajectory, i.e. for all origin frames
    ``t0 <= stop_frame`` with
    ``stop_frame = floor((n_frames - window_size) / gap) * gap``
    (frames are counted starting at 1). Like in the legacy tools, a
    ``stop_frame`` of exactly zero (i.e. ``n_frames == window_size``)
    is reset to one, so that for ``gap == 1`` a single origin is
    spawned at the first frame and only the final lag bin of the
    correlation function stays zero. Each origin is normalized by
    its own aggregate squared velocity norm, so that the correlation
    function starts at exactly one:

    ``C(lag) = (1 / n_origins) * sum_i [sum_j v_j(t0_i + lag) . v_j(t0_i)]
    / [sum_j |v_j(t0_i)|^2]``

    where the sum over ``j`` runs over all selected atoms and the sum
    over ``i`` over all time origins. This reproduces the legacy
    ``FreqCalc``/``Fluxfreqcalc`` estimator exactly, including the
    requirement that ``window_size`` is an integer multiple of ``gap``.

    In the charge-flux mode (legacy ``Fluxfreqcalc``) every velocity is
    replaced by ``q_j(t) * v_j(t)`` at read time, either with static
    charges or with a charge trajectory read in lockstep with the
    velocity trajectory.

    Alternatively, a fast estimator based on the Wiener-Khinchin
    theorem can be selected with ``method='fft'``. Note that this is a
    slightly different (denser-origin) estimator: every frame serves as
    a time origin (the ``gap`` parameter is ignored), each lag is
    averaged over its actual number of origins and the normalization is
    performed with the aggregate mean squared velocity instead of
    per-origin norms. The default ``method='direct'`` is legacy-exact.
    """

    _gap_default = 1
    _method_default = "direct"

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        traj: Trajectory | TrajectoryReader,
        window_size: PositiveInt,
        time_step: PositiveReal,
        target_species: SelectionCompatible = None,
        gap: PositiveInt | None = None,
        charges: Np1DNumberArray | None = None,
        charge_traj: Trajectory | TrajectoryReader | None = None,
        method: str | None = None,
        use_full_atom_info: bool | None = False,
    ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The velocity trajectory to analyze. If a TrajectoryReader
            is provided, the trajectory is read lazily frame by frame.
            All frames must provide velocities.
        window_size : PositiveInt
            The correlation window length in frames. The correlation
            function is calculated for the lags ``0..window_size``.
        time_step : PositiveReal
            The time step between two frames in ps. It is only used to
            build the time axis of the results.
        target_species : SelectionCompatible, optional
            The target species of the VACF analysis, by default None
            (all atoms).
        gap : PositiveInt | None, optional
            The spacing between two time origins in frames, by
            default None (1). ``window_size`` must be an integer
            multiple of ``gap``.
        charges : Np1DNumberArray | None, optional
            Static atomic partial charges for the charge-flux mode, one
            charge per atom of the full system, by default None.
        charge_traj : Trajectory | TrajectoryReader | None, optional
            A charge trajectory for the charge-flux mode, read in
            lockstep with the velocity trajectory, by default None.
        method : str | None, optional
            The estimator to use, either ``direct`` (legacy-exact
            sliding origins) or ``fft`` (denser-origin Wiener-Khinchin
            estimator, ignores ``gap``), by default None (``direct``).
        use_full_atom_info : bool | None, optional
            Whether to use the full atom information for the target
            selection, by default None (False).

        Raises
        ------
        VACFError
            If the time step is not positive.
        VACFError
            If the method is unknown.
        VACFError
            If both static charges and a charge trajectory are given.
        VACFError
            If the window size is not an integer multiple of the gap.
        VACFError
            If the trajectory is too short to place a single time
            origin.
        VACFError
            If the target selection does not select any atoms.
        VACFError
            If the number of static charges does not match the number
            of atoms.
        """

        self.window_size = window_size
        self.time_step = time_step
        self.gap = gap if gap is not None else self._gap_default
        self.method = (
            method if method is not None else self._method_default
        ).lower()

        if self.time_step <= 0.0:
            self.logger.error(
                "The time_step must be a positive real number.",
                exception=VACFError,
            )

        if self.method not in ("direct", "fft"):
            self.logger.error(
                (
                    f"Unknown method '{self.method}'. "
                    "Possible methods are: direct, fft."
                ),
                exception=VACFError,
            )

        if charges is not None and charge_traj is not None:
            self.logger.error(
                (
                    "Only one charge source can be used for the "
                    "charge-flux mode: either static charges or a "
                    "charge trajectory."
                ),
                exception=VACFError,
            )

        if self.window_size % self.gap != 0:
            self.logger.error(
                (
                    f"The window_size {self.window_size} must be an "
                    f"integer multiple of the gap {self.gap} for the "
                    "sliding-origin machinery."
                ),
                exception=VACFError,
            )

        if use_full_atom_info is None:
            use_full_atom_info = False
        self.use_full_atom_info = use_full_atom_info

        ############################################
        # Initialize trajectory iterator/generator #
        ############################################

        self._raw_reader = None
        self._raw_charge_reader = None
        self._charge_frame_generator = None
        self._charge_value_stream = None

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
                self._frame_generator = None
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
                exception=VACFError,
            )

        if self.method == "fft":
            # every lag 0..window_size needs at least one origin
            trajectory_too_short = (
                self.n_frames < self.window_size + self.gap
            )
        else:
            # legacy FreqCalc semantics: stop_frame == 0 is reset to 1
            # (see the stop_frame property), so a trajectory of exactly
            # window_size frames still spawns one origin for gap == 1
            trajectory_too_short = self.stop_frame < self.gap

        if trajectory_too_short:
            self.logger.error(
                (
                    f"The trajectory contains only {self.n_frames} "
                    "frame(s), but at least window_size + gap = "
                    f"{self.window_size + self.gap} frames are needed "
                    "to place a single time origin (or exactly "
                    "window_size frames for the direct method with "
                    "gap == 1)."
                ),
                exception=VACFError,
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
                exception=VACFError,
            )

        self._target_indices_intp = np.ascontiguousarray(
            self.target_indices,
            dtype=np.intp,
        )

        ##########################
        # Initialize charge mode #
        ##########################

        self._static_charges = None

        if charges is not None:
            if len(charges) != self.n_atoms:
                self.logger.error(
                    (
                        f"The number of static charges {len(charges)} "
                        "does not match the number of atoms "
                        f"{self.n_atoms} of the system."
                    ),
                    exception=VACFError,
                )

            self._static_charges = np.asarray(
                charges,
                dtype=np.float64,
            )[self.target_indices]

        if charge_traj is not None:
            self._init_charge_traj(charge_traj)

        self.flux = charges is not None or charge_traj is not None

        # result dummy init
        self.vacf = np.zeros(self.window_size + 1)
        self.n_origins = 0

    def _init_charge_traj(
        self,
        charge_traj: Trajectory | TrajectoryReader,
    ) -> None:
        """
        Sets up the lockstep charge trajectory iterator.

        Parameters
        ----------
        charge_traj : Trajectory | TrajectoryReader
            The charge trajectory to read in lockstep with the
            velocity trajectory.

        Raises
        ------
        VACFError
            If the number of charge frames does not match the number
            of velocity frames.
        """
        if (
            isinstance(charge_traj, TrajectoryReader) and
            charge_traj.traj_format == TrajectoryFormat.CHARGE
        ):
            # additive fast path: stream the raw float64 charge values
            # without building an AtomicSystem per frame
            # (bit-identical values)
            self._raw_charge_reader = RawChargeTrajectoryReader(
                charge_traj.filenames,
                md_format=charge_traj.md_format,
            )
            n_charge_frames = self._raw_charge_reader.count_frames()
        elif isinstance(charge_traj, TrajectoryReader):
            n_charge_frames = sum(
                charge_traj.calculate_number_of_frames_per_file()
            )
            self._charge_frame_generator = charge_traj.frame_generator()
        else:
            n_charge_frames = len(charge_traj)
            self._charge_frame_generator = iter(charge_traj)

        if n_charge_frames != self.n_frames:
            self.logger.error(
                (
                    f"The charge trajectory contains {n_charge_frames} "
                    "frame(s), but the velocity trajectory contains "
                    f"{self.n_frames} frame(s). Both trajectories have "
                    "to be in lockstep."
                ),
                exception=VACFError,
            )

        self._charge_value_stream = self._charge_values()

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms of the VACF analysis."""
        return self.topology.n_atoms

    @property
    def stop_frame(self) -> int:
        """int: The last frame (1-based) at which a time origin is spawned."""
        stop_frame = (
            (self.n_frames - self.window_size) // self.gap
        ) * self.gap

        # legacy FreqCalc reset (process.c): a trajectory of exactly
        # window_size frames still spawns a single origin at frame 1
        # for gap == 1; the final lag bin then stays zero
        if stop_frame == 0:
            stop_frame = 1

        return stop_frame

    @timeit_in_class
    def run(self) -> Tuple[Np1DNumberArray, Np1DNumberArray]:
        """
        Runs the VACF analysis.

        Returns
        -------
        time : Np1DNumberArray
            The lag times ``lag * time_step`` for the lags
            ``0..window_size``.
        vacf : Np1DNumberArray
            The normalized (charge-flux weighted) velocity
            auto-correlation function; ``vacf[0]`` is exactly one for
            the direct method.
        """
        if self.method == "fft":
            self.vacf = self._run_fft()
        else:
            self.vacf = self._run_direct()

        time = np.arange(self.window_size + 1) * self.time_step

        return time, self.vacf

    def _velocities(self) -> Generator[Np2DNumberArray, None, None]:
        """
        Returns the (charge weighted) selected velocity stream.

        Dispatches to the raw fast-path stream if the analysis was
        constructed from a velocity TrajectoryReader and to the
        AtomicSystem based stream otherwise. Both streams yield
        bit-identical float64 arrays.

        Returns
        -------
        Generator[Np2DNumberArray, None, None]
            The (charge weighted) selected velocities of all frames.
        """
        if self._raw_reader is not None:
            return self._raw_weighted_velocities()

        return self._weighted_velocities()

    def _raw_weighted_velocities(
        self
    ) -> Generator[Np2DNumberArray, None, None]:
        """
        Yields the (charge weighted) selected velocities of all frames
        from the raw fast-path reader.

        The raw float32 values of every frame are cast to float64,
        reduced to the target selection and weighted with the frame
        charges - bit-identical to the AtomicSystem based stream of
        :py:meth:`_weighted_velocities`.

        Yields
        ------
        Np2DNumberArray
            The selected velocities of one frame with shape
            ``(n_target_atoms, 3)``.

        Raises
        ------
        VACFError
            If a frame does not provide velocities for all atoms.
        VACFError
            If a charge frame does not provide charges for all atoms.
        """
        n_atoms = self.n_atoms
        indices = self._target_indices_intp

        for values, _cell in self._raw_reader.raw_frame_generator():

            if values.shape[0] != n_atoms:
                self.logger.error(
                    (
                        "A frame of the velocity trajectory does not "
                        f"provide velocities for all {n_atoms} "
                        "atoms. Please provide a velocity trajectory "
                        "(e.g. .vel files)."
                    ),
                    exception=VACFError,
                )

            yield weight_frame(values, indices, self._frame_charges())

    def _weighted_velocities(
        self
    ) -> Generator[Np2DNumberArray, None, None]:
        """
        Yields the (charge weighted) selected velocities of all frames.

        The velocities are accumulated in float64 even though the
        underlying trajectory reader parses them as float32.

        Yields
        ------
        Np2DNumberArray
            The selected velocities of one frame with shape
            ``(n_target_atoms, 3)``.

        Raises
        ------
        VACFError
            If a frame does not provide velocities for all atoms.
        VACFError
            If a charge frame does not provide charges for all atoms.
        """
        frames = itertools.chain([self._first_frame], self._frame_generator)

        for frame in frames:
            vel = np.asarray(frame.vel, dtype=np.float64)

            if vel.ndim != 2 or vel.shape[0] != self.n_atoms:
                self.logger.error(
                    (
                        "A frame of the velocity trajectory does not "
                        f"provide velocities for all {self.n_atoms} "
                        "atoms. Please provide a velocity trajectory "
                        "(e.g. .vel files)."
                    ),
                    exception=VACFError,
                )

            vel = vel[self.target_indices]

            frame_charges = self._frame_charges()

            if frame_charges is not None:
                vel = vel * frame_charges[:, None]

            yield vel

    def _frame_charges(self) -> Np1DNumberArray | None:
        """
        Returns the selected charges of the current frame.

        Returns
        -------
        Np1DNumberArray | None
            The static charges, the charges of the next frame of the
            lockstep charge trajectory or None outside the charge-flux
            mode.
        """
        if self._static_charges is not None:
            return self._static_charges

        if self._charge_value_stream is not None:
            return self._next_charges()

        return None

    def _charge_values(self) -> Generator[Np1DNumberArray, None, None]:
        """
        Yields the full-system float64 charges of all charge frames.

        Yields
        ------
        Np1DNumberArray
            The charges of one frame of the charge trajectory.
        """
        if self._raw_charge_reader is not None:
            for values, _cell in self._raw_charge_reader.raw_frame_generator():
                yield values
        else:
            for charge_frame in self._charge_frame_generator:
                yield np.asarray(charge_frame.charges, dtype=np.float64)

    def _next_charges(self) -> Np1DNumberArray:
        """
        Reads the next charge frame of the lockstep charge trajectory.

        Returns
        -------
        Np1DNumberArray
            The selected charges of the next charge frame.

        Raises
        ------
        VACFError
            If the charge trajectory is exhausted or a charge frame
            does not provide charges for all atoms.
        """
        charge = next(self._charge_value_stream, None)

        if charge is None:
            self.logger.error(
                (
                    "The charge trajectory provides fewer frames than "
                    "the velocity trajectory."
                ),
                exception=VACFError,
            )

        if charge.ndim != 1 or charge.shape[0] != self.n_atoms:
            self.logger.error(
                (
                    "A frame of the charge trajectory does not provide "
                    f"charges for all {self.n_atoms} atoms. Please "
                    "provide a charge trajectory (e.g. .chrg files)."
                ),
                exception=VACFError,
            )

        return charge[self.target_indices]

    def _run_direct(self) -> Np1DNumberArray:
        """
        Runs the legacy-exact sliding-origin estimator.

        A new time origin is spawned every ``gap`` frames up to
        ``stop_frame``. Every frame each active origin ``i``
        contributes ``sum_j v_j(t) . v_j(t0_i) / sum_j |v_j(t0_i)|^2``
        to the lag ``t - t0_i``. An origin is retired after it has
        contributed to the lag ``window_size``. Finally the correlation
        function is divided by the total number of spawned origins.

        Returns
        -------
        Np1DNumberArray
            The normalized auto-correlation function for the lags
            ``0..window_size``.

        Raises
        ------
        VACFError
            If the aggregate squared velocity norm of an origin is zero.
        """
        window_size = self.window_size
        gap = self.gap
        stop_frame = self.stop_frame

        n_slots = window_size // gap + 1
        n_target = len(self.target_indices)

        corr = np.zeros(window_size + 1, dtype=np.float64)
        origin_vel = np.zeros((n_slots, n_target, 3), dtype=np.float64)
        origin_norm = np.zeros(n_slots, dtype=np.float64)
        origin_frame = np.zeros(n_slots, dtype=np.longlong)
        n_active = 0
        n_origins = 0

        for frame_number, vel in enumerate(self._velocities(), 1):

            spawn = frame_number % gap == 0 and frame_number <= stop_frame

            n_active = accumulate_frame(
                corr,
                origin_vel,
                origin_norm,
                origin_frame,
                n_active,
                vel,
                frame_number,
                spawn,
                window_size,
            )

            if n_active < 0:
                self.logger.error(
                    (
                        "The aggregate squared velocity norm of the "
                        f"time origin at frame {frame_number} is "
                        "zero. The normalized VACF is not defined."
                    ),
                    exception=VACFError,
                )

            if spawn:
                n_origins += 1

        self.n_origins = n_origins

        return corr / n_origins

    def _run_fft(self) -> Np1DNumberArray:
        """
        Runs the Wiener-Khinchin (FFT) estimator.

        This is a denser-origin estimator: every frame serves as a time
        origin, the raw aggregate auto-correlation
        ``S(lag) = sum_t sum_j v_j(t) . v_j(t + lag)`` is calculated
        per atom and component via FFT, each lag is divided by its
        number of origins ``n_frames - lag`` and the result is
        normalized with its lag-zero value, so that ``C(0) = 1``.

        In contrast to the direct method, the normalization uses the
        aggregate mean squared velocity instead of per-origin norms and
        the ``gap`` parameter is ignored. All velocities are kept in
        memory.

        Returns
        -------
        Np1DNumberArray
            The normalized auto-correlation function for the lags
            ``0..window_size``.

        Raises
        ------
        VACFError
            If the aggregate squared velocity norm of the trajectory
            is zero.
        """
        vel = np.stack(list(self._velocities()))
        n_frames = vel.shape[0]

        n_fft = 2 * n_frames
        spectrum = np.fft.rfft(vel, n=n_fft, axis=0)
        autocorr = np.fft.irfft(
            spectrum * np.conj(spectrum),
            n=n_fft,
            axis=0,
        )[:self.window_size + 1]

        raw = np.sum(autocorr.real, axis=(1, 2))
        counts = n_frames - np.arange(self.window_size + 1)
        raw = raw / counts

        if raw[0] == 0.0:
            self.logger.error(
                (
                    "The aggregate squared velocity norm of the "
                    "trajectory is zero. The normalized VACF is "
                    "not defined."
                ),
                exception=VACFError,
            )

        self.n_origins = n_frames

        return raw / raw[0]
