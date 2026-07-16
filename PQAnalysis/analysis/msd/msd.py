"""
A module containing the MSD class. The MSD class is used to
calculate the mean square displacement (MSD) of a selection
of atoms using multiple time origins on a sliding window.
The MSD is a measure of the average squared distance
particles travel within a given correlation time and gives
access to the self-diffusion coefficient via the Einstein
relation.
"""

import dataclasses
import itertools
import logging

# 3rd party imports
import numpy as np

from beartype.typing import Dict, Tuple
from scipy.stats import linregress
from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis import config
from PQAnalysis.core import Cell
from PQAnalysis.types import (
    Np1DNumberArray,
    Np2DNumberArray,
    PositiveInt,
    PositiveReal,
)
from PQAnalysis.traj import Trajectory
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import TrajectoryReader
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import MSDError

#: float: Conversion factor from Angstrom^2/ps to m^2/s.
ANGSTROM2_PER_PS_TO_M2_PER_S = 1.0e-8



@dataclasses.dataclass(frozen=True)
class MSDDiffusionFit:

    """
    A container for the result of a linear diffusion fit of an
    MSD component.

    The slope related quantities are given in Angstrom^2/ps,
    the diffusion coefficient and its standard error are given
    in m^2/s.
    """

    label: str
    slope: float
    slope_stderr: float
    intercept: float
    r_squared: float
    diffusion_coefficient: float
    diffusion_coefficient_stderr: float



class MSD:

    """
    A class for calculating the mean square displacement (MSD)
    of a target selection using multiple simultaneously active
    time origins.

    A new time origin is spawned every ``gap`` frames until
    ``window // gap`` origins are active. Every frame, each
    active origin accumulates the squared displacement of every
    selected atom relative to the origin position into the lag
    bin given by the frame distance to the origin. Displacements
    are unwrapped with a running minimum image convention applied
    to the per-frame displacement vectors, so that trajectories
    wrapped into the simulation box are handled correctly. When
    the oldest origin has covered the full window, it accumulates
    its final lag term and is replaced by a new origin (or
    dropped without replacement once no new origins may spawn
    anymore).

    This is a port of the ``Diffcalc`` tool of thh_tools and
    reproduces its results exactly, including its normalization
    convention: every lag bin is divided by the number of
    selected atoms times ``total_origins``, where
    ``total_origins = stop_frame // gap`` and
    ``stop_frame = (n_frames - window) // gap * gap``. Note that
    this normalization intentionally follows the legacy code
    also for ``n_start > 0``. The legacy boundary case of a
    trajectory of exactly ``window`` frames with ``gap == 1`` is
    also kept: a single time origin spawns at the first frame
    and the final lag bin (``lag == window``), which can never
    be sampled, is written as exactly 0.0 (a warning is
    emitted). Any shorter trajectory raises an
    :py:class:`~PQAnalysis.analysis.msd.exceptions.MSDError`.

    If a time step is given, the trailing ``fit_window`` points
    of the MSD components are fitted linearly to extract
    self-diffusion coefficients via the Einstein relation
    (slope / (2 * dimensionality), converted to m^2/s).

    The MSD class can be initialized with either a trajectory
    object or via a TrajectoryReader object. The
    TrajectoryReader allows for lazy loading of the trajectory,
    which is useful for large trajectories that do not fit into
    memory.
    """

    _use_full_atom_default = False
    _window_default = 1000
    _gap_default = 10
    _n_start_default = 0

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        traj: Trajectory | TrajectoryReader,
        target_species: SelectionCompatible,
        use_full_atom_info: bool | None = False,
        window: PositiveInt | None = None,
        gap: PositiveInt | None = None,
        n_start: int | None = None,
        time_step: PositiveReal | None = None,
        fit_window: PositiveInt | None = None,
    ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory to analyze. If a TrajectoryReader is
            provided, the trajectory is read frame by frame via
            a frame_generator.
        target_species : SelectionCompatible
            The target species of the MSD analysis.
        use_full_atom_info : bool | None, optional
            Whether to use the full atom information of the
            trajectory or not, by default None (False).
        window : PositiveInt | None, optional
            The correlation window size in frames,
            by default None (1000).
        gap : PositiveInt | None, optional
            The gap between two time origins in frames,
            by default None (10).
        n_start : int | None, optional
            The first frame (1-based frame counter) at which
            processing starts, by default None (0). Frames
            before n_start are read (to keep the unwrapping
            continuous) but do not contribute to the MSD.
        time_step : PositiveReal | None, optional
            The time step between two frames in ps. If given,
            diffusion coefficients are calculated from a linear
            fit of the MSD tail, by default None.
        fit_window : PositiveInt | None, optional
            The number of trailing MSD points used for the
            diffusion fit, by default None (last 20% of the
            window).

        Raises
        ------
        MSDError
            If n_start is negative.
        MSDError
            If the window is not a multiple of the gap.
        MSDError
            If time_step is not positive.
        MSDError
            If fit_window is smaller than 2.
        MSDError
            If fit_window is larger than window + 1.
        MSDError
            If the trajectory is empty.
        MSDError
            If the target selection is empty.
        MSDError
            If the trajectory is too short to establish at
            least one full window (n_frames < window + gap and
            not the legacy single-origin case of exactly
            window frames with gap == 1).
        MSDError
            If n_start is larger than stop_frame, so that no
            time origin could spawn.

        See Also
        --------
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        :py:class:`~PQAnalysis.topology.selection.Selection`
        :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`
        """

        ##############
        # dummy init #
        ##############

        self.lags = np.array([])
        self.msd_x = np.array([])
        self.msd_y = np.array([])
        self.msd_z = np.array([])
        self.msd_tot = np.array([])
        self.fit_results = None
        self._msd_accumulator = None

        #####################################################
        # Initialize parameters with default values if None #
        #####################################################

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info

        self.window = window if window is not None else self._window_default
        self.gap = gap if gap is not None else self._gap_default
        self.n_start = n_start if n_start is not None else self._n_start_default
        self.time_step = time_step

        if fit_window is None:
            self.fit_window = max(2, self.window // 5)
        else:
            self.fit_window = fit_window

        self._check_parameters()

        ############################################
        # Initialize Trajectory iterator/generator #
        ############################################

        if isinstance(traj, TrajectoryReader):
            # lazy loading of trajectory from file(s)
            self.n_frames = sum(traj.calculate_number_of_frames_per_file())
            self.frame_generator = traj.frame_generator()
        elif len(traj) > 0:
            # use trajectory object as iterator
            self.n_frames = len(traj)
            self.frame_generator = iter(traj)
        else:
            self.logger.error(
                "Trajectory cannot be of length 0.",
                exception=MSDError
            )

        self.first_frame = next(self.frame_generator)

        if traj.topology is not None:
            self.topology = traj.topology
        else:
            self.topology = self.first_frame.topology

        ################################
        # Initialize Selection objects #
        ################################

        self.target_species = target_species
        self.target_selection = Selection(target_species)

        self.target_indices = self.target_selection.select(
            self.topology,
            self.use_full_atom_info
        )

        if len(self.target_indices) == 0:
            self.logger.error(
                "The target selection does not select any atoms.",
                exception=MSDError
            )

        ##########################################
        # Setup time origin bookkeeping (legacy) #
        ##########################################

        self._setup_origin_bookkeeping()

    def _check_parameters(self):
        """
        Checks the consistency of the setup parameters.

        This method is called by the __init__ method of the MSD
        class after all parameters have been initialized.

        Raises
        ------
        MSDError
            If n_start is negative.
        MSDError
            If the window is not a multiple of the gap.
        MSDError
            If time_step is not positive.
        MSDError
            If fit_window is smaller than 2.
        MSDError
            If fit_window is larger than window + 1.
        """

        if self.n_start < 0:
            self.logger.error(
                "n_start must be a non-negative integer.",
                exception=MSDError
            )

        if self.window % self.gap != 0:
            self.logger.error(
                (
                    f"The window size {self.window} has to be an "
                    f"integer multiple of the gap {self.gap}."
                ),
                exception=MSDError
            )

        if self.time_step is not None and self.time_step <= 0.0:
            self.logger.error(
                "The time_step must be a positive real number.",
                exception=MSDError
            )

        if self.fit_window < 2:
            self.logger.error(
                (
                    f"The fit_window {self.fit_window} must be at "
                    "least 2 to perform a linear diffusion fit."
                ),
                exception=MSDError
            )

        if self.fit_window > self.window + 1:
            self.logger.error(
                (
                    f"The fit_window {self.fit_window} cannot be larger "
                    f"than window + 1 = {self.window + 1}."
                ),
                exception=MSDError
            )

    def _setup_origin_bookkeeping(self):
        """
        Sets up the time origin bookkeeping of the MSD analysis.

        This method is called by the __init__ method of the MSD
        class. It follows the legacy Diffcalc conventions: origins
        may only spawn until stop_frame, so that every origin can
        cover the full window, and the total number of origins used
        for the normalization is stop_frame // gap.

        One boundary case is kept from the legacy Diffcalc code:
        for gap == 1 a trajectory of exactly window frames spawns
        a single time origin at the first frame (stop_frame is
        clamped from 0 to 1). In this case the final lag bin
        (lag == window) can never be sampled and is written as
        exactly 0.0, matching the legacy tool. A warning is
        emitted, as this zero bin would bias a diffusion fit that
        includes it.

        Raises
        ------
        MSDError
            If the trajectory is too short to establish at
            least one full window (n_frames < window + gap and
            not the legacy single-origin case of exactly
            window frames with gap == 1).
        MSDError
            If n_start is larger than stop_frame, so that no
            time origin could spawn.
        """

        self.n_origins_max = self.window // self.gap

        self.stop_frame = (
            (self.n_frames - self.window) // self.gap * self.gap
        )

        if self.stop_frame == 0 and self.gap == 1:
            # legacy Diffcalc branch: with gap == 1 a trajectory of
            # exactly window frames spawns a single time origin at
            # the first frame instead of raising the
            # too-short-trajectory error below
            self.stop_frame = 1

            self.logger.warning(
                (
                    f"The trajectory contains exactly window = "
                    f"{self.window} frames with a gap of 1. Following "
                    "the legacy Diffcalc convention a single time "
                    "origin is spawned at the first frame. The final "
                    f"lag bin (lag = {self.window}) can never be "
                    "sampled and is written as exactly 0.0; a "
                    "diffusion fit including this bin would be biased."
                )
            )

        self.total_origins = self.stop_frame // self.gap

        if self.total_origins < 1:
            self.logger.error(
                (
                    f"The trajectory with {self.n_frames} frames is too "
                    f"short to establish a window of {self.window} frames "
                    f"with a gap of {self.gap} frames. At least "
                    f"window + gap = {self.window + self.gap} frames "
                    "are required (or exactly window frames for "
                    "gap == 1)."
                ),
                exception=MSDError
            )

        if self.n_start > self.stop_frame:
            self.logger.error(
                (
                    f"The starting frame {self.n_start} is too large: "
                    "time origins only spawn at multiples of the gap "
                    f"{self.gap} up to stop_frame = (n_frames - window) "
                    f"// gap * gap = {self.stop_frame} (n_frames = "
                    f"{self.n_frames}, window = {self.window}), so no "
                    "time origin could spawn."
                ),
                exception=MSDError
            )

    @timeit_in_class
    def run(
        self
    ) -> Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]:
        """
        Runs the MSD analysis.

        This method runs the MSD analysis and returns the lag
        indices and the MSD components. If a time step was
        given, the diffusion coefficients are calculated from a
        linear fit of the MSD tail and stored in the
        fit_results attribute.

        This method will display a progress bar by default.
        This can be disabled by setting with_progress_bar to
        False.

        Returns
        -------
        lags : Np1DNumberArray
            The lag indices (in frames) of the MSD analysis,
            ranging from 0 to window.
        msd_x : Np1DNumberArray
            The x-component of the MSD in Angstrom^2.
        msd_y : Np1DNumberArray
            The y-component of the MSD in Angstrom^2.
        msd_z : Np1DNumberArray
            The z-component of the MSD in Angstrom^2.
        msd_tot : Np1DNumberArray
            The total MSD (x + y + z) in Angstrom^2.
        """

        self._calculate_msd()
        return self._finalize_run()

    def _calculate_msd(self):
        """
        Calculates the raw (unnormalized) MSD accumulators.

        This method is called by the run method of the MSD
        class. It streams over the frames of the trajectory and
        accumulates the per-axis squared displacements of all
        selected atoms for all active time origins into the
        corresponding lag bins. The per-frame displacements are
        unwrapped with a running minimum image convention so
        that box-wrapped trajectories are handled correctly.

        Raises
        ------
        MSDError
            If a frame does not provide positions for all atoms
            of the topology.
        """

        gap = self.gap

        msd = np.zeros((self.window + 1, 3))

        # unwrapped origin coordinates, oldest origin at index 0
        origins = np.zeros(
            (self.n_origins_max, len(self.target_indices), 3)
        )
        displacements = np.zeros_like(origins)
        n_active = 0
        last = 0

        # cumulative unwrapping shift of the selected atoms
        shift = np.zeros((len(self.target_indices), 3))
        prev_pos = None
        counter = 0

        for frame in tqdm(
            itertools.chain([self.first_frame], self.frame_generator),
            total=self.n_frames,
            disable=not config.with_progress_bar):

            counter += 1

            pos = np.asarray(frame.pos, dtype=np.float64)

            if pos.ndim != 2 or pos.shape[0] != self.n_atoms:
                self.logger.error(
                    (
                        f"Frame {counter} of the trajectory does not "
                        f"provide positions for all {self.n_atoms} "
                        "atoms of the topology. Please provide a "
                        "position trajectory (e.g. .xyz files) with "
                        "a consistent number of atoms."
                    ),
                    exception=MSDError
                )

            pos = pos[self.target_indices]

            if prev_pos is not None:
                shift += self._unwrap_shift(pos - prev_pos, frame.cell)

            prev_pos = pos

            unwrapped = pos + shift

            if counter < self.n_start:
                continue

            if counter % gap == 0:

                if (
                    n_active != self.n_origins_max
                    and counter <= self.stop_frame
                ):
                    # starting stage - add new origin
                    if n_active == 0:
                        last = counter

                    origins[n_active] = unwrapped
                    n_active += 1

                elif n_active > 0 and last + self.window == counter:
                    # oldest origin reached the full window:
                    # accumulate its final lag term
                    disp = unwrapped - origins[0]
                    msd[self.window] += np.einsum('ax,ax->x', disp, disp)

                    origins[:n_active - 1] = origins[1:n_active]

                    if counter > self.stop_frame:
                        # stopping stage - drop without replacement
                        n_active -= 1
                    else:
                        # running stage - replace by a new origin
                        origins[n_active - 1] = unwrapped

                    last += gap

            if n_active > 0:
                lags = counter - last - gap * np.arange(n_active)

                disp = np.subtract(
                    unwrapped,
                    origins[:n_active],
                    out=displacements[:n_active]
                )

                msd[lags] += np.einsum('oax,oax->ox', disp, disp)

        self._msd_accumulator = msd

    def _finalize_run(
        self
    ) -> Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]:
        """
        Finalizes the MSD analysis after running.

        This method is called by the run method of the MSD
        class. It normalizes the raw MSD accumulators by the
        number of selected atoms and the total number of time
        origins (legacy Diffcalc convention) and performs the
        diffusion fit if a time step was given.

        Returns
        -------
        lags : Np1DNumberArray
            The lag indices (in frames) of the MSD analysis.
        msd_x : Np1DNumberArray
            The x-component of the MSD in Angstrom^2.
        msd_y : Np1DNumberArray
            The y-component of the MSD in Angstrom^2.
        msd_z : Np1DNumberArray
            The z-component of the MSD in Angstrom^2.
        msd_tot : Np1DNumberArray
            The total MSD (x + y + z) in Angstrom^2.
        """

        norm = float(len(self.target_indices) * self.total_origins)

        self.lags = np.arange(self.window + 1)
        self.msd_x = self._msd_accumulator[:, 0] / norm
        self.msd_y = self._msd_accumulator[:, 1] / norm
        self.msd_z = self._msd_accumulator[:, 2] / norm
        self.msd_tot = self.msd_x + self.msd_y + self.msd_z

        if self.time_step is not None:
            self.fit_results = self._fit_diffusion(
                self.lags,
                self.msd_x,
                self.msd_y,
                self.msd_z,
                self.msd_tot,
                self.time_step,
                self.fit_window
            )

        return (
            self.lags,
            self.msd_x,
            self.msd_y,
            self.msd_z,
            self.msd_tot
        )

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms of the MSD analysis."""
        return self.topology.n_atoms

    @staticmethod
    def _unwrap_shift(
        displacement: Np2DNumberArray,
        cell: Cell
    ) -> Np2DNumberArray:
        """
        Calculates the change of the unwrapping shift vectors
        for the given per-frame displacements.

        The displacement vectors are folded back into the
        minimum image convention with respect to the given cell
        by subtracting the box matrix multiplied with the
        rounded (round-half-even) fractional displacement. For
        orthorhombic cells this reduces exactly to the legacy
        Diffcalc expression -box * rint(displacement / box).

        Parameters
        ----------
        displacement : Np2DNumberArray
            The per-frame displacement vectors of shape (n, 3).
        cell : Cell
            The unit cell of the current frame.

        Returns
        -------
        Np2DNumberArray
            The change of the unwrapping shift vectors of
            shape (n, 3).
        """

        if cell.is_vacuum:
            return np.zeros_like(displacement)

        fractional = displacement @ cell.inverse_box_matrix.T

        return -np.rint(fractional) @ cell.box_matrix.T

    @classmethod
    def _fit_diffusion(
        cls,
        lags: Np1DNumberArray,
        msd_x: Np1DNumberArray,
        msd_y: Np1DNumberArray,
        msd_z: Np1DNumberArray,
        msd_tot: Np1DNumberArray,
        time_step: PositiveReal,
        fit_window: PositiveInt
    ) -> Dict[str, MSDDiffusionFit]:
        """
        Fits the trailing part of the MSD components linearly
        and calculates diffusion coefficients.

        The last fit_window points of every MSD component are
        fitted with scipy.stats.linregress. The diffusion
        coefficient is calculated from the slope via the
        Einstein relation D = slope / (2 * dim) with dim = 1
        for the per-axis components and dim = 3 for the total
        MSD. The slopes (Angstrom^2/ps) are converted to m^2/s.

        Parameters
        ----------
        lags : Np1DNumberArray
            The lag indices (in frames) of the MSD analysis.
        msd_x : Np1DNumberArray
            The x-component of the MSD in Angstrom^2.
        msd_y : Np1DNumberArray
            The y-component of the MSD in Angstrom^2.
        msd_z : Np1DNumberArray
            The z-component of the MSD in Angstrom^2.
        msd_tot : Np1DNumberArray
            The total MSD in Angstrom^2.
        time_step : PositiveReal
            The time step between two frames in ps.
        fit_window : PositiveInt
            The number of trailing MSD points used for the fit.

        Returns
        -------
        Dict[str, MSDDiffusionFit]
            The fit results for the keys "x", "y", "z" and
            "total".
        """

        times = lags * time_step

        results = {}

        for label, series, dimension in (
            ("x", msd_x, 1),
            ("y", msd_y, 1),
            ("z", msd_z, 1),
            ("total", msd_tot, 3),
        ):
            fit = linregress(times[-fit_window:], series[-fit_window:])

            factor = ANGSTROM2_PER_PS_TO_M2_PER_S / (2.0 * dimension)

            results[label] = MSDDiffusionFit(
                label=label,
                slope=float(fit.slope),
                slope_stderr=float(fit.stderr),
                intercept=float(fit.intercept),
                r_squared=float(fit.rvalue)**2,
                diffusion_coefficient=float(fit.slope) * factor,
                diffusion_coefficient_stderr=float(fit.stderr) * factor,
            )

        return results
