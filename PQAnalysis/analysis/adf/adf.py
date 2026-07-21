"""
A module containing the ADF class. The ADF class is used to calculate
the angular distribution function of a reference (center) selection
and two ligand selections. The angular distribution function (ADF) is
a measure of the probability density of finding the ``j-i-k`` angle
(spanned at a center atom ``i`` by two neighbouring atoms ``j`` and
``k``) at a given value.
"""

import itertools
import logging

# 3rd party imports
import numpy as np

# 3rd party imports
from beartype.typing import List, Tuple
from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis import config
from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import Cell, Cells
from PQAnalysis.traj import Trajectory, TrajectoryFormat
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import TrajectoryReader, RawTrajectoryReader
from PQAnalysis.io.traj_file.exceptions import TrajectoryReaderError
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import ADFError
from . import _adf_kernel_py

try:
    from ._adf_kernel import adf_frame_histogram  # pylint: disable=import-error
except ModuleNotFoundError:
    from ._adf_kernel_py import adf_frame_histogram



class ADF:  # pylint: disable=too-many-instance-attributes

    """
    A class for calculating the angular distribution function (ADF) of
    a reference (center) selection and two ligand selections.

    The ADF is a histogram of the ``j-i-k`` angle at the apex ``i``,
    where ``i`` is a center atom of the reference selection, ``j`` is
    an atom of the (ligand-1) target selection and ``k`` is an atom of
    the (ligand-2) second target selection. For every ordered ligand
    pair ``(j, k)`` (with ``j != i``, ``k != i`` and ``k != j``) the
    minimum image vectors ``v_ij`` and ``v_ik`` are computed and the
    angle ``degrees(acos(v_ij . v_ik / (r_ij * r_ik)))`` is binned.
    Optional radial gates restrict the ``i-j`` and ``i-k`` distances of
    the counted triplets.

    The ADF analysis is initialized with the provided parameters and is
    run by calling the :py:meth:`run` method, which returns the middle
    points of the angle bins, the normalized ADF (a probability density
    whose integral over the angle range is one), the raw angle counts
    and the sine-corrected ADF (the plain ADF divided by the solid
    angle factor ``sin(theta)``, again normalized to unit integral).

    The ADF class can be initialized with either a trajectory object or
    a TrajectoryReader object. If a trajectory object is given, it is
    assumed to have a constant topology over all frames! The main
    difference between the two is that the TrajectoryReader object
    allows for lazy loading of the trajectory, meaning that the
    trajectory is only loaded frame by frame when needed. This can be
    useful for large trajectories that do not fit into memory.

    For xyz trajectory files the frames are streamed via the raw-frame
    fast path
    (:py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`)
    and accumulated with a compiled angle-histogram kernel, which is
    considerably faster and produces identical results.

    This is a port of the legacy ``thh_tools`` ADF analysis. Two
    documented deviations from the legacy tool make the calculation
    numerically robust without changing the physics: the acos argument
    is clamped to ``[-1, 1]`` before ``acos`` (the legacy tool did not,
    a latent NaN bug for nearly (anti-)parallel vectors) and an exact
    ``180`` degree collinear triplet, which lands in bin ``n_bins``, is
    discarded instead of written out of bounds. Ligand pairs are
    counted as ordered pairs ``(j, k)``, exactly as the legacy tool
    did; when the second target selection defaults to the first, every
    unordered pair is therefore counted twice (the ``j-i-k`` and the
    ``k-i-j`` orientation), which leaves the shape of the normalized ADF
    unchanged.
    """

    _use_full_atom_default = False
    _n_angle_bins_default = 180
    _max_angle = 180.0

    #: The chunk size (in bytes) of the buffered header-only cell scan.
    _CELL_SCAN_CHUNK_SIZE = 16 * 1024 * 1024
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        traj: Trajectory | TrajectoryReader,
        reference_species: SelectionCompatible,
        target_species: SelectionCompatible,
        target_species_2: SelectionCompatible | None = None,
        use_full_atom_info: bool | None = False,
        n_angle_bins: PositiveInt | None = None,
        delta_angle: PositiveReal | None = None,
        r_min_1: PositiveReal | None = None,
        r_max_1: PositiveReal | None = None,
        r_min_2: PositiveReal | None = None,
        r_max_2: PositiveReal | None = None,
    ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory to analyze. If a TrajectoryReader is
            provided, the trajectory is read frame by frame via a
            frame_generator.
        reference_species : SelectionCompatible
            The reference (center ``i``) species of the ADF analysis.
        target_species : SelectionCompatible
            The (ligand-1, ``j``) target species of the ADF analysis.
        target_species_2 : SelectionCompatible | None, optional
            The (ligand-2, ``k``) second target species of the ADF
            analysis. If None (default) the first target species is
            reused, so that the ``j-i-k`` angle is taken within one
            ligand set.
        use_full_atom_info : bool | None, optional
            Whether to use the full atom information of the trajectory
            or not, by default None (False).
        n_angle_bins : PositiveInt | None, optional
            The number of angle bins spanning ``0`` to ``180`` degrees,
            by default None (180).
        delta_angle : PositiveReal | None, optional
            The width (in degrees) of the angle bins, by default None.
            Mutually exclusive with n_angle_bins.
        r_min_1 : PositiveReal | None, optional
            The lower ``i-j`` gate radius (inclusive), by default None
            (no lower bound). Enabling either ``i-j`` gate bound
            activates the ``i-j`` radial gate.
        r_max_1 : PositiveReal | None, optional
            The upper ``i-j`` gate radius (exclusive), by default None
            (no upper bound).
        r_min_2 : PositiveReal | None, optional
            The lower ``i-k`` gate radius (inclusive), by default None
            (no lower bound). Enabling either ``i-k`` gate bound
            activates the ``i-k`` radial gate.
        r_max_2 : PositiveReal | None, optional
            The upper ``i-k`` gate radius (exclusive), by default None
            (no upper bound).

        Raises
        ------
        ADFError
            If the trajectory is empty.
        ADFError
            If n_angle_bins and delta_angle are both specified. This
            would lead to ambiguous results.

        Notes
        -----
        To set up the angle bins at least one of n_angle_bins or
        delta_angle can be specified, but not both. If neither is
        given, n_angle_bins defaults to 180 (one degree wide bins). If
        only n_angle_bins is given, delta_angle is computed as
        ``180 / n_angle_bins``. If only delta_angle is given,
        n_angle_bins is computed as ``int(180 / delta_angle)`` and the
        covered angle range is ``n_angle_bins * delta_angle``.

        See Also
        --------
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        :py:class:`~PQAnalysis.topology.selection.Selection`
        :py:class:`~PQAnalysis.io.trajectoryReader.TrajectoryReader`
        """

        ##############
        # dummy init #
        ##############

        self.angle_bins = np.array([])
        self.counts = np.array([])
        self.normalized_angle_bins = np.array([])
        self.sin_corrected_angle_bins = np.array([])
        self._reference_indices = np.array([], dtype=np.int64)
        self._target_indices = np.array([], dtype=np.int64)
        self._target_indices_2 = np.array([], dtype=np.int64)

        #####################################################
        # Initialize parameters with default values if None #
        #####################################################

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info

        self._setup_radial_gates(r_min_1, r_max_1, r_min_2, r_max_2)

        ################################
        # Initialize Selection objects #
        ################################

        self.reference_species = reference_species
        self.target_species = target_species

        if target_species_2 is None:
            self.target_species_2 = target_species
        else:
            self.target_species_2 = target_species_2

        self.reference_selection = Selection(reference_species)
        self.target_selection = Selection(target_species)
        self.target_selection_2 = Selection(self.target_species_2)

        ############################################
        # Initialize Trajectory iterator/generator #
        ############################################

        self._raw_reader = None
        self.frame_generator = None

        if self._use_raw_fast_path(traj):
            # fast path: lazy loading of the raw frame data from
            # file(s) without per-frame AtomicSystem construction; the
            # cells are collected with a cheap header-only scan
            self._raw_reader = RawTrajectoryReader(
                traj.filenames,
                traj_format=traj.traj_format,
                md_format=traj.md_format,
            )
            self.cells, self._setup_cells = self._scan_cells(traj.filenames)
            self.first_frame = self._raw_reader.read_first_frame()
        else:
            self.cells = traj.cells
            self._setup_cells = self.cells

            if isinstance(traj, TrajectoryReader):
                # lazy loading of trajectory from file(s)
                self.frame_generator = traj.frame_generator()
            elif len(traj) > 0:
                # use trajectory object as iterator
                self.frame_generator = iter(traj)
            else:
                self.logger.error(
                    "Trajectory cannot be of length 0.",
                    exception=ADFError
                )

            self.first_frame = next(self.frame_generator)

        if traj.topology is not None:
            self.topology = traj.topology
        else:
            self.topology = self.first_frame.topology

        self._setup_angle_bins(n_angle_bins=n_angle_bins, delta_angle=delta_angle)

        self.reference_indices = self.reference_selection.select(
            self.topology,
            self.use_full_atom_info
        )

        self.target_indices = self.target_selection.select(
            self.topology,
            self.use_full_atom_info
        )

        self.target_indices_2 = self.target_selection_2.select(
            self.topology,
            self.use_full_atom_info
        )

    def _use_raw_fast_path(self, traj: Trajectory | TrajectoryReader) -> bool:
        """
        Whether the raw-frame fast path is used for the given
        trajectory input.

        The fast path streams the raw per-frame coordinates via
        :py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
        and accumulates the angle histogram with a compiled kernel,
        producing bit-identical results. It is only used for a
        TrajectoryReader input with an xyz trajectory format. All other
        inputs take the original path.

        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory input of the ADF analysis.

        Returns
        -------
        bool
            True if the raw-frame fast path is used.
        """

        return (
            isinstance(traj, TrajectoryReader)
            and traj.traj_format == TrajectoryFormat.XYZ
        )

    def _setup_radial_gates(
        self,
        r_min_1: PositiveReal | None,
        r_max_1: PositiveReal | None,
        r_min_2: PositiveReal | None,
        r_max_2: PositiveReal | None,
    ):
        """
        Sets up the optional ``i-j`` and ``i-k`` radial gates.

        A gate is activated if either of its bounds is given. Omitted
        bounds default to ``0`` (lower) and infinity (upper), so that
        the gate only restricts the specified side.

        Parameters
        ----------
        r_min_1 : PositiveReal | None
            The lower ``i-j`` gate radius (inclusive) or None.
        r_max_1 : PositiveReal | None
            The upper ``i-j`` gate radius (exclusive) or None.
        r_min_2 : PositiveReal | None
            The lower ``i-k`` gate radius (inclusive) or None.
        r_max_2 : PositiveReal | None
            The upper ``i-k`` gate radius (exclusive) or None.
        """

        self.r_min_1 = r_min_1
        self.r_max_1 = r_max_1
        self.r_min_2 = r_min_2
        self.r_max_2 = r_max_2

        self._gate_1 = r_min_1 is not None or r_max_1 is not None
        self._gate_2 = r_min_2 is not None or r_max_2 is not None

        self._r_min_1_eff = r_min_1 if r_min_1 is not None else 0.0
        self._r_max_1_eff = r_max_1 if r_max_1 is not None else np.inf
        self._r_min_2_eff = r_min_2 if r_min_2 is not None else 0.0
        self._r_max_2_eff = r_max_2 if r_max_2 is not None else np.inf

    @classmethod
    def _scan_cells(cls, filenames: List[str]) -> Tuple[Cells, Cells]:
        """
        Collects the cells of the trajectory with a header-only scan.

        This is the fast-path counterpart of the cells full scan of
        :py:attr:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.cells`
        with identical semantics: the number of atoms is taken from the
        first line of the first file and determines the frame stride of
        all files, frames without box information inherit the cell of
        the last frame that had one (also across file boundaries) and
        invalid box headers raise the same error. Instead of building
        one Cell object per frame, textually identical box headers share
        a single (deduplicated) Cell object, which makes the scan
        considerably cheaper for constant-box trajectories.

        Parameters
        ----------
        filenames : List[str]
            The names of the trajectory files to scan.

        Returns
        -------
        cells : Cells
            The cells of all frames of the trajectory (with shared Cell
            objects for textually identical boxes).
        unique_cells : Cells
            The unique cells of the trajectory in order of first
            appearance.

        Raises
        ------
        TrajectoryReaderError
            If the box information of a header line is invalid.
        """

        with open(filenames[0], "r", encoding="utf-8") as file:
            n_atoms = int(file.readline().split()[0])

        # +2 for the cell/atom_count + comment lines
        stride = n_atoms + 2

        # mutable scan state shared across files (deduplicated cells,
        # the running list and the inherited last/vacuum cells)
        state = {
            "cells": [],
            "unique_cells": [],
            "cell_cache": {},
            "last_cell": None,
            "vacuum_cell": None,
        }

        for filename in filenames:
            cls._scan_file_cells(filename, stride, state)

        return state["cells"], state["unique_cells"]

    @classmethod
    def _scan_file_cells(cls, filename: str, stride: int, state: dict):
        """
        Scans the cell headers of a single trajectory file.

        Reads the file in large byte chunks and resolves the box header
        of every frame (each ``stride`` lines) into a Cell object,
        appending it to ``state["cells"]``. The scan state is updated in
        place so that the inherited last/vacuum cell carries across file
        boundaries.

        Parameters
        ----------
        filename : str
            The trajectory file to scan.
        stride : int
            The number of lines per frame (atoms + 2 header lines).
        state : dict
            The mutable scan state (see :py:meth:`_scan_cells`).
        """

        line_number = 0
        offset = 0
        pending = b""

        with open(filename, "rb") as file:
            while True:
                chunk = file.read(cls._CELL_SCAN_CHUNK_SIZE)
                at_eof = chunk == b""

                lines = (pending + chunk).split(b"\n")

                if at_eof:
                    # a trailing line without a final newline counts
                    # as a line
                    if lines[-1] == b"":
                        lines.pop()
                else:
                    pending = lines.pop()

                index = offset

                while index < len(lines):
                    line_number += 1

                    stripped_line = lines[index].decode("utf-8").strip()

                    cell = cls._resolve_header_cell(
                        stripped_line, filename, line_number, state
                    )

                    state["cells"].append(cell)
                    state["last_cell"] = cell

                    index += stride

                offset = index - len(lines)

                if at_eof:
                    break

    @classmethod
    def _resolve_header_cell(
        cls,
        stripped_line: str,
        filename: str,
        line_number: int,
        state: dict,
    ) -> Cell:
        """
        Resolves one frame header line into a (deduplicated) Cell.

        A single-token header (no box information) inherits the last
        cell or a shared vacuum cell; a 4- or 7-token header is turned
        into a Cell (textually identical boxes share one object); any
        other header raises. The unique-cell and cache state is updated
        in place.

        Parameters
        ----------
        stripped_line : str
            The stripped header line of the frame.
        filename : str
            The trajectory file the line belongs to.
        line_number : int
            The 1-based line number of the header line.
        state : dict
            The mutable scan state (see :py:meth:`_scan_cells`).

        Returns
        -------
        Cell
            The resolved cell of the frame.

        Raises
        ------
        TrajectoryReaderError
            If the box information of the header line is invalid.
        """

        splitted_line = stripped_line.split()

        if len(splitted_line) == 1:
            if state["last_cell"] is not None:
                return state["last_cell"]

            if state["vacuum_cell"] is None:
                state["vacuum_cell"] = Cell()
                state["unique_cells"].append(state["vacuum_cell"])

            return state["vacuum_cell"]

        if len(splitted_line) in {4, 7}:
            key = tuple(splitted_line[1:])
            cell = state["cell_cache"].get(key)

            if cell is None:
                cell = Cell(*(float(value) for value in splitted_line[1:]))
                state["cell_cache"][key] = cell
                state["unique_cells"].append(cell)

            return cell

        cls.logger.error(
            (
                "Invalid number of arguments for box:"
                f" {len(splitted_line)} encountered in"
                f" file {filename}:{line_number}"
                f" = {stripped_line}"
            ),
            exception=TrajectoryReaderError,
        )

        return None

    def _setup_angle_bins(
        self,
        n_angle_bins: PositiveInt | None = None,
        delta_angle: PositiveReal | None = None,
    ):
        """
        Sets up the angle bins of the ADF analysis.

        The angle range is fixed to ``0`` to ``180`` degrees. If both
        n_angle_bins and delta_angle are given an error is raised. If
        neither is given, n_angle_bins defaults to 180. If only
        n_angle_bins is given, delta_angle is ``180 / n_angle_bins``. If
        only delta_angle is given, n_angle_bins is ``int(180 /
        delta_angle)``.

        Parameters
        ----------
        n_angle_bins : PositiveInt | None, optional
            The number of angle bins, by default None.
        delta_angle : PositiveReal | None, optional
            The width (in degrees) of the angle bins, by default None.

        Raises
        ------
        ADFError
            If n_angle_bins and delta_angle are both specified.
        """

        if n_angle_bins is not None and delta_angle is not None:
            self.logger.error(
                (
                "It is not possible to specify both n_angle_bins and "
                "delta_angle in the same ADF analysis as this would "
                "lead to ambiguous results."
                ),
                exception=ADFError
            )

        if n_angle_bins is None and delta_angle is None:
            self.n_angle_bins = self._n_angle_bins_default
            self.delta_angle = self._max_angle / self.n_angle_bins
        elif n_angle_bins is not None:
            self.n_angle_bins = n_angle_bins
            self.delta_angle = self._max_angle / self.n_angle_bins
        else:
            self.n_angle_bins = int(self._max_angle / delta_angle)
            self.delta_angle = delta_angle

            if self.n_angle_bins < 1:
                self.logger.error(
                    (
                    "The provided delta_angle is larger than the angle "
                    "range of 180 degrees, resulting in no angle bins."
                    ),
                    exception=ADFError
                )

        self.max_angle = self.n_angle_bins * self.delta_angle

        self.angle_bin_middle_points = self._setup_bin_middle_points(
            self.n_angle_bins,
            self.delta_angle
        )

        self.angle_bins = np.zeros(self.n_angle_bins)

    @timeit_in_class
    def run(
        self
    ) -> Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]:
        """
        Runs the ADF analysis.

        This method runs the ADF analysis and returns the middle points
        of the angle bins, the normalized ADF (a probability density
        whose integral over the angle range is one), the raw angle
        counts and the sine-corrected ADF.

        This method will display a progress bar by default. This can be
        disabled by setting with_progress_bar to False.

        Returns
        -------
        angle_bin_middle_points : Np1DNumberArray
            The middle points of the angle bins in degrees.
        normalized_angle_bins : Np1DNumberArray
            The normalized ADF (integral over the angle range is one).
        counts : Np1DNumberArray
            The raw angle counts of the ADF analysis.
        sin_corrected_angle_bins : Np1DNumberArray
            The sine-corrected ADF (integral over the angle range is
            one).
        """

        self._initialize_run()

        if self._raw_reader is not None:
            self._calculate_angles_raw()
        else:
            self._calculate_angles()

        return self._finalize_run()

    def _initialize_run(self):
        """
        Initializes the ADF analysis for running.

        This method is called by the run method of the ADF class. It
        resets the angle histogram accumulator and materializes the
        contiguous int64 index arrays consumed by the kernels.
        """

        self.angle_bins = np.zeros(self.n_angle_bins)

        self._reference_indices = np.ascontiguousarray(
            self.reference_indices, dtype=np.int64
        )
        self._target_indices = np.ascontiguousarray(
            self.target_indices, dtype=np.int64
        )
        self._target_indices_2 = np.ascontiguousarray(
            self.target_indices_2, dtype=np.int64
        )

    @staticmethod
    def _cell_box_arrays(
        cell: Cell
    ) -> Tuple[Np1DNumberArray, np.ndarray, np.ndarray, int]:
        """
        Extracts the loop-invariant box data of a cell.

        Parameters
        ----------
        cell : Cell
            The unit cell of a trajectory frame.

        Returns
        -------
        box_lengths : Np1DNumberArray
            The box lengths of the cell (float64).
        box : np.ndarray
            The box matrix of the cell (float64, C-contiguous).
        inv_box : np.ndarray
            The inverse box matrix of the cell (float64, C-contiguous).
        is_orthorhombic : int
            Whether all box angles of the cell are exactly 90 degrees.
        """

        is_orthorhombic = 1 if (
            cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90
        ) else 0

        box_lengths = np.ascontiguousarray(cell.box_lengths, dtype=np.float64)
        box = np.ascontiguousarray(cell.box_matrix, dtype=np.float64)
        inv_box = np.ascontiguousarray(
            cell.inverse_box_matrix, dtype=np.float64
        )

        return box_lengths, box, inv_box, is_orthorhombic

    def _calculate_angles(self):
        """
        Calculates the angle histogram of the ADF analysis.

        This method is called by the run method of the ADF class when
        the original (non fast-path) input is used. It iterates over the
        frames of the trajectory and accumulates the per-frame angle
        histogram with the numpy fallback kernel (which works for any
        position dtype). The result is bit-identical to the fast path
        for float32 trajectories.
        """

        hist = np.zeros(self.n_angle_bins, dtype=np.int64)

        for frame in tqdm(
            itertools.chain([self.first_frame], self.frame_generator),
            total=self.n_frames,
            disable=not config.with_progress_bar):

            box_lengths, box, inv_box, is_orthorhombic = self._cell_box_arrays(
                frame.cell
            )

            _adf_kernel_py.adf_frame_histogram(
                frame.pos,
                self._reference_indices,
                self._target_indices,
                self._target_indices_2,
                box_lengths,
                box,
                inv_box,
                is_orthorhombic,
                1 if self._gate_1 else 0,
                self._r_min_1_eff,
                self._r_max_1_eff,
                1 if self._gate_2 else 0,
                self._r_min_2_eff,
                self._r_max_2_eff,
                self.delta_angle,
                self.n_angle_bins,
                hist,
            )

        self.angle_bins += hist

    def _calculate_angles_raw(self):
        """
        Calculates the angle histogram of the ADF analysis using the
        raw-frame fast path.

        This method is the fast-path counterpart of
        :py:meth:`_calculate_angles` used when the trajectory is read
        from xyz file(s): it streams the raw per-frame coordinates via
        :py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
        (no per-frame AtomicSystem construction) and accumulates the
        angle histogram with the ADF frame kernel. The box data of the
        cell are loop invariants that are only recomputed when the Cell
        object yielded by the reader changes identity.

        Raises
        ------
        ADFError
            If a frame does not provide positions for all atoms
            referenced by the selections.
        """

        hist = np.zeros(self.n_angle_bins, dtype=np.int64)

        max_index = int(
            max(
                self._reference_indices.max(initial=-1),
                self._target_indices.max(initial=-1),
                self._target_indices_2.max(initial=-1),
            )
        )

        # loop-invariant cell data, recomputed only when the yielded
        # Cell object changes identity
        last_cell = None
        box_lengths = np.ones(3)
        box = np.eye(3)
        inv_box = np.eye(3)
        is_orthorhombic = 1

        counter = 0

        for values, cell in tqdm(
            self._raw_reader.raw_frame_generator(),
            total=self.n_frames,
            disable=not config.with_progress_bar):

            counter += 1

            if values.shape[0] <= max_index:
                self.logger.error(
                    (
                    f"Frame {counter} of the trajectory provides "
                    f"only {values.shape[0]} atoms, but the "
                    "selections reference the atom index "
                    f"{max_index}. Please provide a trajectory "
                    "with a consistent number of atoms."
                    ),
                    exception=ADFError
                )

            if cell is not last_cell:
                last_cell = cell
                box_lengths, box, inv_box, is_orthorhombic = (
                    self._cell_box_arrays(cell)
                )

            adf_frame_histogram(
                values,
                self._reference_indices,
                self._target_indices,
                self._target_indices_2,
                box_lengths,
                box,
                inv_box,
                is_orthorhombic,
                1 if self._gate_1 else 0,
                self._r_min_1_eff,
                self._r_max_1_eff,
                1 if self._gate_2 else 0,
                self._r_min_2_eff,
                self._r_max_2_eff,
                self.delta_angle,
                self.n_angle_bins,
                hist,
            )

        self.angle_bins += hist

    def _finalize_run(
        self
    ) -> Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]:
        """
        Finalizes the ADF analysis after running.

        This method is called by the run method of the ADF class. It
        normalizes the raw angle counts into a probability density whose
        integral over the angle range is one and additionally computes
        the sine-corrected density (the plain ADF divided by the solid
        angle factor ``sin(theta)`` at the bin center), again normalized
        to unit integral.

        Returns
        -------
        angle_bin_middle_points : Np1DNumberArray
            The middle points of the angle bins in degrees.
        normalized_angle_bins : Np1DNumberArray
            The normalized ADF (integral over the angle range is one).
        counts : Np1DNumberArray
            The raw angle counts of the ADF analysis.
        sin_corrected_angle_bins : Np1DNumberArray
            The sine-corrected ADF (integral over the angle range is
            one).
        """

        self.counts = self.angle_bins.astype(np.float64)

        self.normalized_angle_bins = self._normalize(
            self.counts,
            self.delta_angle
        )

        self.sin_corrected_angle_bins = self._sin_correct(
            self.counts,
            self.angle_bin_middle_points,
            self.delta_angle
        )

        return (
            self.angle_bin_middle_points,
            self.normalized_angle_bins,
            self.counts,
            self.sin_corrected_angle_bins
        )

    @property
    def n_frames(self) -> int:
        """int: The number of frames of the ADF analysis."""
        return len(self.cells)

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms of the ADF analysis."""
        return self.topology.n_atoms

    @classmethod
    def _setup_bin_middle_points(
        cls,
        n_angle_bins: PositiveInt,
        delta_angle: PositiveReal
    ) -> Np1DNumberArray:
        """
        Sets up the middle points of the angle bins of the ADF
        analysis for outputting the ADF analysis.

        Parameters
        ----------
        n_angle_bins : PositiveInt
            The number of angle bins.
        delta_angle : PositiveReal
            The width (in degrees) of the angle bins.

        Returns
        -------
        Np1DNumberArray
            The middle points of the angle bins in degrees.
        """

        return (np.arange(n_angle_bins) + 0.5) * delta_angle

    @classmethod
    def _normalize(
        cls,
        counts: Np1DNumberArray,
        delta_angle: PositiveReal
    ) -> Np1DNumberArray:
        """
        Normalizes the raw angle counts into a probability density.

        The density is normalized such that its integral over the angle
        range (the sum of the density times ``delta_angle``) equals one.

        Parameters
        ----------
        counts : Np1DNumberArray
            The raw angle counts of the ADF analysis.
        delta_angle : PositiveReal
            The width (in degrees) of the angle bins.

        Returns
        -------
        Np1DNumberArray
            The normalized ADF.
        """

        total = counts.sum()

        if total <= 0:
            return np.zeros_like(counts)

        return counts / (total * delta_angle)

    @classmethod
    def _sin_correct(
        cls,
        counts: Np1DNumberArray,
        bin_middle_points: Np1DNumberArray,
        delta_angle: PositiveReal
    ) -> Np1DNumberArray:
        """
        Computes the sine-corrected angular density.

        Dividing the raw counts by ``sin(theta)`` at the bin center
        removes the solid-angle factor that biases an isotropic
        distribution towards ``90`` degrees, yielding the true angular
        density. The result is normalized to unit integral over the
        angle range.

        Parameters
        ----------
        counts : Np1DNumberArray
            The raw angle counts of the ADF analysis.
        bin_middle_points : Np1DNumberArray
            The middle points of the angle bins in degrees.
        delta_angle : PositiveReal
            The width (in degrees) of the angle bins.

        Returns
        -------
        Np1DNumberArray
            The sine-corrected, unit-integral ADF.
        """

        sine = np.sin(np.radians(bin_middle_points))

        weights = np.divide(
            counts,
            sine,
            out=np.zeros_like(counts),
            where=sine != 0
        )

        total = weights.sum()

        if total <= 0 or not np.isfinite(total):
            return np.zeros_like(counts)

        return weights / (total * delta_angle)
