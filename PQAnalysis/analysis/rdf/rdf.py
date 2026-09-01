"""
A module containing the RDF class. The RDF class is used
to calculate the radial distribution of a reference
selection to a target selection. The radial distribution
function (RDF) is a measure of the probability density 
of finding a particle at a distance r from another particle. 
"""

import itertools
import logging
import math

from concurrent.futures import ThreadPoolExecutor
from os import cpu_count

# 3rd party imports
import numpy as np

# 3rd party imports
from beartype.typing import List, Tuple
from PQAnalysis.utils.progress import tqdm

# local absolute imports
from PQAnalysis import config
from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import distance, Cell, Cells
from PQAnalysis.traj import (
    Trajectory,
    TrajectoryFormat,
    check_trajectory_pbc,
    check_trajectory_vacuum,
)
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import TrajectoryReader, RawTrajectoryReader
from PQAnalysis.io.traj_file.exceptions import TrajectoryReaderError
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import RDFError

try:
    from ._rdf_kernel import (  # pylint: disable=import-error
        legacy_rdf_batch_histogram,
        legacy_rdf_frame_histogram,
        rdf_frame_histogram,
    )
except ModuleNotFoundError:
    from ._rdf_kernel_py import (
        legacy_rdf_batch_histogram,
        legacy_rdf_frame_histogram,
        rdf_frame_histogram,
    )



# The flat argument list keeps the threaded worker setup compact.
# pylint: disable-next=too-many-arguments,too-many-positional-arguments,too-many-locals
def _raw_histogram_range(
    values,
    cells,
    start,
    stop,
    reference_indices,
    target_indices,
    r_min,
    delta_r,
    n_bins,
    use_legacy,
    legacy_box_lengths,
):
    """Calculate a private exact histogram over a frame range."""
    hist = np.zeros(n_bins, dtype=np.int64)

    if use_legacy:
        legacy_rdf_batch_histogram(
            values[start:stop],
            reference_indices,
            target_indices,
            legacy_box_lengths[start:stop],
            delta_r,
            n_bins,
            hist,
        )

        return hist

    last_cell = None
    box_lengths = np.ones(3)
    box = np.eye(3)
    inv_box = np.eye(3)
    is_orthorhombic = 1

    for frame_index in range(start, stop):
        cell = cells[frame_index]

        if cell is not last_cell:
            last_cell = cell
            is_orthorhombic = int(
                cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90
            )
            box_lengths = np.ascontiguousarray(
                cell.box_lengths,
                dtype=np.float64,
            )
            box = np.ascontiguousarray(
                cell.box_matrix,
                dtype=np.float64,
            )
            inv_box = np.ascontiguousarray(
                cell.inverse_box_matrix,
                dtype=np.float64,
            )

        rdf_frame_histogram(
            values[frame_index],
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
        )

    return hist



class RDF:

    """
    A class for calculating the radial distribution of a 
    reference selection to a target selection. The radial
    distribution function (RDF) is a measure of the 
    probability density of finding a particle at a 
    distance r from another particle. 

    The RDF analysis is initialized with the provided 
    parameters. The RDF analysis can be run by calling
    the run method. The run method returns the middle 
    points of the bins of the RDF analysis, the normalized
    bins of the RDF analysis based on the spherical shell
    model, the integrated bins of the RDF analysis, the
    normalized bins of the RDF analysis based on the
    number of atoms in the system and the differential 
    bins of the RDF analysis based on the spherical 
    shell model.

    The RDF class can be initialized with either a 
    trajectory object or via a TrajectoryReader object.
    If a trajectory object is given, it is assumed to 
    have a constant topology over all frames! The main 
    difference between the two is that the
    TrajectoryReader object allows for lazy loading of
    the trajectory, meaning that the trajectory is only
    loaded frame by frame when needed. This can be useful
    for large trajectories that do not fit into memory.

    For xyz trajectory files without intra-molecular exclusion the
    frames are streamed in float64 via the raw-frame
    fast path
    (:py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`)
    and accumulated with a compiled distance-histogram
    kernel. A periodic orthorhombic calculation specified by
    ``delta_r`` alone, with the default ``r_min = 0``, reproduces the
    numeric operation order of the legacy ``thh_tools/RDF`` program.
    Other geometries and explicit radial ranges use the general
    PQAnalysis RDF semantics.
    """

    _use_full_atom_default = False
    _no_intra_molecular_default = False
    _r_min_default = 0.0

    #: The chunk size (in bytes) of the buffered header-only cell scan.
    _CELL_SCAN_CHUNK_SIZE = 16 * 1024 * 1024
    _batch_max_bytes = 512 * 1024 * 1024
    _batch_max_workers = 16
    _batch_pairs_per_worker = 250_000
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(  # pylint: disable=too-many-branches
        self,
        traj: Trajectory | TrajectoryReader,
        reference_species: SelectionCompatible,
        target_species: SelectionCompatible,
        use_full_atom_info: bool | None = False,
        no_intra_molecular: bool | None = False,
        n_bins: PositiveInt | None = None,
        delta_r: PositiveReal | None = None,
        r_max: PositiveReal | None = None,
        r_min: PositiveReal | None = 0.0,
    ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory to analyze. If a TrajectoryReader is provided,
            the trajectory frame by frame via a frame_generator
        reference_species : SelectionCompatible
            The reference species of the RDF analysis.
        target_species : SelectionCompatible
            The target species of the RDF analysis.
        use_full_atom_info : bool | None, optional
            Whether to use the full atom information of the trajectory
            or not, by default None (False).
        no_intra_molecular : bool | None, optional
            Whether to exclude intra-molecular distances or not, by default None (False).
        n_bins : PositiveInt | None, optional
            number of bins, by default None
        delta_r : PositiveReal | None, optional
            delta r between bins, by default None
        r_max : PositiveReal | None, optional
            maximum radius from reference species of the RDF analysis,
            by default None
        r_min : PositiveReal, optional
            minimum (starting) radius from reference species of the 
            RDF analysis, by default 0.0 (equals to None)

        Raises
        ------
        RDFError
            If the trajectory is not fully periodic or fully in vacuum.
            Meaning that some frames are in vacuum and others are periodic.
        RDFError
            If the trajectory is empty.
        RDFError
            If n_bins and delta_r are both not specified.
        RDFError
            If n_bins, delta_r and r_max are all specified.
            This would lead to ambiguous results.

        Notes
        -----        
        Furthermore, to initialize the RDF analysis object at least
        one of the parameters n_bins or delta_r must be specified. 
        If n_bins and delta_r are both specified, r_max is calculated
        from these parameters. If n_bins and r_max are both specified,
        delta_r is calculated from these parameters. If delta_r and 
        r_max are both specified, n_bins is calculated from 
        these parameters.

        It is not possible to specify all of n_bins, delta_r and r_max
        in the same RDF analysis as this would lead to ambiguous results.

        It is also possible to initialize a non-vacuum trajectory by
        only using n_bins or delta_r. In this case, r_max is calculated
        from the provided parameters and the box vectors of the 
        trajectory. If the trajectory is in vacuum, an RDFError is
        raised as the maximum radius cannot be inferred from the 
        box vectors.

        See Also
        --------
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        :py:class:`~PQAnalysis.topology.selection.Selection`
        :py:class:`~PQAnalysis.io.trajectoryReader.TrajectoryReader`
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        """

        ##############
        # dummy init #
        ##############

        self._average_volume = 0.0
        self._reference_density = 0.0
        self.delta_r, self.n_bins, self.r_max = 0.0, 0, 0.0
        self.bin_middle_points, self.bins = np.array([]), np.array([])
        self.target_index_combinations = []
        self.normalized_bins = np.array([])
        self.integrated_bins = np.array([])
        self.normalized_bins2 = np.array([])
        self.differential_bins = np.array([])

        #####################################################
        # Initialize parameters with default values if None #
        #####################################################

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info

        if no_intra_molecular is None:
            self.no_intra_molecular = self._no_intra_molecular_default
        else:
            self.no_intra_molecular = no_intra_molecular

        if r_min is None:
            self.r_min = self._r_min_default
        else:
            self.r_min = r_min

        ################################
        # Initialize Selection objects #
        ################################

        self.reference_species = reference_species
        self.target_species = target_species

        self.reference_selection = Selection(reference_species)
        self.target_selection = Selection(target_species)

        ############################################
        # Initialize Trajectory iterator/generator #
        ############################################

        self._raw_reader = None
        self._raw_batch = None
        self._raw_batch_attempted = False
        self.frame_generator = None

        if self._use_raw_fast_path(traj):
            self._initialize_raw_fast_path(traj)
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
                    "Trajectory cannot be of length 0.", exception=RDFError
                )

            self.first_frame = next(self.frame_generator)
        if traj.topology is not None:
            self.topology = traj.topology
        else:
            self.topology = self.first_frame.topology

        self._setup_analysis_bins(
            n_bins=n_bins,
            delta_r=delta_r,
            r_max=r_max,
            r_min=self.r_min,
        )

        self.reference_indices = self.reference_selection.select(
            self.topology, self.use_full_atom_info
        )

        self.target_indices = self.target_selection.select(
            self.topology, self.use_full_atom_info
        )

    def _use_raw_fast_path(self, traj: Trajectory | TrajectoryReader) -> bool:
        """
        Whether the raw-frame fast path is used for the given
        trajectory input.

        The fast path streams the raw per-frame coordinates via
        :py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
        and accumulates the distance histogram with a compiled
        kernel, producing bit-identical results. It is only used for
        the plain case: a TrajectoryReader input with an xyz
        trajectory format and no intra-molecular exclusion. All
        other inputs take the original path.

        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory input of the RDF analysis.

        Returns
        -------
        bool
            True if the raw-frame fast path is used.
        """

        return (
            isinstance(traj, TrajectoryReader) and
            traj.traj_format == TrajectoryFormat.XYZ and
            not self.no_intra_molecular
        )

    def _initialize_raw_fast_path(self, traj: TrajectoryReader):
        """Initialize raw input and reuse a compatible batch for setup."""
        self._raw_reader = RawTrajectoryReader(
            traj.filenames,
            traj_format=traj.traj_format,
            md_format=traj.md_format,
            dtype="float64",
        )
        self.first_frame = self._raw_reader.read_first_frame()

        if not rdf_frame_histogram.__module__.endswith("_rdf_kernel_py"):
            self._raw_batch_attempted = True
            self._raw_batch = self._raw_reader.try_read_all_frames(
                expected_n_atoms=self.first_frame.n_atoms,
                max_bytes=self._batch_max_bytes,
            )

        if self._raw_batch is None:
            self.cells, self._setup_cells = self._scan_cells(traj.filenames)
            return

        _values, self.cells = self._raw_batch
        self._setup_cells = list({
            id(cell): cell for cell in self.cells
        }.values())

    def _use_legacy_rdf(
        self,
        n_bins: PositiveInt | None,
        delta_r: PositiveReal | None,
        r_max: PositiveReal | None,
    ) -> bool:
        """Whether the input can reproduce the legacy RDF calculation."""
        return (
            self._raw_reader is not None and n_bins is None and
            delta_r is not None and r_max is None and self.r_min == 0.0 and
            all(not cell.is_vacuum for cell in self._setup_cells) and all(
                np.array_equal(cell.box_angles, np.array([90, 90, 90]))
                for cell in self._setup_cells
            )
        )

    def _setup_analysis_bins(
        self,
        n_bins: PositiveInt | None,
        delta_r: PositiveReal | None,
        r_max: PositiveReal | None,
        r_min: PositiveReal,
    ):
        """Selects legacy-compatible or general RDF bin semantics."""
        self._legacy_rdf = self._use_legacy_rdf(
            n_bins=n_bins,
            delta_r=delta_r,
            r_max=r_max,
        )

        if self._legacy_rdf:
            self._setup_legacy_bins(delta_r)
            return

        self._setup_bins(
            n_bins=n_bins,
            delta_r=delta_r,
            r_max=r_max,
            r_min=r_min,
        )

    @classmethod
    def _scan_cells(cls, filenames: List[str]) -> Tuple[Cells, Cells]:
        """
        Collects the cells of the trajectory with a header-only scan.

        This is the fast-path counterpart of the cells full scan of
        :py:attr:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.cells`
        with identical semantics: the number of atoms is taken from
        the first line of the first file and determines the frame
        stride of all files, frames without box information inherit
        the cell of the last frame that had one (also across file
        boundaries) and invalid box headers raise the same error.
        Instead of building one Cell object per frame, textually
        identical box headers share a single (deduplicated) Cell
        object, which makes the scan considerably cheaper for
        constant-box trajectories.

        Parameters
        ----------
        filenames : List[str]
            The names of the trajectory files to scan.

        Returns
        -------
        cells : Cells
            The cells of all frames of the trajectory (with shared
            Cell objects for textually identical boxes).
        unique_cells : Cells
            The unique cells of the trajectory in order of first
            appearance. Every computation of the RDF setup over the
            cells is deduplication invariant (all/any/min over the
            cells), so the unique cells serve as a cheap stand-in
            for the full cell list in the setup checks.

        Raises
        ------
        TrajectoryReaderError
            If the box information of a header line is invalid.
        """

        with open(filenames[0], "r", encoding="utf-8") as file:
            n_atoms = int(file.readline().split()[0])

        # +2 for the cell/atom_count + comment lines
        stride = n_atoms + 2

        cells = []
        unique_cells = []
        cell_cache = {}
        last_cell = None
        vacuum_cell = None

        for filename in filenames:
            line_number = 0
            offset = 0
            pending = b""

            with open(filename, "rb") as file:
                while True:
                    chunk = file.read(cls._CELL_SCAN_CHUNK_SIZE)
                    at_eof = chunk == b""

                    lines = (pending + chunk).split(b"\n")

                    if at_eof:
                        # a trailing line without a final newline
                        # counts as a line
                        if lines[-1] == b"":
                            lines.pop()
                    else:
                        pending = lines.pop()

                    index = offset

                    while index < len(lines):
                        line_number += 1

                        stripped_line = (lines[index].decode("utf-8").strip())
                        splitted_line = stripped_line.split()

                        if len(splitted_line) == 1:

                            if last_cell is not None:
                                cell = last_cell
                            else:
                                if vacuum_cell is None:
                                    vacuum_cell = Cell()
                                    unique_cells.append(vacuum_cell)

                                cell = vacuum_cell

                        elif len(splitted_line) in (4, 7):

                            key = tuple(splitted_line[1:])
                            cell = cell_cache.get(key)

                            if cell is None:
                                cell = Cell(
                                    *(
                                        float(value)
                                        for value in splitted_line[1:]
                                    )
                                )
                                cell_cache[key] = cell
                                unique_cells.append(cell)

                        else:

                            cls.logger.error(
                                (
                                    "Invalid number of arguments for box:"
                                    f" {len(splitted_line)} encountered in"
                                    f" file {filename}:{line_number}"
                                    f" = {stripped_line}"
                                ),
                                exception=TrajectoryReaderError,
                            )

                        cells.append(cell)
                        last_cell = cell

                        index += stride

                    offset = index - len(lines)

                    if at_eof:
                        break

        return cells, unique_cells

    def _setup_legacy_bins(self, delta_r: PositiveReal):
        """Sets bins with the scalar semantics of ``thh_tools/RDF``."""
        self.r_min = 0.0
        self.delta_r = float(np.float32(delta_r))
        self.n_bins = max(
            int(max(cell.box_lengths) / 2.0 / self.delta_r)
            for cell in self.cells
        )
        self.r_max = self.delta_r * self.n_bins
        self.bin_middle_points = np.array(
            [
                (float(np.float32(index)) + 0.5) * self.delta_r
                for index in range(self.n_bins)
            ]
        )
        self.bins = np.zeros(self.n_bins)

    def _setup_bins(
        self,
        n_bins: PositiveInt | None = None,
        delta_r: PositiveReal | None = None,
        r_max: PositiveReal | None = None,
        r_min: PositiveReal = 0.0
    ):
        """
        Sets up the bins of the RDF analysis.

        This method is called by the __init__ method of the 
        RDF class, but can also be called manually to 
        re-initialize the bins of the RDF analysis. It sets
        up the bins of the RDF analysis based on the
        provided parameters. If n_bins and delta_r are both
        specified, r_max is calculated from these parameters.
        If n_bins and r_max are both specified, delta_r is
        calculated from these parameters. If delta_r and 
        r_max are both specified, n_bins is calculated from
        these parameters.

        Parameters
        ----------
        n_bins : PositiveInt | None, optional
            number of bins, by default None
        delta_r : PositiveReal | None, optional
            delta r between bins, by default None
        r_max : PositiveReal | None, optional
            maximum radius of the RDF analysis, by default None
        r_min : PositiveReal, optional
            minimum (starting) radius of the RDF analysis, by default 0.0

        Raises
        ------
        RDFError
            If the trajectory is not fully periodic or fully in vacuum.
            Meaning that some frames are in vacuum and others are periodic.
        RDFError
            If n_bins and delta_r are both not specified.
        RDFError
            If n_bins, delta_r and r_max are all specified. This would lead to ambiguous results.
        """

        self.r_min = r_min

        self._check_trajectory_conditions()

        # check if n_bins and delta_r are both not specified
        if n_bins is None and delta_r is None:
            self.logger.error(
                "Either n_bins or delta_r must be specified.",
                exception=RDFError
            )

        # check if n_bins, delta_r and r_max are all specified
        elif all([n_bins, delta_r, r_max]):
            self.logger.error(
                (
                    "It is not possible to specify all of n_bins, "
                    "delta_r and r_max in the same RDF analysis "
                    "as this would lead to ambiguous results."
                ),
                exception=RDFError
            )

        # set r_max based on the provided parameters n_bins and delta_r
        if n_bins is not None and delta_r is not None:

            self.n_bins = n_bins

            self.delta_r = delta_r

            self.r_max = self._calculate_r_max(
                n_bins, delta_r, r_min, self._setup_cells
            )

            self.n_bins, self.r_max = self._calculate_n_bins(
                delta_r,
                self.r_max,
                r_min
            )

        else:

            self.r_max = r_max if r_max is not None else self._infer_r_max(
                self._setup_cells
            )

            self.r_max = self._check_r_max(self.r_max, self._setup_cells)

            if n_bins is None:

                self.delta_r = delta_r

                self.n_bins, self.r_max = self._calculate_n_bins(
                    delta_r,
                    self.r_max,
                    r_min
                )

            else:
                self.n_bins = n_bins
                self.delta_r = (self.r_max - self.r_min) / self.n_bins

        self.bin_middle_points = self._setup_bin_middle_points(
            self.n_bins, self.r_min, self.r_max, self.delta_r
        )

        self.bins = np.zeros(self.n_bins)

    def _check_trajectory_conditions(self):
        """
        Checks if the trajectory is fully periodic or fully in vacuum.

        Raises
        ------
        RDFError
            If the trajectory is not fully periodic or fully in vacuum.
            Meaning that some frames are in vacuum and others are periodic.
        """

        if not check_trajectory_pbc(
            self._setup_cells
        ) and not check_trajectory_vacuum(self._setup_cells):
            self.logger.error(
                (
                    "The provided trajectory is not fully periodic or "
                    "in vacuum, meaning that some frames are in vacuum "
                    "and others are periodic. This is not supported by "
                    "the RDF analysis."
                ),
                exception=RDFError
            )

    @property
    def average_volume(self) -> PositiveReal:
        """PositiveReal: The average volume of the trajectory."""
        return np.mean([cell.volume for cell in self.cells])

    def _calculate_average_volume(self) -> PositiveReal:
        """
        Calculates the average volume of the trajectory.

        The legacy-compatible path sums orthorhombic ``x*y*z`` volumes
        in frame order, preserving the C implementation's rounding.
        For the other raw-frame paths the volume of every unique cell
        is computed only once and broadcast to the full (shared
        object) cell list before averaging, which is bit-identical
        to (but much cheaper than) the plain mean over the volumes
        of all cells of the :py:attr:`average_volume` property. For
        the original path the property is evaluated as before.

        Returns
        -------
        PositiveReal
            The average volume of the trajectory.
        """

        if self._legacy_rdf:
            volume = 0.0

            for cell in self.cells:
                length_x = float(cell.box_lengths[0])
                length_y = float(cell.box_lengths[1])
                length_z = float(cell.box_lengths[2])
                volume = volume + (length_x * length_y) * length_z

            return volume / float(self.n_frames)

        if self._raw_reader is None:
            return self.average_volume

        volumes = {id(cell): cell.volume for cell in self._setup_cells}

        return np.mean([volumes[id(cell)] for cell in self.cells])

    @timeit_in_class
    def run(
        self
    ) -> Tuple[Np1DNumberArray,
               Np1DNumberArray,
               Np1DNumberArray,
               Np1DNumberArray,
               Np1DNumberArray]:
        """
        Runs the RDF analysis.

        This method runs the RDF analysis and returns the five arrays
        written by :py:class:`~PQAnalysis.analysis.rdf.rdf_output_file_writer.RDFDataWriter`.
        See :ref:`RDF output files <analysis-output-rdf>` for their exact
        definitions and normalization formulas.

        This method will display a progress bar by default.
        This can be disabled by setting with_progress_bar to False.

        Returns
        -------
        bin_middle_points : Np1DNumberArray
            The bin-center distances in Angstrom.
        normalized_bins : Np1DNumberArray
            The dimensionless radial distribution function ``g(r)``.
        integrated_bins : Np1DNumberArray
            The cumulative coordination number.
        normalized_bins2 : Np1DNumberArray
            The density-normalized shell population in Angstrom^3.
        differential_bins : Np1DNumberArray
            The pair-count residual relative to the ideal-gas shell count.
        """

        self._initialize_run()

        if self._raw_reader is not None:
            self._calculate_bins_raw()
        else:
            self._calculate_bins()

        return self._finalize_run()

    def _initialize_run(self):
        """
        Initializes the RDF analysis for running.

        This method is called by the run method of the RDF class.
        It initializes the RDF analysis for running by calculating
        the average volume of the trajectory, the reference 
        density of the RDF analysis and the target index 
        combinations of the RDF analysis.
        """

        self._average_volume = self._calculate_average_volume()

        _ref_indices_len = len(self.reference_indices)

        self._reference_density = _ref_indices_len / self._average_volume

        self._initialize_target_index_combinations()

    def _initialize_target_index_combinations(self):
        """
        Initializes the target index combinations of the RDF analysis.

        This method is called by the _initialize_run method of 
        the RDF class. It initializes the target index combinations
        of the RDF analysis by calculating the target index 
        combinations of the RDF analysis based on the 
        no_intra_molecular parameter.
        """

        self.target_index_combinations = []

        if self.no_intra_molecular:
            for reference_index in self.reference_indices:
                residue_number = self.topology.residue_numbers[reference_index]
                residue_indices = self.topology.residue_atom_indices[
                    residue_number]
                self.target_index_combinations.append(
                    np.setdiff1d(self.target_indices, residue_indices)
                )

    def _calculate_bins(self):
        """
        Calculates the bins of the RDF analysis.

        This method is called by the run method of the RDF class.
        It calculates the bins of the RDF analysis by iterating 
        over the frames of the trajectory and calculating the 
        distances between the reference and target indices of the
        RDF analysis. The bins of the RDF analysis are then 
        calculated from these distances.
        """

        for frame in tqdm(
            itertools.chain([self.first_frame], self.frame_generator),
            total=self.n_frames,
            disable=not config.with_progress_bar
        ):
            for i, reference_index in enumerate(self.reference_indices):

                if self.no_intra_molecular:
                    target_indices = self.target_index_combinations[i]
                else:
                    target_indices = self.target_indices

                target_indices = target_indices[target_indices !=
                                                reference_index]

                reference_position = frame.pos[reference_index]
                target_positions = frame.pos[target_indices]

                distances = distance(
                    reference_position, target_positions, frame.cell
                ).ravel()

                self.bins += self._add_to_bins(
                    distances, self.r_min, self.delta_r, self.n_bins
                )

    def _calculate_bins_raw(self):
        """
        Calculates the bins of the RDF analysis using the raw-frame
        fast path.

        This method is the fast-path counterpart of
        :py:meth:`_calculate_bins` used when the trajectory is read
        from xyz file(s) and no intra-molecular exclusion is active:
        it streams the raw per-frame coordinates via
        :py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
        (no per-frame AtomicSystem construction) and accumulates the
        distance histogram with the RDF frame kernel, which
        replicates the numeric semantics of the original loop bit
        for bit. The box data of the cell are loop invariants that
        are only recomputed when the Cell object yielded by the
        reader changes identity (the raw reader caches Cell objects,
        so constant-box trajectories compute them exactly once).

        Raises
        ------
        RDFError
            If a frame does not provide positions for all atoms
            referenced by the selections.
        """

        if self._calculate_bins_raw_batch():
            return

        hist = np.zeros(self.n_bins, dtype=np.int64)

        reference_indices = np.ascontiguousarray(
            self.reference_indices, dtype=np.int64
        )
        target_indices = np.ascontiguousarray(
            self.target_indices, dtype=np.int64
        )

        max_index = int(
            max(
                reference_indices.max(initial=-1),
                target_indices.max(initial=-1),
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
                    exception=RDFError
                )

            if cell is not last_cell:
                last_cell = cell
                is_orthorhombic = 1 if (
                    cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90
                ) else 0

                box_lengths = np.ascontiguousarray(
                    cell.box_lengths, dtype=np.float64
                )
                box = np.ascontiguousarray(cell.box_matrix, dtype=np.float64)
                inv_box = np.ascontiguousarray(
                    cell.inverse_box_matrix, dtype=np.float64
                )

            if self._legacy_rdf:
                legacy_rdf_frame_histogram(
                    values,
                    reference_indices,
                    target_indices,
                    box_lengths,
                    self.delta_r,
                    self.n_bins,
                    hist,
                )
            else:
                rdf_frame_histogram(
                    values,
                    reference_indices,
                    target_indices,
                    box_lengths,
                    box,
                    inv_box,
                    is_orthorhombic,
                    self.r_min,
                    self.delta_r,
                    self.n_bins,
                    hist,
                )

        self.bins += hist

    # This method owns batch sizing, worker partitioning and exact reduction.
    # pylint: disable-next=too-many-locals
    def _calculate_bins_raw_batch(self) -> bool:
        """Calculate independent frame histograms from one raw batch."""
        if rdf_frame_histogram.__module__.endswith("_rdf_kernel_py"):
            return False

        if self._batch_max_bytes <= 0:
            self._raw_batch = None
            self._raw_batch_attempted = False
            return False

        batch = self._raw_batch
        self._raw_batch = None

        if batch is not None:
            # A later calculation must read a fresh trajectory batch.
            self._raw_batch_attempted = False

        if batch is None and not self._raw_batch_attempted:
            batch = self._raw_reader.try_read_all_frames(
                expected_n_atoms=self.n_atoms,
                expected_n_frames=self.n_frames,
                max_bytes=self._batch_max_bytes,
                include_cells=False,
            )

        if batch is None:
            return False

        values, _cells = batch
        cells = self.cells
        legacy_box_lengths = None

        if self._legacy_rdf:
            legacy_box_lengths = np.ascontiguousarray(
                [cell.box_lengths for cell in cells],
                dtype=np.float64,
            )

        reference_indices = np.ascontiguousarray(
            self.reference_indices,
            dtype=np.int64,
        )
        target_indices = np.ascontiguousarray(
            self.target_indices,
            dtype=np.int64,
        )
        total_pairs = (
            self.n_frames * len(reference_indices) * len(target_indices)
        )
        useful_workers = max(
            1,
            total_pairs // self._batch_pairs_per_worker,
        )
        n_workers = min(
            self._batch_max_workers,
            cpu_count() or 1,
            self.n_frames,
            useful_workers,
        )

        boundaries = np.linspace(
            0,
            self.n_frames,
            n_workers + 1,
            dtype=np.intp,
        )
        arguments = [
            (
                values,
                cells,
                int(start),
                int(stop),
                reference_indices,
                target_indices,
                self.r_min,
                self.delta_r,
                self.n_bins,
                self._legacy_rdf,
                legacy_box_lengths,
            )
            for start, stop in zip(boundaries[:-1], boundaries[1:])
            if start < stop
        ]

        if n_workers == 1:
            histograms = [_raw_histogram_range(*arguments[0])]
        else:
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = [
                    executor.submit(_raw_histogram_range, *args)
                    for args in arguments
                ]
                histograms = [future.result() for future in futures]

        hist = np.zeros(self.n_bins, dtype=np.int64)
        for private_histogram in histograms:
            hist += private_histogram

        self.bins += hist

        return True

    def _finalize_run(
        self
    ) -> Tuple[Np1DNumberArray,
               Np1DNumberArray,
               Np1DNumberArray,
               Np1DNumberArray,
               Np1DNumberArray]:
        """
        Finalizes the RDF analysis after running.

        This method is called by the run method of the RDF class.
        It finalizes the RDF analysis after running by calculating the
        radial distribution function, cumulative coordination number,
        density-normalized shell population and ideal-gas pair-count
        residual. See :ref:`RDF output files <analysis-output-rdf>` for
        their exact definitions.

        Returns
        -------
        bin_middle_points : Np1DNumberArray
            The bin-center distances in Angstrom.
        normalized_bins : Np1DNumberArray
            The dimensionless radial distribution function ``g(r)``.
        integrated_bins : Np1DNumberArray
            The cumulative coordination number.
        normalized_bins2 : Np1DNumberArray
            The density-normalized shell population in Angstrom^3.
        differential_bins : Np1DNumberArray
            The pair-count residual relative to the ideal-gas shell count.
        """

        if self._legacy_rdf:
            return self._finalize_legacy_run()

        if self.no_intra_molecular:
            target_density = len(self.target_index_combinations[0])
        else:
            target_density = len(self.target_indices)
        target_density /= self._average_volume

        norm = self._norm(
            self.n_bins,
            self.delta_r,
            target_density,
            len(self.reference_indices),
            self.n_frames,
            self.r_min,
        )

        self.normalized_bins = self.bins / norm

        self.integrated_bins = self._integration(
            self.bins, len(self.reference_indices), self.n_frames
        )

        self.normalized_bins2 = self.bins / target_density
        self.normalized_bins2 /= len(self.reference_indices)
        self.normalized_bins2 /= self.n_frames

        self.differential_bins = self.bins - norm

        return (
            self.bin_middle_points,
            self.normalized_bins,
            self.integrated_bins,
            self.normalized_bins2,
            self.differential_bins
        )

    def _finalize_legacy_run(self):
        """Finalizes with the scalar operation order of legacy RDF."""
        n_reference = float(len(self.reference_indices))
        n_target = float(len(self.target_indices))
        n_frames = float(self.n_frames)
        target_density = n_target / self._average_volume
        pi = 3.1415926535897932384626433832795028841971693993751058209749445923078

        self.normalized_bins = np.empty(self.n_bins, dtype=np.float64)
        self.integrated_bins = np.empty(self.n_bins, dtype=np.float64)
        self.normalized_bins2 = np.empty(self.n_bins, dtype=np.float64)
        self.differential_bins = np.empty(self.n_bins, dtype=np.float64)

        integration = 0.0

        for index in range(self.n_bins):
            inner_radius = float(index) * self.delta_r
            outer_radius = (float(index) + 1.0) * self.delta_r
            shell = (math.pow(outer_radius, 3.0) - math.pow(inner_radius, 3.0))
            norm = 4.0 / 3.0
            norm = norm * pi
            norm = norm * target_density
            norm = norm * shell
            norm = norm * n_reference
            norm = norm * n_frames

            count = float(self.bins[index])
            integration = integration + count / n_reference / n_frames

            self.normalized_bins[index] = count / norm
            self.integrated_bins[index] = integration
            self.normalized_bins2[index] = (
                count / target_density / n_reference / n_frames
            )
            self.differential_bins[index] = count - norm

        return (
            self.bin_middle_points,
            self.normalized_bins,
            self.integrated_bins,
            self.normalized_bins2,
            self.differential_bins,
        )

    @property
    def n_frames(self) -> int:
        """int: The number of frames of the RDF analysis."""
        return len(self.cells)

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms of the RDF analysis."""
        return self.topology.n_atoms

    @classmethod
    def _add_to_bins(
        cls,
        distances: Np1DNumberArray,
        r_min: PositiveReal,
        delta_r: PositiveReal,
        n_bins: PositiveInt
    ) -> Np1DNumberArray:
        """
        Returns the bins of the RDF analysis based on the provided distances.

        Parameters
        ----------
        distances : Np1DNumberArray
            The distances to add to the bins of the RDF analysis.
        r_min : PositiveReal
            minimum (starting) radius of the RDF analysis
        delta_r : PositiveReal
            spacing between bins
        n_bins : PositiveInt
            number of bins

        Returns
        -------
        Np1DNumberArray
            The bins of the RDF analysis.
        """

        distances = np.floor_divide(distances - r_min, delta_r).astype(int)

        distances = distances[(distances < n_bins) & (distances >= 0)]

        return np.bincount(distances, minlength=n_bins)

    @classmethod
    def _setup_bin_middle_points(
        cls,
        n_bins: PositiveInt,
        r_min: PositiveReal,
        r_max: PositiveReal,
        delta_r: PositiveReal
    ) -> Np1DNumberArray:
        """
        Sets up the middle points of the bins of the RDF analysis for outputting the RDF analysis.

        Parameters
        ----------
        n_bins : PositiveInt
            number of bins
        r_min : PositiveReal
            minimum (starting) radius of the RDF analysis
        r_max : PositiveReal
            maximum radius of the RDF analysis
        delta_r : PositiveReal
            spacing between bins

        Returns
        -------
        Np1DNumberArray
            The middle points of the bins of the RDF analysis.
        """

        bin_middle_points = np.arange(r_min + delta_r / 2, r_max, delta_r)

        assert len(bin_middle_points) == n_bins

        return bin_middle_points

    @classmethod
    def _calculate_r_max(
        cls,
        n_bins: PositiveInt,
        delta_r: PositiveReal,
        r_min: PositiveReal,
        cells: Cells
    ) -> PositiveReal:
        """
        Calculates the maximum radius of the RDF analysis from the provided parameters.

        Parameters
        ----------
        n_bins : PositiveInt
            number of bins
        delta_r : PositiveReal
            spacing between bins
        r_min : PositiveReal
            minimum (starting) radius of the RDF analysis
        cells : Cells
            The cells of the trajectory to calculate the
            maximum radius of the RDF analysis from.

        Returns
        -------
        PositiveReal
            maximum radius of the RDF analysis
        """

        r_max = delta_r * n_bins + r_min
        r_max = cls._check_r_max(r_max, cells)

        return r_max

    @classmethod
    def _check_r_max(cls, r_max: PositiveReal, cells: Cells) -> PositiveReal:
        """
        Checks if the provided maximum radius is larger than the 
        maximum allowed radius according to the box vectors of the trajectory.

        Parameters
        ----------
        r_max : PositiveReal
            maximum radius of the RDF analysis
        cells : Cells
            The cells of the trajectory to check the
            maximum radius of the RDF analysis against.

        Returns
        -------
        PositiveReal
            maximum radius of the RDF analysis if it is smaller
            than the maximum allowed radius according to the box
            vectors of the trajectory, than the maximum allowed 
            radius according to the box vectors of the trajectory.

        Raises
        ------
        RDFWarning
            If the calculated r_max is larger than the maximum
            allowed radius according to the box vectors of the trajectory.
        """

        if check_trajectory_pbc(cells) and r_max > cls._infer_r_max(cells):
            cls.logger.warning(
                (
                    f"The calculated r_max {r_max} is larger "
                    "than the maximum allowed radius according "
                    "to the box vectors of the trajectory "
                    f"{cls._infer_r_max(cells)}. r_max will be "
                    "set to the maximum allowed radius."
                ),
            )

            r_max = cls._infer_r_max(cells)

        return r_max

    @classmethod
    def _calculate_n_bins(
        cls, delta_r: PositiveReal, r_max: PositiveReal, r_min: PositiveReal
    ) -> Tuple[PositiveInt, PositiveReal]:
        """
        Calculates the number of bins of the RDF analysis from the provided parameters.

        The number of bins is calculated as the number of bins that 
        fit in the range between r_min and r_max. The maximum radius
        is re-calculated from the number of bins and delta_r to ensure
        that the maximum radius is a multiple of delta_r.

        Parameters
        ----------
        delta_r : PositiveReal
            spacing between bins
        r_max : PositiveReal
            maximum radius of the RDF analysis
        r_min : PositiveReal
            minimum (starting) radius of the RDF analysis

        Returns
        -------
        PositiveInt
            number of bins of the RDF analysis
        PositiveReal
            maximum radius of the RDF analysis
        """

        n_bins = int((r_max - r_min) / delta_r)
        r_max = delta_r * n_bins + r_min

        return n_bins, r_max

    @classmethod
    def _infer_r_max(cls, cells: Cells) -> PositiveReal:
        """
        Infers the maximum radius of the RDF analysis from 
        the box vectors of the trajectory.

        If the trajectory is in vacuum, an RDFError is raised as
        the maximum radius cannot be inferred from the box vectors.

        Parameters
        ----------
        cells : Cells
            The cells of the trajectory to infer the maximum
            radius of the RDF analysis from.

        Returns
        -------
        r_max: PositiveReal
            The maximum radius of the RDF analysis.

        Raises
        ------
        RDFError
            If the trajectory is in vacuum.
        """

        if not check_trajectory_pbc(cells):
            cls.logger.error(
                (
                    "To infer r_max of the RDF analysis, "
                    "the trajectory cannot be a vacuum trajectory. "
                    "Please specify r_max manually or use the "
                    "combination n_bins and delta_r."
                ),
                exception=RDFError
            )

        return np.min([cell.box_lengths for cell in cells]) / 2.0

    @classmethod
    def _norm(
        cls,
        n_bins: int,
        delta_r: PositiveReal,
        target_density: PositiveReal,
        n_reference_indices: int,
        n_frames: int,
        r_min: PositiveReal = 0.0,
    ) -> Np1DNumberArray:
        """
        Calculates the normalization of the RDF analysis 
        based on a spherical shell model.

        Parameters
        ----------
        n_bins : int
            The number of bins of the RDF analysis.
        delta_r : PositiveReal
            The spacing between bins of the RDF analysis.
        target_density : PositiveReal
            The target density of the RDF analysis.
        n_reference_indices : int
            The number of reference indices of the RDF analysis.
        n_frames : int
            The number of frames of the RDF analysis.
        r_min : PositiveReal, optional
            The lower edge of the first shell in Angstrom, by default 0.0.

        Returns
        -------
        Np1DNumberArray
            The normalization of the RDF analysis.
        """

        volume_prefactor = 4.0 / 3.0 * np.pi

        small_radius_range = r_min + np.arange(n_bins) * delta_r
        large_radius_range = small_radius_range + delta_r
        delta_volume = large_radius_range**3 - small_radius_range**3
        volume = volume_prefactor * delta_volume

        return volume * target_density * n_reference_indices * n_frames

    @classmethod
    def _integration(
        cls, bins: Np1DNumberArray, n_reference_indices: int, n_frames: int
    ) -> Np1DNumberArray:
        """
        Calculates the integrated RDF analysis. 
        The integral is calculated using a cumulative sum.

        Parameters
        ----------
        bins : Np1DNumberArray
            The bins of the RDF analysis.
        n_reference_indices : int
            The number of reference indices of the RDF analysis.
        n_frames : int
            The number of frames of the RDF analysis.

        Returns
        -------
        Np1DNumberArray
            The integrated RDF analysis.
        """

        return np.cumsum(bins) / (n_reference_indices * n_frames)
