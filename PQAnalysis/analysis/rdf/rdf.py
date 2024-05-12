"""
A module containing the RDF class. The RDF class is used
to calculate the radial distribution of a reference
selection to a target selection. The radial distribution
function (RDF) is a measure of the probability density 
of finding a particle at a distance r from another particle. 
"""

import logging

# 3rd party imports
import numpy as np

# 3rd party imports
from beartype.typing import Tuple
from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis.config import with_progress_bar
from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import distance, Cells
from PQAnalysis.traj import Trajectory, check_trajectory_pbc, check_trajectory_vacuum
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import TrajectoryReader
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import RDFError



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
    """

    _use_full_atom_default = False
    _no_intra_molecular_default = False
    _r_min_default = 0.0
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        traj: Trajectory | TrajectoryReader,
        reference_species: SelectionCompatible,
        target_species: SelectionCompatible,
        use_full_atom_info: bool = False,
        no_intra_molecular: bool = False,
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
        use_full_atom_info : bool, optional
            Whether to use the full atom information of the trajectory
            or not, by default None (False).
        no_intra_molecular : bool, optional
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

        self.cells = traj.cells

        if isinstance(traj, TrajectoryReader):
            # lazy loading of trajectory from file(s)
            self.frame_generator = traj.frame_generator()
        elif len(traj) > 0:
            # use trajectory object as iterator
            self.frame_generator = iter(traj)
        else:
            self.logger.error(
                "Trajectory cannot be of length 0.",
                exception=RDFError
            )

        self.first_frame = next(self.frame_generator)
        self.topology = traj.topology

        self._setup_bins(
            n_bins=n_bins,
            delta_r=delta_r,
            r_max=r_max,
            r_min=self.r_min
        )

        self.reference_indices = self.reference_selection.select(
            self.topology,
            self.use_full_atom_info
        )

        self.target_indices = self.target_selection.select(
            self.topology,
            self.use_full_atom_info
        )

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
                n_bins,
                delta_r,
                r_min,
                self.cells
            )

            self.n_bins, self.r_max = self._calculate_n_bins(
                delta_r,
                self.r_max,
                r_min
            )

        else:

            self.r_max = r_max if r_max is not None else self._infer_r_max(
                self.cells
            )

            self.r_max = self._check_r_max(self.r_max, self.cells)

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
            self.n_bins,
            self.r_min,
            self.r_max,
            self.delta_r
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
                self.cells) and not check_trajectory_vacuum(self.cells):
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

        This method runs the RDF analysis and returns the 
        middle points of the bins of the RDF analysis, the
        normalized bins of the RDF analysis based on the 
        spherical shell model, the integrated bins of the
        RDF analysis, the normalized bins of the RDF 
        analysis based on the number of atoms in the 
        system and the differential bins of the RDF 
        analysis based on the spherical shell model.

        This method will display a progress bar by default.
        This can be disabled by setting with_progress_bar to False.

        Returns
        -------
        bin_middle_points : Np1DNumberArray
            The middle points of the bins of the RDF analysis.
        normalized_bins : Np1DNumberArray
            The normalized bins of the RDF analysis based 
            on the spherical shell model.
        integrated_bins : Np1DNumberArray
            The integrated bins of the RDF analysis.
        normalized_bins2 : Np1DNumberArray
            The normalized bins of the RDF analysis based
            on the number of atoms in the system.
        differential_bins : Np1DNumberArray
            The differential bins of the RDF analysis based
            on the spherical shell model.
        """

        self._initialize_run()
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

        self._average_volume = self.average_volume

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
                residue_indices = self.topology.residue_atom_indices[
                    reference_index]
                self.target_index_combinations.append(
                    np.setdiff1d(self.target_indices,
                    residue_indices)
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

        for frame in tqdm(self.frame_generator,
            total=self.n_frames,
            disable=not with_progress_bar):
            for i, reference_index in enumerate(self.reference_indices):

                if self.no_intra_molecular:
                    target_indices = self.target_index_combinations[i]
                else:
                    target_indices = self.target_indices

                reference_position = frame.pos[reference_index]
                target_positions = frame.pos[target_indices]

                distances = distance(
                    reference_position,
                    target_positions,
                    frame.cell
                )

                self.bins += self._add_to_bins(
                    distances,
                    self.r_min,
                    self.delta_r,
                    self.n_bins
                )

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
        It finalizes the RDF analysis after running by calculating
        the normalized bins of the RDF analysis based on the 
        spherical shell model, the integrated bins of the RDF 
        analysis, the normalized bins of the RDF analysis based
        on the number of atoms in the system and the differential
        bins of the RDF analysis based on the spherical shell model.

        Returns
        -------
        bin_middle_points : Np1DNumberArray
            The middle points of the bins of the RDF analysis.
        normalized_bins : Np1DNumberArray
            The normalized bins of the RDF analysis based on the spherical shell model.
        integrated_bins : Np1DNumberArray
            The integrated bins of the RDF analysis.
        normalized_bins2 : Np1DNumberArray
            The normalized bins of the RDF analysis based on the number of atoms in the system.
        differential_bins : Np1DNumberArray
            The differential bins of the RDF analysis based on the spherical shell model.
        """

        target_density = len(self.target_index_combinations[0])
        target_density /= self._average_volume

        norm = self._norm(
            self.n_bins,
            self.delta_r,
            target_density,
            len(self.reference_indices),
            self.n_frames
        )

        self.normalized_bins = self.bins / norm

        self.integrated_bins = self._integration(
            self.bins,
            len(self.reference_indices),
            self.n_frames
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
        cls,
        delta_r: PositiveReal,
        r_max: PositiveReal,
        r_min: PositiveReal
    ) -> Tuple[PositiveInt,
        PositiveReal]:
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
        n_frames: int
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

        Returns
        -------
        Np1DNumberArray
            The normalization of the RDF analysis.
        """

        surface_prefactor = 4.0 / 3.0 * np.pi

        small_radius_range = np.arange(0, n_bins)
        large_radius_range = np.arange(1, n_bins + 1)
        delta_radius_range = large_radius_range**3 - small_radius_range**3

        delta_volume = delta_r**3 * delta_radius_range
        volume = surface_prefactor * delta_volume

        return volume * target_density * n_reference_indices * n_frames

    @classmethod
    def _integration(
        cls,
        bins: Np1DNumberArray,
        n_reference_indices: int,
        n_frames: int
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
