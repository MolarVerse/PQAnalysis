"""
A module containing the RDF class. The RDF class is used to calculate the radial distribution function of a system.
"""

from __future__ import annotations

import numpy as np
import warnings

from beartype.typing import List, Tuple
from tqdm.auto import tqdm

from . import RDFError, RDFWarning
from ..types import Np1DIntArray, Np1DNumberArray, PositiveInt, PositiveReal
from ..core import Atom, distance
from ..traj import Trajectory


class RadialDistributionFunction:
    """
    A class for calculating the radial distribution function of a system.
    """

    def __init__(self,
                 traj: Trajectory,
                 reference_species: List[Atom] | List[str] | Np1DIntArray,
                 target_species: List[Atom] | List[str] | Np1DIntArray,
                 use_full_atom_info: bool = False,
                 n_bins: PositiveInt | None = None,
                 delta_r: PositiveReal | None = None,
                 r_max: PositiveReal | None = None,
                 r_min: PositiveReal = 0.0,
                 ):
        """
        Initializes the RDF analysis.

        It sets up the bins of the RDF analysis based on the provided parameters.

        Furthermore, it sets up the reference and target indices based on the provided reference and target species.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to perform the RDF analysis on.
        reference_species : List[Atom] | List[str] | Np1DIntArray
            The reference species of the RDF analysis.
        target_species : List[Atom] | List[str] | Np1DIntArray
            The target species of the RDF analysis.
        use_full_atom_info : bool, optional
            Whether to use the full atom information of the trajectory or not, by default False. For more information, see the documentation of the indices_from_atoms method of the System class.
        n_bins : PositiveInt | None, optional
            number of bins, by default None
        delta_r : PositiveReal | None, optional
            delta r between bins, by default None
        r_max : PositiveReal | None, optional
            maximum radius of the RDF analysis, by default None
        r_min : PositiveReal, optional
            minimum (starting) radius of the RDF analysis, by default 0.0
        """

        self.traj = traj
        self.reference_species = reference_species
        self.target_species = target_species

        self.reference_indices = self.traj[0].system.indices_from_atoms(
            atoms=self.reference_species, use_full_atom_info=use_full_atom_info)
        self.target_indices = self.traj[0].system.indices_from_atoms(
            atoms=self.target_species, use_full_atom_info=use_full_atom_info)

        self.setup_bins(n_bins=n_bins, delta_r=delta_r,
                        r_max=r_max, r_min=r_min)

    def setup_bins(self,
                   n_bins: PositiveInt | None = None,
                   delta_r: PositiveReal | None = None,
                   r_max: PositiveReal | None = None,
                   r_min: PositiveReal = 0.0
                   ):
        """
        Sets up the bins of the RDF analysis.

        This method is called by the __init__ method of the RadialDistributionFunction class.
        It sets up the bins of the RDF analysis based on the provided parameters.

        if n_bins and delta_r are both specified, r_max is calculated from these parameters.

        if n_bins and r_max are both specified, delta_r is calculated from these parameters.

        if delta_r and r_max are both specified, n_bins is calculated from these parameters.

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
            If the trajectory is not fully periodic or fully in vacuum. Meaning that some frames are in vacuum and others are periodic.
        RDFError
            If n_bins and delta_r are both not specified.
        RDFError
            If n_bins, delta_r and r_max are all specified. This would lead to ambiguous results.
        """

        self.r_min = r_min

        # check if the trajectory is fully periodic or fully in vacuum
        if not self.traj.check_PBC and not self.traj.check_vacuum():
            raise RDFError(
                "The provided trajectory is not fully periodic or in vacuum, meaning that some frames are in vacuum and others are periodic. This is not supported by the RDF analysis.")

        # check if n_bins and delta_r are both not specified
        if n_bins is None and delta_r is None:
            raise RDFError(
                "Either n_bins or delta_r must be specified.")

        # check if n_bins, delta_r and r_max are all specified
        elif n_bins is not None and delta_r is not None and r_max is not None:
            raise RDFError(
                "It is not possible to specify all of n_bins, delta_r and r_max in the same RDF analysis as this would lead to ambiguous results.")

        # set r_max based on the provided parameters n_bins and delta_r
        elif n_bins is not None and delta_r is not None:
            self.n_bins = n_bins
            self.delta_r = delta_r
            self.r_max = self._calculate_r_max(n_bins, delta_r, r_min)

        else:
            self.r_max = r_max

            if r_max is None:
                self.r_max = _infer_r_max(self.traj)

            self.r_max = self._check_r_max(self.r_max)

            if n_bins is None:
                self.delta_r = delta_r
                self.n_bins, self.r_max = _calculate_n_bins(
                    delta_r, self.r_max, r_min)

            else:
                self.n_bins = n_bins
                self.delta_r = (self.r_max - self.r_min) / self.n_bins

        self._setup_bin_middle_points()
        self.bins = np.zeros(self.n_bins)

    def _calculate_r_max(self, n_bins: PositiveInt, delta_r: PositiveReal, r_min: PositiveReal) -> PositiveReal:
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

        Returns
        -------
        PositiveReal
            maximum radius of the RDF analysis
        """
        r_max = delta_r * n_bins + r_min
        r_max = self._check_r_max(r_max)

        return r_max

    def _check_r_max(self, r_max: PositiveReal) -> PositiveReal:
        """
        Checks if the provided maximum radius is larger than the maximum allowed radius according to the box vectors of the trajectory.

        Parameters
        ----------
        r_max : PositiveReal
            maximum radius of the RDF analysis

        Returns
        -------
        PositiveReal
            maximum radius of the RDF analysis if it is smaller than the maximum allowed radius
            according to the box vectors of the trajectory, otherwise the maxassert r_max == 101.0
        Raises
        ------
        RDFWarning
            If the calculated r_max is larger than the maximum allowed radius according to the box vectors of the trajectory.
        """
        if self.traj.check_PBC() and r_max > _infer_r_max(self.traj):
            warnings.warn(
                f"The calculated r_max {r_max} is larger than the maximum allowed radius \
                according to the box vectors of the trajectory {_infer_r_max(self.traj)}. \
                r_max will be set to the maximum allowed radius.", RDFWarning)

            r_max = _infer_r_max(self.traj)

        return r_max

    def _setup_bin_middle_points(self):
        """
        Sets up the middle points of the bins of the RDF analysis for outputting the RDF analysis.
        """
        self.bin_middle_points = np.arange(
            self.r_min + self.delta_r / 2, self.r_max, self.delta_r)

        assert len(self.bin_middle_points) == self.n_bins

    def run(self) -> Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]:
        """
        Runs the RDF analysis.

        Returns
        -------assert r_max == 101.0
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
        self._average_volume = np.mean(self.traj.box_volumes)
        self._reference_density = len(
            self.reference_indices) / self._average_volume
        self._target_density = len(self.target_indices) / self._average_volume

        for frame in tqdm(self.traj):
            reference_positions = frame.pos[self.reference_indices]
            target_positions = frame.pos[self.target_indices]

            for reference_position in reference_positions:
                distances = distance(reference_position,
                                     target_positions, frame.cell)

                self._add_to_bins(distances)

        normalized_bins = self.bins / self._norm()
        integrated_bins = _integration(self.bins, len(
            self.reference_indices), len(self.traj))
        normalized_bins2 = self.bins / self._target_density / \
            len(self.reference_indices) / len(self.traj)
        differential_bins = self.bins - self._norm()

        return self.bin_middle_points, normalized_bins, integrated_bins, normalized_bins2, differential_bins

    def _add_to_bins(self, distances: Np1DNumberArray):
        """
        Adds the provided distances to the bins of the RDF analysis based on the provided parameters.

        Parameters
        ----------
        distances : Np1DNumberArray
            The distances to add to the bins of the RDF analysis.
        """
        distances = np.floor_divide(
            distances - self.r_min, self.delta_r).astype(int)

        distances = distances[(distances < self.n_bins) & (distances >= 0)]

        self.bins += np.bincount(distances, minlength=self.n_bins)

    def _norm(self) -> Np1DNumberArray:
        """
        Calculates the normalization of the RDF analysis based on a spherical shell model.

        Returns
        -------
        Np1DNumberArray
            The normalization of the RDF analysis.
        """

        volume = 4.0 / 3.0 * np.pi * \
            (np.arange(1, self.n_bins + 1)**3 -
             np.arange(0, self.n_bins) ** 3) * self.delta_r ** 3

        return volume * self._target_density * len(self.reference_indices) * len(self.traj)


def _calculate_n_bins(delta_r: PositiveReal, r_max: PositiveReal, r_min: PositiveReal) -> Tuple[PositiveInt, PositiveReal]:
    """
    Calculates the number of bins of the RDF analysis from the provided parameters.

    The number of bins is calculated as the number of bins that fit in the range between r_min and r_max.
    The maximum radius is re-calculated from the number of bins and delta_r to ensure that the maximum radius is a multiple of delta_r.

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


def _infer_r_max(traj: Trajectory):
    """
    Infers the maximum radius of the RDF analysis from the box vectors of the trajectory.

    If the trajectory is in vacuum, an RDFError is raised as the maximum radius cannot be inferred from the box vectors.

    Returns
    -------
    r_max: Real
        The maximum radius of the RDF analysis.

    Raises
    ------
    RDFError
        If the trajectory is in vacuum.
    """
    if not traj.check_PBC():
        raise RDFError(
            "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r.")

    return np.min(traj.box_lengths) / 2.0


def _integration(bins: Np1DNumberArray, n_reference_indices: int, n_frames: int) -> Np1DNumberArray:
    """
    Calculates the integrated RDF analysis. The integral is calculated using a cumulative sum.

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
