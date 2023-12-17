"""
A module containing the RDF class. The RDF class is used to calculate the radial distribution function of a system.
"""

import numpy as np

from beartype.typing import List
from tqdm.auto import tqdm

from . import RDFError
from ..types import Np1DIntArray, Np1DNumberArray, PositiveInt, PositiveReal
from ..core import Atom, distance
from ..traj import Trajectory


class RadialDistributionFunction:
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

        self.r_min = r_min

        if not self.traj.check_PBC and not self.traj.check_vacuum():
            raise RDFError(
                "The provided trajectory is not fully periodic or in vacuum, meaning that some frames are in vacuum and others are periodic. This is not supported by the RDF analysis.")

        if n_bins is None and delta_r is None:
            raise RDFError(
                "Either n_bins or delta_r must be specified.")

        elif n_bins is not None and delta_r is not None and r_max is not None:
            raise RDFError(
                "It is not possible to specify all of n_bins, delta_r and r_max in the same RDF analysis as this would lead to ambiguous results.")

        elif n_bins is not None and delta_r is not None:
            self.n_bins = n_bins
            self.delta_r = delta_r
            self.r_max = self.delta_r * self.n_bins + self.r_min

        else:
            if r_max is None:
                self.r_max = self._infer_r_max()
            else:
                self.r_max = r_max

            if n_bins is None:
                self.delta_r = delta_r
                self.n_bins = int((self.r_max - self.r_min) / self.delta_r)
                self.r_max = self.delta_r * self.n_bins + self.r_min

            else:
                self.n_bins = n_bins
                self.delta_r = (self.r_max - self.r_min) / self.n_bins

        self._setup_bin_middle_points()
        self.bins = np.zeros(self.n_bins)

    def _infer_r_max(self):
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
        if self.traj.check_vacuum():
            raise RDFError(
                "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r.")

        return np.min(self.traj.box_lengths) / 2.0

    def _setup_bin_middle_points(self):
        """
        Sets up the middle points of the bins of the RDF analysis for outputting the RDF analysis.
        """
        self.bin_middle_points = np.arange(
            self.r_min + self.delta_r / 2, self.r_max, self.delta_r)

        assert len(self.bin_middle_points) == self.n_bins

    def run(self):
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
        integrated_bins = self._integration()
        normalized_bins2 = self.bins / self._target_density / \
            len(self.reference_indices) / len(self.traj)
        differential_bins = self.bins - self._norm()

        return self.bin_middle_points, normalized_bins, integrated_bins, normalized_bins2, differential_bins

    def _add_to_bins(self, distances: Np1DNumberArray):
        distances = np.floor_divide(
            distances - self.r_min, self.delta_r).astype(int)
        distances = distances[(distances < self.n_bins) & (distances >= 0)]

        self.bins += np.bincount(distances, minlength=self.n_bins)

    def _norm(self) -> Np1DNumberArray:
        return 4.0 / 3.0 * np.pi * self._target_density * (np.arange(1, self.n_bins + 1)**3 - np.arange(0, self.n_bins) ** 3) * self.delta_r ** 3 * len(self.reference_indices) * len(self.traj)

    def _integration(self) -> Np1DNumberArray:
        return np.cumsum(self.bins / len(self.reference_indices) / len(self.traj))
