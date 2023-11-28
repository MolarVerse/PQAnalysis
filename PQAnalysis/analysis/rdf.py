import numpy as np

from beartype.typing import List

from . import RDFError
from ..types import Np1DIntArray, PositiveInt, PositiveReal
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

        ####### LOGGER FOR BIN INFO #######

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
                self._infer_r_max()
            else:
                self.r_max = r_max

            if n_bins is None:
                self.delta_r = delta_r
                self.n_bins = (self.r_max - self.r_min) // self.delta_r

            else:
                self.n_bins = n_bins
                self.delta_r = (self.r_max - self.r_min) / self.n_bins

        self._setup_bin_middle_points()
        self.bins = np.zeros(self.n_bins)

    def _infer_r_max(self):
        if self.traj.check_vacuum():
            raise RDFError(
                "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r.")

        self.r_max = np.min(self.traj.box_lengths)

    def _setup_bin_middle_points(self):

        self.bin_middle_points = np.arange(
            self.r_min + self.delta_r / 2, self.r_max, self.delta_r)

        assert len(self.bin_middle_points) == self.n_bins

    def run(self):
        for frame in self.traj:
            reference_positions = frame.pos[self.reference_indices]
            target_positions = frame.pos[self.target_indices]

            for reference_position in reference_positions:
                distances = distance(reference_position,
                                     target_positions, frame.cell)

                distances -= self.r_min
                distances //= self.delta_r

                for distance in distances:
                    if distance < self.n_bins and distance >= 0:
                        self.bins[distance] += 1

    def normalise(self):
        self.bins /= self.traj.n_frames * \
            self.reference_indices.shape[0] * self.delta_r * \
            self.target_indices.shape[0] * 4 * np.pi * \
            self.bin_middle_points**2
