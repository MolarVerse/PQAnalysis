"""
A module containing the ShakeTopologyGenerator class.
"""

import logging
import numpy as np

from beartype.typing import List

from PQAnalysis import __package_name__
from PQAnalysis.traj import Trajectory
from PQAnalysis.types import Np1DIntArray, Np2DIntArray
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.io import (
    FileWritingMode,
    TopologyFileWriter,
)

from .selection import SelectionCompatible, Selection
from .bonded_topology import BondedTopology, Bond



class ShakeTopologyGenerator:

    """
    A class for generating the shake topology for a given trajectory
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        selection: SelectionCompatible = None,
        use_full_atom_info: bool = False
    ) -> None:
        """
        Parameters
        ----------
        selection : SelectionCompatible, optional
            Selection is either a selection object or any object that can be
            initialized via 'Selection(selection)'. default None (all atoms)
        use_full_atom_info : bool, optional
            If True, the full atom information (name, index, mass) is used
            for the selection, by default False
            Is always ignored if atoms is not a list of atom objects.
        """

        self._use_full_atom_info = use_full_atom_info
        self.selection = Selection(selection)

        self.indices = None
        self.target_indices = None
        self.distances = None
        self._topology = None
        self.line_comments = None

    @runtime_type_checking
    def generate_topology(self, trajectory: Trajectory) -> None:
        """
        Generates a tuple of indices, target_indices, and distances for the given trajectory.

        The generated numpy arrays represent all important information about the shake
        topology for the given trajectory.

            - indices: The indices of the atoms to use for the topology.
            - target_indices: The indices of the target atoms for the shaked atoms.
            - distances: The average distances between the shaked atoms and the target atoms.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to generate the shake topology for.
        """

        atomic_system = trajectory[0]
        self._topology = trajectory.topology

        indices = self.selection.select(
            self._topology, self._use_full_atom_info
        )

        target_indices, distances = atomic_system.nearest_neighbours(
            n=1, selection=indices, use_full_atom_info=self._use_full_atom_info)

        target_indices = target_indices.flatten()
        distances = [distances.flatten()]

        for frame in trajectory[1:]:
            pos = frame.pos[indices]
            target_pos = frame.pos[target_indices]

            delta_pos = pos - target_pos

            delta_pos = frame.cell.image(delta_pos)

            distances.append(np.linalg.norm(delta_pos, axis=1))

        self.indices = indices
        self.target_indices = target_indices
        self.distances = np.mean(np.array(distances), axis=0)

    @runtime_type_checking
    def average_equivalents(
        self,
        indices: List[Np1DIntArray] | Np2DIntArray,
        comments: List[str] | None = None
    ) -> None:
        """
        Averages the distances for equivalent atoms.

        The parameter indices is a numpy 2d array of equivalent atom indices.
        Each row of the array contains the indices of equivalent atoms.
        All of the equivalent atoms will be averaged to a single distance.

        Parameters
        ----------
        indices : List[Np1DIntArray] | Np2DIntArray
            The indices of the equivalent atoms.
        comments : List[str], optional
            The line comments for the averaged distances, by default None
        """

        for equivalent_indices in indices:
            _indices = np.nonzero(np.isin(self.indices, equivalent_indices))[0]

            mean_distance = np.mean(self.distances[_indices])

            self.distances[_indices] = mean_distance

        if comments is not None:
            if len(comments) != len(indices):
                self.logger.error(
                    "The number of comments does not match the number of indices.",
                    exception=PQValueError
                )

            self.line_comments = [""] * len(self.indices)

            for i, equivalent_indices in enumerate(indices):
                _indices = np.nonzero(
                    np.isin(self.indices, equivalent_indices)
                )[0]

                for index in _indices:
                    self.line_comments[index] = comments[i]

    @runtime_type_checking
    def add_comments(self, comments: List[str]) -> None:
        """
        Adds comments to the topology.

        Parameters
        ----------
        comments : List[str]
            The comments to add to each line of the shake topology.
        """

        if self.indices is None or len(comments) != len(self.indices):
            self.logger.error(
                "The number of comments does not match the number of indices.",
                exception=PQValueError
            )

        if self.line_comments is None:
            self.line_comments = [""] * len(self.indices)

        for i, comment in enumerate(comments):
            self.line_comments[i] = comment

    @runtime_type_checking
    def write_topology(
        self,
        filename: str | None = None,
        mode: FileWritingMode | str = "w"
    ) -> None:
        """
        Writes the topology to a file.

        The topology is written to a file with the given filename.
        The file will contain the indices, target_indices, and distances of the topology.
        If no filename is given, the topology will be written to stdout.

        Parameters
        ----------
        filename : str
            The filename to write the topology to.
        mode : FileWritingMode | str, optional
            The writing mode, by default "w". The following modes are available:
            - "w": write
            - "a": append
            - "o": overwrite
        """

        shake_bonds = []
        if self.line_comments is None:
            comments = [None] * len(self.indices)
        else:
            comments = self.line_comments

        for i, index in enumerate(self.indices):
            index = index + 1
            target_index = self.target_indices[i] + 1
            distance = self.distances[i]
            comment = comments[i]

            bond = Bond(
                index1=index,
                index2=target_index,
                equilibrium_distance=distance,
                is_shake=True,
                comment=comment
            )

            shake_bonds.append(bond)

        topology = BondedTopology(shake_bonds=shake_bonds)

        writer = TopologyFileWriter(filename, mode)
        writer.write(topology)

    @property
    def selection_object(self) -> SelectionCompatible:
        """SelectionCompatible: The selection object."""
        return self.selection.selection_object
