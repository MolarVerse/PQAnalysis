"""
A module containing the ShakeTopologyGenerator class.
"""

import numpy as np

from beartype.typing import List

from PQAnalysis.traj import Trajectory
from PQAnalysis.types import Np1DIntArray, Np2DIntArray
from PQAnalysis.io import BaseWriter, FileWritingMode
from PQAnalysis.type_checking import runtime_type_checking

from .selection import SelectionCompatible, Selection



class ShakeTopologyGenerator:

    """
    A class for generating the shake topology for a given trajectory
    """

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
            self._topology,
            self._use_full_atom_info
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
        indices: List[Np1DIntArray] | Np2DIntArray
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
        """

        for equivalent_indices in indices:
            _indices = np.nonzero(np.in1d(self.indices, equivalent_indices))[0]

            mean_distance = np.mean(self.distances[_indices])

            self.distances[_indices] = mean_distance

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

        writer = BaseWriter(filename, mode=mode)
        writer.open()

        print(
            (
            f"SHAKE {len(self.indices)}  "
            f"{len(np.unique(self.target_indices))}  0"
            ),
            file=writer.file
        )

        for i, index in enumerate(self.indices):
            target_index = self.target_indices[i]
            distance = self.distances[i]

            print(f"{index+1} {target_index+1} {distance}", file=writer.file)

        print("END", file=writer.file)

    @property
    def selection_object(self) -> SelectionCompatible:
        """SelectionCompatible: The selection object."""
        return self.selection.selection_object
