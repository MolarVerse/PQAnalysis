"""
A module containing the ShakeTopologyGenerator class.

...

Classes
-------
ShakeTopologyGenerator
    A class for generating the shake topology for a given trajectory.
"""
import numpy as np

from beartype.typing import List

from ..core import Atom
from ..traj import Trajectory
from ..types import Np1DIntArray, Np2DIntArray, Np1DNumberArray
from ..io import BaseWriter


class ShakeTopologyGenerator:
    """
    A class for generating the shake topology for a given trajectory.

    Attributes
    ----------
    indices : Np1DIntArray
        The indices of the atoms to use for the topology.
    target_indices : Np1DIntArray
        The indices of the target atoms for the shaked atoms.
    distances : Np1DNumberArray
        The average distances between the shaked atoms and the target atoms.
    """

    def __init__(self,
                 atoms: List[Atom] | List[str] | Np1DIntArray | None = None,
                 use_full_atom_info: bool = False
                 ) -> None:
        """
        Initializes the ShakeTopologyGenerator with the given parameters.

        It can be initialized with a list of atoms, a list of atom names, or a numpy 1d array of atom indices or None.

            - If atoms is None shake distances for all atoms in the system will be computed.

            - If atoms is a list of atoms, shake distances for the given atoms will be computed.
              For the list of atoms, and additional parameter use_full_atom_info can be given, which
              determines if the full atom information (name, index, mass) is used for the selection or 
              only the element type.

            - If atoms is a list of atom names, shake distances for the given atom type names will be computed.

            - If atoms is a numpy 1d array of atom indices, shake distances for the given atom indices will be computed.

        Parameters
        ----------
        atoms : List[Atom] | List[str] | Np1DIntArray | None, optional
            The atoms to use for the topology, by default None
        use_full_atom_info : bool, optional
            If True, the full atom information (name, index, mass) is used for the selection, by default False
            Is always ignored if atoms is not a list of atom objects.
        """

        self._use_full_atom_info = False

        if atoms is None:
            self.atoms = None
        elif isinstance(atoms[0], Atom):
            self.atoms = atoms
            self._use_full_atom_info = use_full_atom_info
        elif isinstance(atoms[0], str):
            self.atoms = [Atom(name) for name in atoms]
        else:
            self.atoms = atoms

    def generate_topology(self, trajectory: Trajectory) -> None:
        """
        Generates a tuple of indices, target_indices, and distances for the given trajectory.

        The generated numpy arrays represent all important information about the shake topology for the given trajectory.

            - indices: The indices of the atoms to use for the topology.
            - target_indices: The indices of the target atoms for the shaked atoms.
            - distances: The average distances between the shaked atoms and the target atoms.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to generate the shake topology for.
        """

        start_frame = trajectory[0]
        target_indices, distances = start_frame.system.nearest_neighbours(
            n=1, atoms=self.atoms, use_full_atom_info=self._use_full_atom_info)

        target_indices = target_indices.flatten()
        distances = distances.flatten()

        indices = start_frame.system.indices_from_atoms(
            self.atoms, use_full_atom_info=self._use_full_atom_info)

        distances = [distances]

        for frame in trajectory[1:]:
            pos = frame.pos[indices]
            target_pos = frame.pos[target_indices]

            delta_pos = pos - target_pos

            delta_pos = frame.cell.image(delta_pos)

            distances.append(np.linalg.norm(delta_pos, axis=1))

        self.indices = indices
        self.target_indices = target_indices
        self.distances = np.mean(np.array(distances), axis=0)

    def average_equivalents(self, indices: List[Np1DIntArray] | Np2DIntArray) -> None:
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

    def write_topology(self, filename: str | None = None) -> None:
        """
        Writes the topology to a file.

        The topology is written to a file with the given filename.
        The file will contain the indices, target_indices, and distances of the topology.
        If no filename is given, the topology will be written to stdout.

        Parameters
        ----------
        filename : str
            The filename to write the topology to.
        """

        writer = BaseWriter(filename)
        writer.open()

        print(
            f"SHAKE {len(self.indices)}  {len(np.unique(self.target_indices))}  0", file=writer.file)
        for i, index in enumerate(self.indices):
            target_index = self.target_indices[i]
            distance = self.distances[i]

            print(f"{index+1} {target_index+1} {distance}", file=writer.file)

        print("END", file=writer.file)

    @property
    def atoms(self) -> List[Atom] | Np1DIntArray | None:
        """
        The atoms to use for the topology.
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms: List[Atom] | Np1DIntArray | None) -> None:
        """
        Sets the atoms to use for the topology.
        """
        self._atoms = atoms
