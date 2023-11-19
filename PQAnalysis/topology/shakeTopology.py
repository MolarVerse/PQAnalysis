import numpy as np

from beartype.typing import Tuple, List

from ..core.atom import Atom
from ..traj.trajectory import Trajectory
from ..types import Np1DIntArray, Np1DNumberArray


class ShakeTopologyGenerator:
    def __init__(self,
                 atoms: List[Atom] | List[str] | Np1DIntArray | None = None,
                 use_full_atom_info: bool = False
                 ) -> None:

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

    def generate_topology(self, trajectory: Trajectory) -> Tuple[Np1DIntArray, Np1DIntArray, Np1DNumberArray]:
        """
        Generates the topology.
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

            frame_distances = np.linalg.norm(delta_pos, axis=1)

            distances.append(frame_distances)

        return indices, target_indices, np.mean(np.array(distances), axis=1)

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
