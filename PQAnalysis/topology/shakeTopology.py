import numpy as np

from multimethod import multimethod
from beartype.typing import Tuple

from ..core.atom import Atom
from ..traj.trajectory import Trajectory
from ..types import Np1DIntArray, Np1DNumberArray


class ShakeTopologyGenerator:
    @multimethod
    def __init__(self) -> None:
        self._atoms = None
        self._use_full_atom_info = False

    @multimethod
    def __init__(self, atoms: list[Atom], use_full_atom_info: bool = False) -> None:
        self._atoms = atoms
        self._use_full_atom_info = use_full_atom_info

    @multimethod
    def __init__(self, atoms: list[str]):
        self._atoms = [Atom(name) for name in atoms]
        self._use_full_atom_info = False

    @multimethod
    def __init__(self, indices: Np1DIntArray):
        self._atoms = indices
        self._use_full_atom_info = False

    @multimethod
    def generate_topology(self, trajectory: Trajectory) -> Tuple[Np1DIntArray, Np1DIntArray, Np1DNumberArray]:
        """
        Generates the topology.
        """

        start_frame = trajectory[0]
        target_indices, distances = start_frame.nearest_neighbours(
            n=1, atoms=self.atoms, use_full_atom_info=self._use_full_atom_info)

        target_indices = target_indices.flatten()
        distances = distances.flatten()

        indices = start_frame.indices_from_atoms(
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
