from multimethod import multimethod

from ..core.atom import Atom
from ..traj.trajectory import Trajectory
from ..types import Np1DIntArray


class ShakeTopology:
    def __init__(self, atoms: list[Atom] | None = None) -> None:
        if atoms is None:
            atoms = []

        self.atoms = atoms

    @multimethod
    def generate_topology(self, trajectory: Trajectory, type: Atom) -> None:
        """
        Generates the topology.
        """

        for frame in trajectory:
            pass

        raise NotImplementedError

    @multimethod
    def generate_topology(self, trajectory: Trajectory, indices: Np1DIntArray) -> None:
        """
        Generates the topology.
        """
        raise NotImplementedError

    ##############
    # Properties #
    ##############

    @property
    def atoms(self) -> list[Atom]:
        """
        The atoms in the system.

        Returns
        -------
        list[Atom]
            The atoms in the system.
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms: list[Atom]) -> None:
        """
        The atoms in the system.

        Parameters
        ----------
        atoms : list[Atom]
            The atoms in the system.
        """
        self._atoms = atoms
