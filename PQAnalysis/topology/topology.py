from beartype.typing import Any

from ..core import Atoms
from ..topology import Moltypes


class Topology:

    def __init__(self,
                 atoms: Atoms | None = None,
                 moltypes: Moltypes | None = None,
                 ) -> None:

        if atoms is None:
            self._atoms = []
            self._atomtypes = []
        else:
            self._atoms = atoms
            self._atomtypes = [atom.name for atom in atoms]

        if moltypes is None:
            self._moltypes = []
        else:
            self._moltypes = moltypes

            if sum([len(moltype.elements) for moltype in self.moltypes]) != len(self.atoms):
                raise ValueError(
                    "The number of atoms does not match the number of mol_types.")

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the Topology is equal to another Topology.

        Parameters
        ----------
        other : Topology
            The other Topology to compare to.

        Returns
        -------
        bool
            Whether the Topology is equal to the other Topology.
        """

        is_equal = True

        if self.n_atoms == 0:
            return is_equal
        if self.n_atoms != other.n_atoms:
            return False

        is_equal &= isinstance(other, Topology)
        is_equal &= self.atoms == other.atoms
        is_equal &= self.moltypes == other.moltypes

        return is_equal

    @property
    def moltypes(self) -> Moltypes:
        return self._moltypes

    @moltypes.setter
    def moltypes(self, value: Moltypes):
        self._moltypes = value

    @property
    def atoms(self) -> Atoms:
        return self._atoms

    @atoms.setter
    def atoms(self, value: Atoms):
        self._atoms = value

    @property
    def atomtypes(self) -> list[str]:
        return self._atomtypes

    @property
    def n_atoms(self) -> int:
        """
        Returns the number of atoms in the topology.

        Returns
        -------
        int
            The number of atoms in the topology.
        """
        return len(self.atoms)
