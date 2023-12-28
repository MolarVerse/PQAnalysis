from __future__ import annotations

import numpy as np

from beartype.typing import Any
from numbers import Integral

from ..core import Atoms
from ..types import Np1DIntArray
from . import Moltypes, MolType


class Topology:

    def __init__(self,
                 atoms: Atoms | None = None,
                 moltype_ids: Np1DIntArray | None = None,
                 moltypes: Moltypes | None = None,
                 ) -> None:

        if atoms is None:
            self._atoms = []
            self._atomtype_names = []
        else:
            self._atoms = atoms
            self._atomtype_names = [atom.name for atom in atoms]

        if moltypes is None:
            self._moltypes = []
        else:
            self._moltypes = moltypes

        if moltype_ids is None:
            moltype_ids = np.zeros(len(self.atoms), dtype=int)
        if len(self.atoms) != len(moltype_ids):
            raise ValueError(
                "The number of atoms does not match the number of mol_type_ids.")

        self._moltype_ids = moltype_ids

        # TODO: _check_moltype_ids

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

        if not isinstance(other, Topology):
            return False

        if self.n_atoms != other.n_atoms:
            return False

        if self.n_atoms == 0:
            return is_equal

        is_equal &= self.atoms == other.atoms
        is_equal &= self.moltypes == other.moltypes
        is_equal &= np.all(self.moltype_ids == other.moltype_ids)

        return bool(is_equal)

    def __getitem__(self, indices: Np1DIntArray) -> Topology:

        moltypes = self.moltypes

        if len(moltypes) == 0:
            moltypes = None

        if self.n_atoms == 0:
            return Topology(moltypes=moltypes)

        atoms = [self.atoms[index] for index in indices]

        if len(self.moltype_ids) == 0:
            return Topology(atoms=atoms, moltypes=moltypes)

        atom_counter = 0
        moltype_ids = self.moltype_ids[indices]
        while atom_counter < len(indices):
            moltype_id = moltype_ids[atom_counter]

            if moltype_id == 0:
                atom_counter += 1
                continue

            moltype = self._find_moltype_by_id(moltype_id)
            for i in range(moltype.n_atoms-1) + atom_counter:
                if moltype_ids[i] != moltype_id:
                    raise ValueError(
                        f"The moltype ids are not contiguous. Problems with residue {moltype.name} with indices {indices[atom_counter]}-{indices[atom_counter + moltype.n_atoms-1]}")

            atom_counter += moltype.n_atoms

        return Topology(atoms=atoms, moltypes=self.moltypes, moltype_ids=moltype_ids)

    def _find_moltype_by_id(self, id: Integral) -> MolType:
        bool_array = np.array([moltype.id == id for moltype in self.moltypes])

        moltype = np.argwhere(bool_array)

        if len(moltype) > 1:
            raise ValueError(f"The moltype id {moltype.id} is not unique.")

        if len(moltype) == 0:
            raise ValueError(f"The moltype id {moltype.id} was not found.")

        return moltype[0]

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
    def atomtype_names(self) -> list[str]:
        return self._atomtype_names

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

    @property
    def moltype_ids(self) -> Np1DIntArray:
        return self._moltype_ids
