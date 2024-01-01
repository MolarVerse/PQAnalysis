from __future__ import annotations

import numpy as np

from beartype.typing import Any
from numbers import Integral

from ..core import Atoms
from ..types import Np1DIntArray
from . import Residues, Residue


class Topology:

    def __init__(self,
                 atoms: Atoms | None = None,
                 residue_ids: Np1DIntArray | None = None,
                 residues: Residues | None = None,
                 ) -> None:

        if atoms is None:
            self._atoms = []
            self._atomtype_names = []
        else:
            self._atoms = atoms
            self._atomtype_names = [atom.name for atom in atoms]

        if residues is None:
            self._residues = []
        else:
            self._residues = residues

        if residue_ids is None:
            residue_ids = np.zeros(len(self.atoms), dtype=int)
        if len(self.atoms) != len(residue_ids):
            raise ValueError(
                "The number of atoms does not match the number of residue ids.")

        self._residue_ids = residue_ids

        # TODO: _check_residue_ids

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
        is_equal &= np.all(self.residue_ids == other.residue_ids)

        return bool(is_equal)

    def __getitem__(self, indices: Np1DIntArray) -> Topology:

        residues = self.residues

        if len(residues) == 0:
            residues = None

        if self.n_atoms == 0:
            return Topology(residues=residues)

        atoms = [self.atoms[index] for index in indices]

        if len(self.residue_ids) == 0:
            return Topology(atoms=atoms, residues=residues)

        atom_counter = 0
        residue_ids = self.residue_ids[indices]
        while atom_counter < len(indices):
            residue_id = residue_ids[atom_counter]

            if residue_id == 0:
                atom_counter += 1
                continue

            residue = self._find_residue_by_id(residue_id)
            for i in range(residue.n_atoms-1) + atom_counter:
                if residue_ids[i] != residue_id:
                    raise ValueError(
                        f"The residue ids are not contiguous. Problems with residue {residue.name} with indices {indices[atom_counter]}-{indices[atom_counter + residue.n_atoms-1]}")

            atom_counter += residue.n_atoms

        return Topology(atoms=atoms, residues=self.residues, residue_ids=residue_ids)

    def _find_residue_by_id(self, id: Integral) -> Residue:
        bool_array = np.array([residue.id == id for residue in self.residues])

        residue = np.argwhere(bool_array)

        if len(residue) > 1:
            raise ValueError(f"The residue id {residue.id} is not unique.")

        if len(residue) == 0:
            raise ValueError(f"The residue id {residue.id} was not found.")

        return residue[0]

    @property
    def residues(self) -> Residues:
        return self._residues

    @residues.setter
    def residues(self, value: Residues):
        self._residues = value

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
    def residue_ids(self) -> Np1DIntArray:
        return self._residue_ids
