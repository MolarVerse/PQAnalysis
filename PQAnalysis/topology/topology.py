"""
A module containing the Topology class.

...

Classes
-------
Topology
    A class for representing a topology.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import Any
from numbers import Integral

from . import Residues, Residue, QMResidue, ResidueError
from ..core import Atoms, Element
from ..types import Np1DIntArray


class Topology:
    """
    A class for representing a topology.

    A topology is a collection of atoms, residues and residue ids.
    It is used to represent the topology of a system.

    Examples
    --------
    >>> from PQAnalysis.topology import Topology
    >>> import numpy as np
    >>> from PQAnalysis.core import Atom
    >>> from PQAnalysis.topology import Residue
    >>> atoms = [Atom("H", 1), Atom("O", 8), Atom("H", 1)]
    >>> residue_ids = np.array([1, 2, 2])
    >>> reference_residues = [Residue("HO", 1, 0, ["H", "O"], [1, 2], [0, 0]), Residue("H2", 2, 0, ["H", "H"], [1, 2], [0, 0])]
    >>> topology = Topology(atoms=atoms, residue_ids=residue_ids, residues=residues)
    """

    def __init__(self,
                 atoms: Atoms | None = None,
                 residue_ids: Np1DIntArray | None = None,
                 reference_residues: Residues | None = None,
                 check_residues: bool = True
                 ) -> None:
        """
        Initializes a Topology object.

        All of the parameters are optional, if they are not given, they are initialized with empty values.
        It checks if the residue ids are compatible and contiguous regarding the given residues if the list is not empty.

        It is also possible to initialize a Topology object with residues that are not referenced within the residue_ids.
        This means that it is possible to have a Topology object with residues that are not used in the system.

        Parameters
        ----------
        atoms : Atoms | None, optional
            a list of atoms, by default None (empty list)
        residue_ids : Np1DIntArray | None, optional
            a list of residue ids, by default None (empty list)
        reference_residues : Residues | None, optional
            a list of residues, by default None (empty list)
        check_residues : bool, optional
            If the residues should be checked, by default True

        Raises
        ------
        ValueError
            If the number of atoms does not match the number of residue ids.
        """

        self._check_residues = check_residues

        if atoms is None:
            self._atoms = []
            self._atomtype_names = []
        else:
            self._atoms = atoms
            self._atomtype_names = [atom.name for atom in atoms]

        if reference_residues is None:
            self._reference_residues = []
        else:
            self._reference_residues = reference_residues

        if residue_ids is None:
            residue_ids = np.zeros(len(self.atoms), dtype=int)
        if len(self.atoms) != len(residue_ids):
            raise ValueError(
                "The number of atoms does not match the number of residue ids.")

        self._residue_ids = residue_ids
        self._residues = self._setup_residues(self.residue_ids, self.atoms)

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
        """
        Returns a new Topology with the given indices.

        Parameters
        ----------
        indices : Np1DIntArray
            The indices of the atoms to return.

        Returns
        -------
        Topology
            The new Topology with the given indices.
        """
        reference_residues = self.reference_residues

        if len(reference_residues) == 0:
            reference_residues = None

        if self.n_atoms == 0:
            return Topology(reference_residues=reference_residues)

        atoms = [self.atoms[index] for index in indices]
        residue_ids = self.residue_ids[indices]

        return Topology(atoms=atoms, reference_residues=self.reference_residues, residue_ids=residue_ids, check_residues=self.check_residues)

    def _find_residue_by_id(self, id: Integral) -> Residue:
        """
        Finds a residue by its id.

        Parameters
        ----------
        id : Integral
            The id of the residue to find.

        Returns
        -------
        Residue
            The residue with the given id.

        Raises
        ------
        ResidueError
            If the residue id is not unique.
        ResidueError
            If the residue id is not found.
        """
        bool_array = np.array(
            [residue.id == id for residue in self.reference_residues])

        residue = np.argwhere(bool_array)

        if len(residue) > 1:
            raise ResidueError(f"The residue id {id} is not unique.")

        if len(residue) == 0:
            raise ResidueError(f"The residue id {id} was not found.")

        return residue[0]

    def _setup_residues(self, residue_ids: Np1DIntArray, atoms: Atoms) -> Residues:
        """
        Sets up the residues of the topology.

        It checks if the residue ids are contiguous and compatible with the reference residues.
        If check_residues is False, it does not check the residues and just returns an empty list.

        Parameters
        ----------
        residue_ids : Np1DIntArray
            The residue ids of the topology.
        atoms : Atoms
            The atoms of the topology.

        Returns
        -------
        Residues
            The residues of the topology.

        Raises
        ------
        ResidueError
            If the residue ids are not contiguous.
        ResidueError
            If the reference residues are not empty and residue_ids with 0 don't have any element information.
            This problem can be avoided by setting 'check_residues' to False.
        """
        residues = []

        if len(self.reference_residues) == 0 or not self.check_residues:
            return residues

        atom_counter = 0
        while atom_counter < len(residue_ids):
            if residue_ids[atom_counter] == 0:
                if atoms[atom_counter].element == Element():
                    message = f"""
The element of atom {atom_counter} is not set. If any reference residues are given
the program tries to automatically deduce the residues from the residue ids and the reference residues.
This means that any atom with an unknown element raises an error. To avoid deducing residue information
please set 'check_residues' to False"""

                    raise ResidueError(message)
                else:
                    residues.append(QMResidue(atoms[atom_counter].element))
                continue

            residue = self._find_residue_by_id(residue_ids[atom_counter])

            for i in range(residue.n_atoms-1) + atom_counter:
                if residue_ids[i] != residue_ids[atom_counter]:
                    raise ResidueError(
                        f"The residue ids are not contiguous. Problems with residue {residue.name} with indices {atom_counter}-{atom_counter + residue.n_atoms-1}")

            residues.append(residue)

            atom_counter += residue.n_atoms

        return residues

    def __str__(self) -> str:
        """
        Returns a string representation of the Topology.

        Returns
        -------
        str
            The string representation of the Topology.
        """

        return f"Topology with {self.n_atoms} atoms and {self.n_residues} residues ({self.n_QM_residues} QM residues) and {self.n_unique_residues} unique residues."

    def __repr__(self) -> str:
        """
        Returns a string representation of the Topology.

        Returns
        -------
        str
            The string representation of the Topology.
        """
        return self.__str__()

    @property
    def check_residues(self) -> bool:
        """
        Returns whether the residues should be checked.

        Returns
        -------
        bool
            Whether the residues should be checked.
        """
        return self._check_residues

    @check_residues.setter
    def check_residues(self, value: bool) -> None:
        """
        Sets whether the residues should be checked.

        Parameters
        ----------
        value : bool
            Whether the residues should be checked.
        """
        self._check_residues = value
        self._setup_residues(self.residue_ids, self.atoms)

    @property
    def reference_residue_ids(self) -> Np1DIntArray:
        """
        Returns all residue ids of the reference residues.

        Returns
        -------
        Np1DIntArray
            The residue ids of the reference residues.
        """
        return np.array([residue.id for residue in self.reference_residues])

    @property
    def reference_residues(self) -> Residues:
        """
        Returns the residues of the topology.

        Returns
        -------
        Residues
            The residues of the topology.
        """
        return self._reference_residues

    @reference_residues.setter
    def reference_residues(self, value: Residues):
        """
        Sets the residues of the topology.

        Parameters
        ----------
        value : Residues
            The residues of the topology.
        """
        self._reference_residues = value

    @property
    def atoms(self) -> Atoms:
        """
        Returns the atoms of the topology.

        Returns
        -------
        Atoms
            The atoms of the topology.
        """
        return self._atoms

    @atoms.setter
    def atoms(self, value: Atoms):
        """
        Sets the atoms of the topology.

        Parameters
        ----------
        value : Atoms
            The atoms of the topology.
        """
        self._atoms = value

    @property
    def atomtype_names(self) -> list[str]:
        """
        Returns the names of the atomtypes of the topology.

        Returns
        -------
        list[str]
            The names of the atomtypes of the topology.
        """
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
        """
        Returns the residue ids of the topology.

        Returns
        -------
        Np1DIntArray
            The residue ids of the topology.
        """
        return self._residue_ids

    @property
    def residues(self) -> Residues:
        """
        Returns the residues of the topology.

        Returns
        -------
        Residues
            The residues of the topology.
        """
        return self._residues

    @property
    def n_residues(self) -> int:
        """
        Returns the number of residues in the topology.

        Returns
        -------
        int
            The number of residues in the topology.
        """
        return len(self.residues)

    def n_QM_residues(self) -> int:
        """
        Returns the number of QM residues in the topology.

        Returns
        -------
        int
            The number of QM residues in the topology.
        """
        return len([residue for residue in self.residues if isinstance(residue, QMResidue)])

    def n_MM_residues(self) -> int:
        """
        Returns the number of MM residues in the topology.

        Returns
        -------
        int
            The number of MM residues in the topology.
        """
        return self.n_residues - self.n_QM_residues

    def n_unique_residues(self) -> int:
        """
        Returns the number of unique residues in the topology.

        Returns
        -------
        int
            The number of unique residues in the topology.
        """
        return len(set(self.residue_ids))
