"""
A module containing the Topology class and related functions.
"""

# standard library
from numbers import Integral
import logging

# third-party packages
from beartype.typing import Any, Tuple, List
import numpy as np

from PQAnalysis.core.exceptions import ResidueError
from PQAnalysis.core import Residues, Residue, QMResidue, Atoms, Element
from PQAnalysis.types import Np1DIntArray
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import (
    runtime_type_checking,
    runtime_type_checking_setter,
)

from .exceptions import TopologyError
from .bonded_topology.bonded_topology import BondedTopology



class Topology:

    """
    A class for representing a topology.

    A topology is a collection of atoms, residues and residue ids.
    It is used to represent the topology of a system.

    Examples
    --------
    >>> import numpy as np
    >>> import warnings
    >>> from PQAnalysis.topology import Topology
    >>> from PQAnalysis.core import Atom, Residue
    >>> atoms = [Atom("H", 1), Atom("O", 8), Atom("H", 1), Atom("H", 1)]
    >>> residue_ids = np.array([1, 1, 2, 2])
    >>> residue1 = Residue("HO", 1, 0, ["H", "O"], np.array([1, 2]), np.array([0, 0]))
    >>> residue2 = Residue("H2", 2, 0, ["H", "H"], np.array([1, 2]), np.array([0, 0]))
    >>> residues = [residue1, residue2]
    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter("ignore")
    ...     topology = Topology(atoms=atoms, residue_ids=residue_ids, reference_residues=residues)

    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        atoms: Atoms | None = None,
        residue_ids: Np1DIntArray | None = None,
        reference_residues: Residues | None = None,
        check_residues: bool = True,
        bonded_topology: BondedTopology | None = None,
    ) -> None:
        """
        All of the parameters are optional, if they are not given,
        they are initialized with empty values. It checks if the
        residue ids are compatible and contiguous regarding the given
        residues if the list is not empty.

        It is also possible to initialize a Topology object with
        residues that are not referenced within the residue_ids. This
        means that it is possible to have a Topology object with
        residues that are not used in the system.

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
        bonded_topology : BondedTopology | None, optional
            Contains the bonded topology information, by default None

        Raises
        ------
        TopologyError
            If the number of atoms does not match the number of residue ids.

        Warns
        -----
        UserWarning
            If the bonded topology is not None. There is no check yet if the
            bonded topology is compatible with the topology. Please make sure
            that the bonded topology is compatible with the topology!
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

        if residue_ids is None or len(residue_ids) == 0:
            residue_ids = np.zeros(len(self.atoms), dtype=int)
        if len(self.atoms) != len(residue_ids):
            self.logger.error(
                "The number of atoms does not match the number of residue ids.",
                exception=TopologyError
            )

        self.setup_residues(residue_ids)

        self.bonded_topology = bonded_topology
        if self.bonded_topology is not None:
            self.logger.warning(
                "There is no check yet if the bonded topology is compatible "
                "with the topology. Please make sure that the bonded topology "
                "is compatible with the topology!"
            )

    def setup_residues(self, residue_ids: Np1DIntArray) -> None:
        """
        A method to set up the residues of the topology.

        Parameters
        ----------
        residue_ids : Np1DIntArray
            The residue ids of the topology.

        Raises
        ------
        ResidueError
            If the residue ids are not contiguous.
        ResidueError
            If the reference residues are not empty and residue_ids
            with 0 don't have any element information. This problem
            can be avoided by setting 'check_residues' to False.
        """

        self._residue_ids = residue_ids
        self._residues, self._atoms = self._setup_residues(
            self.residue_ids,
            self.atoms
        )

        if not self.residues:
            self._residue_numbers = np.arange(self.n_atoms)
            self._residue_atom_indices = [
                np.arange(i, i + 1) for i in range(self.n_atoms)
            ]
        else:
            _residue_numbers = []
            self._residue_atom_indices = []
            atom_counter = 0
            for i in range(self.n_residues):
                _residue_numbers += [i] * self.residues[i].n_atoms
                self._residue_atom_indices.append(
                    np.arange(
                        atom_counter,
                        atom_counter + self.residues[i].n_atoms,
                    )
                )

                atom_counter += self.residues[i].n_atoms
            self._residue_numbers = np.array(_residue_numbers)

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

        is_equal &= self.atoms == other.atoms
        is_equal &= np.all(self.residue_ids == other.residue_ids)

        return bool(is_equal)

    def __getitem__(self, indices: Np1DIntArray | int) -> "Topology":
        """
        Returns a new Topology with the given indices.

        Parameters
        ----------
        indices : Np1DIntArray | int
            The indices of the atoms to return.

        Returns
        -------
        Topology
            The new Topology with the given indices.
        """
        reference_residues = self.reference_residues

        if isinstance(indices, int):
            indices = np.array([indices])

        if len(reference_residues) == 0:
            reference_residues = None

        if self.n_atoms == 0:
            return Topology(reference_residues=reference_residues)

        atoms = [self.atoms[index] for index in indices]
        residue_ids = self.residue_ids[indices]

        return Topology(
            atoms=atoms,
            reference_residues=self.reference_residues,
            residue_ids=residue_ids,
            check_residues=self.check_residues
        )

    def get_atom_indices_from_residue_names(
        self,
        residue_name: str,
    ) -> Np1DIntArray:
        """
        Returns the atom indices for the given residue name.

        Parameters
        ----------
        residue_name : str
            The name of the residue to get the atom indices for.

        Returns
        -------
        Np1DIntArray
            The atom indices for the given residue name.
        """
        atom_indices = np.array([], dtype=int)
        for residue, _atom_indices in zip(self.residues,
                                          self.residue_atom_indices):
            if residue.name.lower() == residue_name.lower():
                atom_indices = np.append(atom_indices, _atom_indices)

        return atom_indices

    def get_atom_indices_from_residue_numbers(
        self,
        residue_numbers: Np1DIntArray,
    ) -> Np1DIntArray:
        """
        Returns the atom indices for the given residue numbers.

        Parameters
        ----------
        residue_numbers : Np1DIntArray
            The residue numbers to get the atom indices for.

        Returns
        -------
        Np1DIntArray
            The atom indices for the given residue numbers.
        """
        return np.argwhere(np.isin(
            self.residue_numbers,
            residue_numbers,
        )).flatten()

    def _setup_residues(
        self,
        residue_ids: Np1DIntArray,
        atoms: Atoms,
    ) -> Tuple[Residues, Atoms]:
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
            If the reference residues are not empty and residue_ids with 0
            don't have any element information. This problem can be avoided
            by setting 'check_residues' to False.
        """
        residues = []

        if len(self.reference_residues) == 0 or not self.check_residues:
            return residues, atoms

        bool_array = [
            id in self.reference_residue_ids or id == 0 for id in residue_ids
        ]

        if not np.all(bool_array):
            not_found_residue_ids = np.unique(
                residue_ids[np.argwhere(~np.array(bool_array))]
            )
            self.logger.error(
                (
                    f"Residue ids {not_found_residue_ids} "
                    "have no corresponding reference residue."
                ),
                exception=ResidueError
            )

        atom_counter = 0
        while atom_counter < len(residue_ids):
            if residue_ids[atom_counter] == 0:
                if atoms[atom_counter].element == Element():
                    message = f"""
The element of atom {atom_counter} is not set. If any reference residues are given
the program tries to automatically deduce the residues from the residue ids and the reference residues.
This means that any atom with an unknown element raises an error. To avoid deducing residue information
please set 'check_residues' to False"""

                    self.logger.error(message, exception=ResidueError)
                else:
                    residues.append(QMResidue(atoms[atom_counter].element))
                    atom_counter += 1
                    continue

            residue = _find_residue_by_id(
                residue_ids[atom_counter],
                self.reference_residues,
            )

            for i in np.arange(residue.n_atoms) + atom_counter:
                if (
                    atoms[i].element != Element() and
                    atoms[i].element != residue.elements[i - atom_counter]
                ):
                    self.logger.warning(
                        (
                            f"The element of atom {i} ({atoms[i].element}) "
                            "does not match the element of the reference residue "
                            f"{residue.name} "
                            f"({residue.elements[i - atom_counter]}). "
                            "Therefore the element type of the residue "
                            "description will be used within the topology format!"
                        )
                    )

                    atoms[i].element = residue.elements[i - atom_counter]

                if residue_ids[i] != residue_ids[atom_counter]:
                    self.logger.error(
                        (
                            "The residue ids are not contiguous. Problems with residue "
                            f"{residue.name} with indices {atom_counter}-"
                            f"{atom_counter + residue.n_atoms-1}"
                        ),
                        exception=ResidueError
                    )

            residues.append(residue)

            atom_counter += residue.n_atoms

        return residues, atoms

    def __str__(self) -> str:
        """
        Returns a string representation of the Topology.

        Returns
        -------
        str
            The string representation of the Topology.
        """

        message = f"Topology with {self.n_atoms} atoms "
        message += f"and {self.n_residues} residues "
        message += f"({self.n_qm_residues} QM residues) " if self.n_qm_residues > 0 else ""
        message += f"and {self.n_unique_residues} unique residues."

        return message

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
        """bool: Whether the residues should be checked."""
        return self._check_residues

    @check_residues.setter
    @runtime_type_checking_setter
    def check_residues(self, value: bool) -> None:
        self._check_residues = value
        self._residues, self._atoms = self._setup_residues(
            self.residue_ids,
            self.atoms,
        )

    @property
    def reference_residue_ids(self) -> Np1DIntArray:
        """Np1DIntArray: The residue ids of the reference residues."""
        return np.array([residue.id for residue in self.reference_residues])

    @property
    def reference_residues(self) -> Residues:
        """Residues: The reference residues of the topology."""
        return self._reference_residues

    @reference_residues.setter
    @runtime_type_checking_setter
    def reference_residues(self, value: Residues):
        self._reference_residues = value

    @property
    def atoms(self) -> Atoms:
        """Atoms: The atoms of the topology."""
        return self._atoms

    @property
    def atomtype_names(self) -> List[str]:
        """List[str]: The atomtype names of the topology."""
        return self._atomtype_names

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms in the topology."""
        return len(self.atoms)

    @property
    def residue_ids(self) -> Np1DIntArray:
        """Np1DIntArray: The residue ids of the topology."""
        return self._residue_ids

    @property
    def residues(self) -> Residues:
        """Residues: The residues of the topology."""
        return self._residues

    @property
    def n_residues(self) -> int:
        """int: The number of residues in the topology."""
        return len(self.residues)

    @property
    def n_qm_residues(self) -> int:
        """int: The number of QM residues in the topology."""
        return len(
            [
                residue for residue in self.residues
                if isinstance(residue, QMResidue)
            ]
        )

    @property
    def n_mm_residues(self) -> int:
        """int: The number of MM residues in the topology."""
        return self.n_residues - self.n_qm_residues

    @property
    def n_unique_residues(self) -> int:
        """int: The number of unique residues in the topology."""
        return len(_unique_residues_(self.residues))

    @property
    def residue_numbers(self) -> Np1DIntArray:
        """Np1DIntArray: The residue numbers of the topology."""
        return self._residue_numbers

    @property
    def residue_atom_indices(self) -> List[Np1DIntArray]:
        """List[Np1DIntArray]: The residue atom indices of the topology."""
        return self._residue_atom_indices

    @property
    def n_atoms_per_residue(self) -> Np1DIntArray:
        """Np1DIntArray: The number of atoms per residue."""
        return np.array(
            [len(indices) for indices in self.residue_atom_indices]
        )

    @property
    def residue_ids_per_residue(self) -> Np1DIntArray:
        """List[Np1DIntArray]: The residue ids per residue."""
        residue_ids_per_residue = []

        if len(self.residue_ids) == 0:
            return residue_ids_per_residue

        for i in np.cumsum(self.n_atoms_per_residue):
            residue_ids_per_residue.append(self.residue_ids[i - 1])

        return np.array(residue_ids_per_residue)



def _find_residue_by_id(res_id: Integral, residues: Residues) -> Residue:
    """
    Finds a residue by its id.

    Parameters
    ----------
    id : Integral
        The id of the residue to find.
    residues : Residues
        The residues to search in - a list of residues.

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
    bool_array = np.array([residue.id == res_id for residue in residues])

    residues = np.array(residues)
    residue = residues[np.argwhere(bool_array)].flatten()

    if len(residue) > 1:
        Topology.logger.error(
            f"The residue id {res_id} is not unique.", exception=ResidueError
        )

    if len(residue) == 0:
        Topology.logger.error(
            f"The residue id {res_id} was not found.", exception=ResidueError
        )

    return residue[0]



def _unique_residues_(residues: Residues) -> Residues:
    """
    Returns a list of unique residues.

    Parameters
    ----------
    residues : Residues
        The residues to make unique.

    Returns
    -------
    Residues
        A list of unique residues.
    """
    unique_residues = []

    for residue in residues:
        if residue not in unique_residues:
            unique_residues.append(residue)

    return unique_residues
