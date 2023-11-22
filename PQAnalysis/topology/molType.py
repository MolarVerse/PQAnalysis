"""
A module containing the MolType class

Classes
-------
MolType
    A class for representing a molecule type.
"""

from numbers import Real
from beartype.typing import List

from ..types import Np1DIntArray, Np1DNumberArray
from ..core.atom import Atom


class MolType:
    """
    A class for representing a molecule type.

    In general mol types are used to represent the different molecules in a system.
    In case a QMCF based simulation is performed, the mol types are defined via 
    the moldescriptor file.

    Attributes
    ----------
    name : str
        The name of the molecule type.
    id : int
        The id of the molecule type.
    total_charge : Real
        The total charge of the molecule type.
    elements : List[Atom]
        The elements of the molecule type.
    atom_types : Np1DIntArray
        The atom types of the molecule type.
    partial_charges : Np1DNumberArray
        The partial charges of the molecule type.
    """

    def __init__(self,
                 name: str,
                 id: int,
                 total_charge: Real,
                 elements: List[Atom],
                 atom_types: Np1DIntArray,
                 partial_charges: Np1DNumberArray,
                 ) -> None:
        """
        Initializes the MolType with the given parameters.

        Parameters
        ----------
        name : str
            The name of the molecule type.
        id : int
            The id of the molecule type.
        total_charge : Real
            The total charge of the molecule type.
        elements : List[Atom]
            The elements of the molecule type.
        atom_types : Np1DIntArray
            The atom types of the molecule type.
        partial_charges : Np1DNumberArray
            The partial charges of the molecule type.

        Raises
        ------
        ValueError
            If the number of elements, atom_types and partial_charges are not the same.
        """

        if not (len(elements) == len(atom_types) == len(partial_charges)):
            raise ValueError(
                "The number of elements, atom_types and partial_charges must be the same.")

        self.name = name
        self.id = id
        self.total_charge = total_charge

        # set here the internal variables to avoid setters
        # (which would check the length of the arrays)
        self._elements = elements
        self._atom_types = atom_types
        self._partial_charges = partial_charges

    @property
    def n_atoms(self) -> int:
        """
        Returns the number of atoms in the molecule type.

        Returns
        -------
        int
            The number of atoms in the molecule type.
        """
        return len(self._elements)

    @property
    def name(self) -> str:
        """
        Returns the name of the molecule type.

        Returns
        -------
        str
            The name of the molecule type.
        """
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        """
        Sets the name of the molecule type.

        Parameters
        ----------
        name : str
            The name of the molecule type.
        """
        self._name = name

    @property
    def id(self) -> int:
        """
        Returns the id of the molecule type.

        Returns
        -------
        int
            The id of the molecule type.
        """
        return self._id

    @id.setter
    def id(self, id: int) -> None:
        """
        Sets the id of the molecule type.

        Parameters
        ----------
        id : int
            The id of the molecule type.
        """
        self._id = id

    @property
    def total_charge(self) -> Real:
        """
        Returns the total charge of the molecule type.

        Returns
        -------
        Real
            The total charge of the molecule type.
        """
        return self._total_charge

    @total_charge.setter
    def total_charge(self, total_charge: Real) -> None:
        """
        Sets the total charge of the molecule type.

        Parameters
        ----------
        total_charge : Real
            The total charge of the molecule type.
        """
        self._total_charge = total_charge

    @property
    def elements(self) -> List[Atom]:
        """
        Returns the elements of the molecule type.

        Returns
        -------
        List[Atom]
            The elements of the molecule type.
        """
        return self._elements

    @elements.setter
    def elements(self, elements: List[Atom]) -> None:
        """
        Sets the elements of the molecule type.

        Parameters
        ----------
        elements : List[Atom]
            The elements of the molecule type.

        Raises
        ------
        ValueError
            If the number of elements is not the same as the number of atoms.
        """
        if len(elements) != self.n_atoms:
            raise ValueError(
                "The number of elements must be the same as the number of atoms.")

        self._elements = elements

    @property
    def atom_types(self) -> Np1DIntArray:
        """
        Returns the atom types of the molecule type.

        Returns
        -------
        Np1DIntArray
            The atom types of the molecule type.
        """
        return self._atom_types

    @atom_types.setter
    def atom_types(self, atom_types: Np1DIntArray) -> None:
        """
        Sets the atom types of the molecule type.

        Parameters
        ----------
        atom_types : Np1DIntArray
            The atom types of the molecule type.

        Raises
        ------
        ValueError
            If the number of atom_types is not the same as the number of atoms.
        """
        if len(atom_types) != self.n_atoms:
            raise ValueError(
                "The number of atom_types must be the same as the number of atoms.")

        self._atom_types = atom_types

    @property
    def partial_charges(self) -> Np1DNumberArray:
        """
        Returns the partial charges of the molecule type.

        Returns
        -------
        Np1DNumberArray
            The partial charges of the molecule type.
        """
        return self._partial_charges

    @partial_charges.setter
    def partial_charges(self, partial_charges: Np1DNumberArray) -> None:
        """
        Sets the partial charges of the molecule type.

        Parameters
        ----------
        partial_charges : Np1DNumberArray
            The partial charges of the molecule type.

        Raises
        ------
        ValueError
            If the number of partial_charges is not the same as the number of atoms.
        """
        if len(partial_charges) != self.n_atoms:
            raise ValueError(
                "The number of partial_charges must be the same as the number of atoms.")

        self._partial_charges = partial_charges
