"""
A module containing the Residue class

Classes
-------
Residue
    A class for representing a residue.
QMResidue
    A class for representing a QM residue.
"""

import logging

# library imports
from numbers import Real

# 3rd party object imports
import numpy as np

from beartype.typing import Any, List, TypeVar

# local imports
from PQAnalysis.types import Np1DIntArray, Np1DNumberArray
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking_setter, runtime_type_checking

from .atom import Elements, Element
from .exceptions import ResidueError

#: A type hint for a list of residues (mol types).
Residues = TypeVar("Residues", bound=List["PQAnalysis.core.Residue"])



class Residue:

    """
    A class for representing a residue type (mol type).

    In general residues are used to represent the different molecules in a system.
    In case a QMCF based simulation is performed, the residues are defined via 
    the moldescriptor file.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        name: str,
        residue_id: int,
        total_charge: Real,
        elements: Element | Elements | str | List[str],
        atom_types: int | Np1DIntArray,
        partial_charges: Real | Np1DNumberArray,
    ) -> None:
        """
        Initializes the Residue with the given parameters.

        Parameters
        ----------
        name : str
            The name of the residue.
        residue_id : int
            The id of the residue.
        total_charge : Real
            The total charge of the residue.
        elements : Element | Elements | str | List[str]
            The elements of the residue.
        atom_types : Np1DIntArray
            The atom types of the residue.
        partial_charges : Np1DNumberArray
            The partial charges of the residue.

        Raises
        ------
        MolTypeError
            If the number of elements, atom_types and partial_charges are not the same.
        """

        self.name = name
        self.id = residue_id
        self.total_charge = total_charge

        # set here the internal variables to avoid setters
        # (which would check the length of the arrays)
        if isinstance(elements, Element):
            self._elements = [elements]
        elif isinstance(elements, str):
            self._elements = [Element(elements)]
        elif isinstance(elements,
            list) and len(elements) > 0 and isinstance(elements[0],
            Element):
            self._elements = elements
        else:
            self._elements = [Element(element) for element in elements]

        self._atom_types = np.atleast_1d(atom_types)
        self._partial_charges = np.atleast_1d(partial_charges)

        if not len(self.elements) == len(self.atom_types) == len(
                self.partial_charges):
            self.logger.error(
                "The number of elements, atom_types and partial_charges must be the same.",
                exception=ResidueError
            )

    @property
    def n_atoms(self) -> int:
        """
        Returns the number of atoms in the residue.

        Returns
        -------
        int
            The number of atoms in the residue.
        """
        return len(self._elements)

    @property
    def name(self) -> str:
        """
        Returns the name of the residue.

        Returns
        -------
        str
            The name of the residue.
        """
        return self._name

    @name.setter
    @runtime_type_checking_setter
    def name(self, name: str) -> None:
        """
        Sets the name of the residue.

        Parameters
        ----------
        name : str
            The name of the residue.
        """
        self._name = name

    @property
    def id(self) -> int:
        """
        Returns the id of the residue.

        Returns
        -------
        int
            The id of the residue.
        """
        return self._id

    @id.setter
    @runtime_type_checking_setter
    def id(self, residue_id: int) -> None:
        """
        Sets the id of the residue.

        Parameters
        ----------
        residue_id : int
            The id of the residue.
        """
        self._id = residue_id

    @property
    def total_charge(self) -> Real:
        """
        Returns the total charge of the residue.

        Returns
        -------
        Real
            The total charge of the residue.
        """
        return self._total_charge

    @total_charge.setter
    @runtime_type_checking_setter
    def total_charge(self, total_charge: Real) -> None:
        """
        Sets the total charge of the residue.

        Parameters
        ----------
        total_charge : Real
            The total charge of the residue.
        """
        self._total_charge = total_charge

    @property
    def elements(self) -> Elements:
        """
        Returns the elements of the residue.

        Returns
        -------
        Elements
            The elements of the residue.
        """
        return self._elements

    @elements.setter
    @runtime_type_checking_setter
    def elements(self, elements: Elements) -> None:
        """
        Sets the elements of the residue.

        Parameters
        ----------
        elements : Elements
            The elements of the residue.

        Raises
        ------
        ResidueError
            If the number of elements is not the same as the number of atoms.
        """
        if len(elements) != self.n_atoms:
            self.logger.error(
                "The number of elements must be the same as the number of atoms.",
                exception=ResidueError
            )

        self._elements = elements

    @property
    def atom_types(self) -> Np1DIntArray:
        """
        Returns the atom types of the residue.

        Returns
        -------
        Np1DIntArray
            The atom types of the residue.
        """
        return self._atom_types

    @atom_types.setter
    @runtime_type_checking_setter
    def atom_types(self, atom_types: Np1DIntArray) -> None:
        """
        Sets the atom types of the residue.

        Parameters
        ----------
        atom_types : Np1DIntArray
            The atom types of the residue.

        Raises
        ------
        ResidueError
            If the number of atom_types is not the same as the number of atoms.
        """
        if len(atom_types) != self.n_atoms:
            self.logger.error(
                "The number of atom_types must be the same as the number of atoms.",
                exception=ResidueError
            )

        self._atom_types = atom_types

    @property
    def partial_charges(self) -> Np1DNumberArray:
        """
        Returns the partial charges of the residue.

        Returns
        -------
        Np1DNumberArray
            The partial charges of the residue.
        """
        return self._partial_charges

    @partial_charges.setter
    @runtime_type_checking_setter
    def partial_charges(self, partial_charges: Np1DNumberArray) -> None:
        """
        Sets the partial charges of the residue.

        Parameters
        ----------
        partial_charges : Np1DNumberArray
            The partial charges of the residue.

        Raises
        ------
        ResidueError
            If the number of partial_charges is not the same as the number of atoms.
        """
        if len(partial_charges) != self.n_atoms:
            self.logger.error(
                "The number of partial_charges must be the same as the number of atoms.",
                exception=ResidueError
            )

        self._partial_charges = partial_charges

    def __str__(self) -> str:
        """
        Returns the string representation of the Residue.

        Returns
        -------
        str
            The string representation of the Residue.
        """

        name = self.name
        total_charge = self.total_charge
        n_atoms = self.n_atoms

        return f"Residue({name=}, id={self.id}, {total_charge=}, {n_atoms=})"

    def __repr__(self) -> str:
        """
        Returns the string representation of the Residue.

        Returns
        -------
        str
            The string representation of the Residue.
        """
        return self.__str__()

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the Residue is equal to another Residue.

        Parameters
        ----------
        other : Any
            The other Residue to compare to.

        Returns
        -------
        bool
            True if the Residue is equal to the other Residue, False otherwise.
        """

        is_equal = True
        is_equal &= self.name.lower() == other.name.lower()
        is_equal &= self.id == other.id
        is_equal &= np.allclose(self.total_charge, other.total_charge)
        is_equal &= self.elements == other.elements
        is_equal &= np.allclose(self.atom_types, other.atom_types)
        is_equal &= np.allclose(self.partial_charges, other.partial_charges)
        return is_equal



class QMResidue(Residue):

    """
    A class for representing a QM residue type (mol type).

    It is a subclass of Residue and is used to represent any atoms that are
    represented solely by via QM methods during a simulation, meaning that
    they cannot leave the QM zone and do not have any MM-like properties.

    Examples
    --------
    >>> QMResidue(Element('C'))
    Residue(name='QM', id=0, total_charge=0.0, n_atoms=1)

    """

    @runtime_type_checking
    def __init__(self, element: Element | str) -> None:
        """
        Initializes the QMResidue with the given parameters.

        For the initialization only the element of the QM residue is required.
        One QMResidue can only contain one element.

        Parameters
        ----------
        element : Element
            The element of the QM residue.
        """
        super().__init__(
            name="QM",
            residue_id=0,
            total_charge=0.0,
            elements=[element],
            atom_types=np.array([0]),
            partial_charges=np.array([element.atomic_number])
        )
