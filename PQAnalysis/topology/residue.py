"""
A module containing the Residue class

Classes
-------
Residue
    A class for representing a residue.
QMResidue
    A class for representing a QM residue.
"""

from __future__ import annotations

# library imports
import numpy as np

# 3rd party object imports
from numbers import Real
from beartype.typing import List
from beartype.vale import Is
from typing import Annotated

# local imports
from .exceptions import ResidueError
from ..types import Np1DIntArray, Np1DNumberArray
from ..core import Elements, Element

"""
A type hint for a list of residues (mol types).
"""
Residues = Annotated[list, Is[lambda list: all(
    isinstance(residue, Residue) for residue in list)]]


class Residue:
    """
    A class for representing a residue type (mol type).

    In general residues are used to represent the different molecules in a system.
    In case a QMCF based simulation is performed, the residues are defined via 
    the moldescriptor file.

    Attributes
    ----------
    name : str
        The name of the residue.
    id : int
        The id of the residue.
    total_charge : Real
        The total charge of the residue.
    elements : Elements
        The elements of the residue.
    atom_types : Np1DIntArray
        The atom types of the residue.
    partial_charges : Np1DNumberArray
        The partial charges of the residue.
    """

    def __init__(self,
                 name: str,
                 id: int,
                 total_charge: Real,
                 elements: Elements,
                 atom_types: Np1DIntArray,
                 partial_charges: Np1DNumberArray,
                 ) -> None:
        """
        Initializes the Residue with the given parameters.

        Parameters
        ----------
        name : str
            The name of the residue.
        id : int
            The id of the residue.
        total_charge : Real
            The total charge of the residue.
        elements : Elements
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

        if not (len(elements) == len(atom_types) == len(partial_charges)):
            raise ResidueError(
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
    def id(self, id: int) -> None:
        """
        Sets the id of the residue.

        Parameters
        ----------
        id : int
            The id of the residue.
        """
        self._id = id

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
    def elements(self, elements: Elements) -> None:
        """
        Sets the elements of the residue.

        Parameters
        ----------
        elements : Elements
            The elements of the residue.

        Raises
        ------
        MolTypeError
            If the number of elements is not the same as the number of atoms.
        """
        if len(elements) != self.n_atoms:
            raise ResidueError(
                "The number of elements must be the same as the number of atoms.")

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
    def atom_types(self, atom_types: Np1DIntArray) -> None:
        """
        Sets the atom types of the residue.

        Parameters
        ----------
        atom_types : Np1DIntArray
            The atom types of the residue.

        Raises
        ------
        MolTypeError
            If the number of atom_types is not the same as the number of atoms.
        """
        if len(atom_types) != self.n_atoms:
            raise ResidueError(
                "The number of atom_types must be the same as the number of atoms.")

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
    def partial_charges(self, partial_charges: Np1DNumberArray) -> None:
        """
        Sets the partial charges of the residue.

        Parameters
        ----------
        partial_charges : Np1DNumberArray
            The partial charges of the residue.

        Raises
        ------
        MolTypeError
            If the number of partial_charges is not the same as the number of atoms.
        """
        if len(partial_charges) != self.n_atoms:
            raise ResidueError(
                "The number of partial_charges must be the same as the number of atoms.")

        self._partial_charges = partial_charges

    def __str__(self) -> str:
        """
        Returns the string representation of the Residue.

        Returns
        -------
        str
            The string representation of the Residue.
        """
        return f"Residue(name={self.name}, id={self.id}, total_charge={self.total_charge}, n_atoms={self.n_atoms})"

    def __repr__(self) -> str:
        """
        Returns the string representation of the Residue.

        Returns
        -------
        str
            The string representation of the Residue.
        """
        return self.__str__()


class QMResidue(Residue):
    """
    A class for representing a QM residue type (mol type).

    It is a subclass of Residue and is used to represent any atoms that are
    represented solely by via QM methods during a simulation, meaning that
    they cannot leave the QM zone and do not have any MM-like properties.

    Examples
    --------
    >>> QMResidue(Element('C'))

    """

    def __init__(self, element: Element) -> None:
        """
        Initializes the QMResidue with the given parameters.

        For the initialization only the element of the QM residue is required.
        One QMResidue can only contain one element.

        Parameters
        ----------
        element : Element
            The element of the QM residue.
        """
        super().__init__(name="QM", id=0, total_charge=0.0, elements=[
            element], atom_types=np.array([0]), partial_charges=np.array([element.atomic_number]))
