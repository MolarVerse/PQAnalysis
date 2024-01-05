"""
A module containing the Atom class and some associated functions.

The Atom class is used to represent an atom in a molecule. It contains the
following attributes:
    
        - name: the name of the atom (e.g. 'C')
        - element: the element of the atom (e.g. 'C')

The atomic number and mass are automatically determined from the name or symbol
of the atom. The name and symbol are automatically determined from the atomic
number. The name is the symbol in lower case.

...

Classes
-------
Atom
    A class used to represent an atom in a molecule.
"""

from __future__ import annotations

from multimethod import multimeta
from beartype.typing import Any
from beartype.vale import Is
from typing import Annotated
from numbers import Real

from . import Element

"""
A type hint for a list of Atom objects.
"""
Atoms = Annotated[list, Is[lambda list: all(
    isinstance(atom, Atom) for atom in list)]]


class Atom(metaclass=multimeta):
    """
    A class used to represent an atom in a molecule.

    There are three ways to initialize an Atom object:

    1) By giving the name of the atom_type (e.g. 'C1')
       If use_guess_element is True (default), the atom_type name has to be
       a valid element symbol (e.g. 'C'). If use_guess_element is False, the
       atom_type name can be anything and an empty element is created.

    2) By giving the name of the atom_type (e.g. 'C1') and the id of the atom_type
       (e.g. 6). The id can be either an integer (atomic number) or a string (element symbol).

    3) By giving the id of the atom_type (e.g. 6). The id can be either an integer (atomic number) 
       or a string (element symbol).

    Examples
    --------
    >>> atom = Atom('C1') # use_guess_element is True by default - raises ElementNotFoundError if the element is not found

    >>> atom = Atom('C1', use_guess_element=False)
    >>> (atom.name, atom.element)
    ('C1', Element())

    >>> atom = Atom('C1', 'C')
    >>> (atom.name, atom.element)
    ('C1', Element('C'))

    >>> atom = Atom('C1', 6)
    >>> (atom.name, atom.element)
    ('C1', Element(6))

    >>> atom = Atom(6)
    >>> (atom.name, atom.element)
    ('c', Element(6))

    >>> atom = Atom('C')
    >>> (atom.name, atom.element)
    ('C', Element('C'))


    """

    def __init__(self, name: str, use_guess_element: bool = True) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        If use_guess_element is True, the symbol, atomic number and mass are
        determined from the name of the atom_type. If use_guess_element is
        False, the symbol, atomic number and mass are set to None, meaning that
        an empty element is created (Element()).

        Parameters
        ----------
        name : str
            The name of the atom_type (e.g. 'C1')
        use_guess_element : bool, optional
            Whether to use the guess_element function to determine the element type of the atom_type 
            by its name, by default True
        """
        self._name = name
        if use_guess_element:
            self._element = Element(name)
        else:
            self._element = Element()

    def __init__(self, name: str, id: int | str) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        Parameters
        ----------
        name : str
            The name of the atom_type (e.g. 'C1')
        id : int | str
            The identifier of the atom_type. If an integer is given, it is
            assumed to be the atomic number of the atom_type. If a string is
            given, it is assumed to be the name of the element.
        """
        self._name = name
        self._element = Element(id)

    def __init__(self, id: int) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        Parameters
        ----------
        id : int
            The identifier of the atom_type. It is
            assumed to be the atomic number of the atom_type.
        """
        self._element = Element(id)
        self._name = self._element.symbol.lower()

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the Atom is equal to another Atom.

        Parameters
        ----------
        other : Any
            The other Atom to compare to.

        Returns
        -------
        bool
            True if the Atom is equal to the other Atom, False otherwise.
        """
        if not isinstance(other, Atom):
            return False

        is_equal = True
        is_equal &= self.name.lower() == other.name.lower()
        is_equal &= self._element == other._element
        return is_equal

    def __str__(self) -> str:
        """
        Returns a string representation of the Atom.

        Returns
        -------
        str
            A string representation of the Atom.
        """
        return f"Atom({self.name}, {self.atomic_number}, {self.symbol}, {self.mass})"

    def __repr__(self) -> str:
        """
        Returns a representation of the Atom.

        Returns
        -------
        str
            A representation of the Atom.
        """
        return self.__str__()

    #######################
    #                     #
    # standard properties #
    #                     #
    #######################

    @property
    def name(self) -> str:
        """
        The name of the atom_type (e.g. 'C1')

        Returns
        -------
        str
            The name of the atom_type (e.g. 'C1')
        """
        return self._name

    @property
    def symbol(self) -> str | None:
        """
        The symbol of the element (e.g. 'c')

        Returns
        -------
        str
            The symbol of the element (e.g. 'c')
        """
        return self._element.symbol

    @property
    def atomic_number(self) -> int | None:
        """
        The atomic number of the element (e.g. 6)

        Returns
        -------
        int
            The atomic number of the element (e.g. 6)
        """
        return self._element.atomic_number

    @property
    def mass(self) -> Real | None:
        """
        The mass of the element (e.g. 12.0107)

        Returns
        -------
        Real
            The mass of the element (e.g. 12.0107)
        """
        return self._element.mass

    @property
    def element(self) -> Element:
        """
        The element of the atom

        Returns
        -------
        Element
            The element of the atom
        """
        return self._element

    @element.setter
    def element(self, element: Element) -> None:
        """
        Sets the element of the atom

        Parameters
        ----------
        element : Element
            The element of the atom
        """
        self._element = element
