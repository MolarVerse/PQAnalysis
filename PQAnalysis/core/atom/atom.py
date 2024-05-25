"""
A module containing the Atom class and some associated functions.

The Atom class is used to represent an atom in a molecule. It contains the
following attributes:
    
        - name: the name of the atom (e.g. 'C')
        - element: the element of the atom (e.g. 'C')

The atomic number and mass are automatically determined from the name or symbol
of the atom. The name and symbol are automatically determined from the atomic
number. The name is the symbol in lower case.
"""

import logging

from numbers import Real

from beartype.typing import Any, List, TypeVar

from PQAnalysis import __package_name__
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.exceptions import PQValueError

from .element import Element
from ..exceptions import AtomError

#: A type hint for a list of atoms
Atoms = TypeVar("Atoms", bound=List["PQAnalysis.core.Atom"])



class Atom():

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
    >>> atom = Atom('C1') # use_guess_element is True by default
    Traceback (most recent call last):
        ...
    ElementNotFoundError: Id C1 is not a valid element identifier.

    >>> atom = Atom('C1', use_guess_element=False)
    >>> (atom.name, atom.element)
    ('C1', Element(None, None, None))

    >>> atom = Atom('C1', 'C')
    >>> (atom.name, atom.element)
    ('C1', Element(c, 6, 12.0107))

    >>> atom = Atom('C1', 6)
    >>> (atom.name, atom.element)
    ('C1', Element(c, 6, 12.0107))

    >>> atom = Atom(6)
    >>> (atom.name, atom.element)
    ('c', Element(c, 6, 12.0107))

    >>> atom = Atom('C')
    >>> (atom.name, atom.element)
    ('C', Element(c, 6, 12.0107))
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        name: str | int | Element,
        element_id: int | str | None = None,
        use_guess_element: bool = True,
        **_kwargs
    ) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        If use_guess_element is True, the symbol, atomic number and mass are
        determined from the name of the atom_type. If use_guess_element is
        False, the symbol, atomic number and mass are set to None, meaning that
        an empty element is created (Element()).

        Parameters
        ----------
        name : str | int | Element
            The name of the atom_type (e.g. 'C1')
            If this parameter is an integer, it is interpreted as the atomic number of the 
            element symbol and cannot be used together with the id parameter. If the parameter
            is an Element object, the element is set to this object and no other parameters
            are used.
        element_id : int | str, optional
            The atomic number or element symbol of the atom_type (e.g. 6 or 'C') 
            If his parameter is not given, the name parameter is used to determine 
            the element type of the atom_type.
        use_guess_element : bool, optional
            Whether to use the guess_element function to determine
            the element type of the atom_type by its name,
            by default True
        **_kwargs
            Additional keyword arguments that are not used by the function but 
            by the runtime type checker.

        Raises
        ------
        PQValueError
            If the name of the atom_type is an Element object and the id is given.
        AtomError
            If the name of the atom_type is an integer and the id is given.
        """

        if isinstance(name, Element):
            if element_id is not None:
                self.logger.error(
                    (
                        "The element of the atom_type cannot be an "
                        "Element object if the id is given."
                    ),
                    exception=PQValueError
                )

            self._name = name.symbol
            self._element = name
            return

        if element_id is not None and isinstance(name, int):

            self.logger.error(
                "The name of the atom_type cannot be an integer if the id is given.",
                exception=AtomError
            )

        if element_id is not None and isinstance(name, str):

            self._name = name
            self._element = Element(element_id)

        elif isinstance(name, int):

            self._element = Element(name)
            self._name = self._element.symbol.lower()

        else:

            self._name = name
            if use_guess_element:
                self._element = Element(name)
            else:
                self._element = Element()

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
        """str: The name of the atom_type (e.g. 'C1')"""
        return self._name

    @property
    def symbol(self) -> str | None:
        """str: The symbol of the element (e.g. 'c')"""
        return self._element.symbol

    @property
    def atomic_number(self) -> int | None:
        """int | None: The atomic number of the element (e.g. 6)"""
        return self._element.atomic_number

    @property
    def mass(self) -> Real | None:
        """Real | None: The mass of the element (e.g. 12.011)"""
        return self._element.mass

    @property
    def element(self) -> Element:
        """Element: The element type of the atom"""
        return self._element

    @element.setter
    def element(self, element: Element) -> None:
        self._element = element

    @property
    def element_name(self) -> str:
        """str: The name of the element (e.g. 'Carbon')"""
        return self._element.symbol
