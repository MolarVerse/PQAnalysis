"""
A module containing the Element class.
"""

import logging

from typing import Annotated
from numbers import Real

from beartype.typing import Any, NewType
from beartype.vale import Is

from PQAnalysis.type_checking import (
    runtime_type_checking,
    runtime_type_checking_setter,
)
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from ..exceptions import ElementNotFoundError



class Element:

    """
    A class representing an element.

    An Element can be initialized in three ways:

    1) By giving the atomic number of the element (e.g. 6).
    2) By giving the symbol of the element (e.g. 'C').
    3) By giving None.

    Examples
    --------

    >>> element = Element(6)
    >>> (element.symbol, element.atomic_number, element.mass)
    ('c', 6, 12.0107)

    >>> element = Element('C')
    >>> (element.symbol, element.atomic_number, element.mass)
    ('c', 6, 12.0107)

    >>> element = Element()
    >>> (element.symbol, element.atomic_number, element.mass)
    (None, None, None)

    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, element_id: int | str | None = None) -> None:
        """
        Initializes the Element with the given parameters.

        Parameters
        ----------
        element_id : int | str | None, optional
            The identifier of the element. If an integer is given, it is
            interpreted as the atomic number. If a string is given, it is
            interpreted as the symbol of the element. If None is given, the
            symbol, atomic number and mass are set to None, by default None

        Raises
        ------
        ElementNotFoundError
            If the given id is not a valid element identifier.
        """

        try:
            # If id is an integer, it is interpreted as the atomic number.
            if isinstance(element_id, int):

                self._atomic_number = element_id
                self._symbol = atomicNumbersReverse[element_id]

            # If id is a string, it is interpreted as the symbol of the element.
            elif isinstance(element_id, str):

                self._symbol = element_id.lower()
                self._atomic_number = atomicNumbers[self._symbol]

            # If id is None, the symbol, atomic number
            # and mass are set to None, meaning an empty element.
            if element_id is None:

                self._symbol = None
                self._atomic_number = None
                self._mass = None

            else:

                self._mass = atomicMasses[self._symbol]

        except KeyError:
            self.logger.error(
                ElementNotFoundError(element_id).message,
                exception=ElementNotFoundError
            )

    def __str__(self) -> str:
        """
        Returns a string representation of the Element.

        Returns
        -------
        str
            A string representation of the Element.
        """
        return f"Element({self.symbol}, {self.atomic_number}, {self.mass})"

    def __repr__(self) -> str:
        """
        Returns a representation of the Element.

        Returns
        -------
        str
            A representation of the Element.
        """
        return self.__str__()

    def __eq__(self, other: Any) -> bool:
        """
        Checks whether the Element is equal to another Element.

        Parameters
        ----------
        other : Any
            The other Element to check for equality.

        Returns
        -------
        bool
            True if the Element is equal to the other Element, False otherwise.
        """

        if not isinstance(other, Element):
            return False

        return self.symbol == other.symbol and self.atomic_number == other.atomic_number

    @property
    def symbol(self) -> str | None:
        """str | None: The symbol of the Element."""
        return self._symbol

    @property
    def atomic_number(self) -> int | None:
        """int | None: The atomic number of the Element."""
        return self._atomic_number

    @property
    def mass(self) -> Real | None:
        """Real | None: The mass of the Element."""
        return self._mass



class CustomElement(Element):

    """
    A class representing a custom element it 
    inherits from the Element class.
    """

    @runtime_type_checking
    def __init__(self, symbol: str, atomic_number: int, mass: Real):  # pylint: disable=super-init-not-called
        """
        Parameters
        ----------
        symbol : str
            the custom atomic symbol
        atomic_number : int
            the custom atomic number
        mass : Real
            the custom atomic mass
        """

        self._symbol = symbol
        self._atomic_number = atomic_number
        self._mass = mass

    @Element.symbol.setter
    @runtime_type_checking_setter
    def symbol(self, value: str) -> None:
        self._symbol = value



#: A type hint for a list of elements
Elements = NewType(
    "Elements",
    Annotated[
        list,
        Is[lambda list: all(isinstance(element, Element) for element in list)]]
)

atomicMasses = {
    "h": 1.00794,
    "d": 2.014101778,
    "t": 3.0160492675,
    "he": 4.002602,
    "li": 6.941,
    "be": 9.012182,
    "b": 10.811,
    "c": 12.0107,
    "n": 14.0067,
    "o": 15.9994,
    "f": 18.9984032,
    "ne": 20.1797,
    "na": 22.989770,
    "mg": 24.3050,
    "al": 26.981538,
    "si": 28.0855,
    "p": 30.973761,
    "s": 32.065,
    "cl": 35.453,
    "ar": 39.948,
    "k": 39.0983,
    "ca": 40.078,
    "sc": 44.955910,
    "ti": 47.880,
    "v": 50.9415,
    "cr": 51.9961,
    "mn": 54.938049,
    "fe": 55.845,
    "co": 58.933200,
    "ni": 58.6934,
    "cu": 63.546,
    "zn": 65.399,
    "ga": 69.723,
    "ge": 72.64,
    "as": 74.92160,
    "se": 78.96,
    "br": 79.904,
    "kr": 83.798,
    "rb": 85.4678,
    "sr": 87.62,
    "y": 88.90585,
    "zr": 91.224,
    "nb": 92.90638,
    "mo": 95.94,
    "tc": 98.9063,
    "ru": 101.07,
    "rh": 102.9055,
    "pd": 106.42,
    "ag": 107.8682,
    "cd": 112.411,
    "in": 114.818,
    "sn": 118.71,
    "sb": 121.76,
    "te": 127.6,
    "i": 126.90447,
    "xe": 131.293,
    "cs": 132.90546,
    "ba": 137.327,
    "la": 138.9055,
    "ce": 140.116,
    "pr": 140.90765,
    "nd": 144.24,
    "pm": 146.9151,
    "sm": 150.36,
    "eu": 151.964,
    "gd": 157.25,
    "tb": 158.92534,
    "dy": 162.5,
    "ho": 164.93032,
    "er": 167.259,
    "tm": 168.93421,
    "yb": 173.04,
    "lu": 174.967,
    "hf": 178.49,
    "ta": 180.9479,
    "w": 183.84,
    "re": 186.207,
    "os": 190.23,
    "ir": 192.217,
    "pt": 195.078,
    "au": 196.96655,
    "hg": 200.59,
    "tl": 204.3833,
    "pb": 207.2,
    "bi": 208.98038,
    "po": 208.9824,
    "at": 209.9871,
    "rn": 222.0176,
    "fr": 223.0197,
    "ra": 226.0254,
    "ac": 227.0278,
    "th": 232.0381,
    "pa": 231.03588,
    "u": 238.0289,
    "np": 237.0482,
    "pu": 244.0642,
    "am": 243.0614,
    "cm": 247.0703,
    "bk": 247.0703,
    "cf": 251.0796,
    "es": 252.0829,
    "fm": 257.0951,
    "md": 258.0986,
    "no": 259.1009,
    "lr": 260.1053,
    "q": 999.00000,
    "x": 999.00000,
    "cav": 1000.00000,
    "sup": 1000000.0,
    "dum": 1.0
}

atomicNumbers = {
    "h": 1,
    "d": 1,
    "t": 1,
    "he": 2,
    "li": 3,
    "be": 4,
    "b": 5,
    "c": 6,
    "n": 7,
    "o": 8,
    "f": 9,
    "ne": 10,
    "na": 11,
    "mg": 12,
    "al": 13,
    "si": 14,
    "p": 15,
    "s": 16,
    "cl": 17,
    "ar": 18,
    "k": 19,
    "ca": 20,
    "sc": 21,
    "ti": 22,
    "v": 23,
    "cr": 24,
    "mn": 25,
    "fe": 26,
    "co": 27,
    "ni": 28,
    "cu": 29,
    "zn": 30,
    "ga": 31,
    "ge": 32,
    "as": 33,
    "se": 34,
    "br": 35,
    "kr": 36,
    "rb": 37,
    "sr": 38,
    "y": 39,
    "zr": 40,
    "nb": 41,
    "mo": 42,
    "tc": 43,
    "ru": 44,
    "rh": 45,
    "pd": 46,
    "ag": 47,
    "cd": 48,
    "in": 49,
    "sn": 50,
    "sb": 51,
    "te": 52,
    "i": 53,
    "xe": 54,
    "cs": 55,
    "ba": 56,
    "la": 57,
    "ce": 58,
    "pr": 59,
    "nd": 60,
    "pm": 61,
    "sm": 62,
    "eu": 63,
    "gd": 64,
    "tb": 65,
    "dy": 66,
    "ho": 67,
    "er": 68,
    "tm": 69,
    "yb": 70,
    "lu": 71,
    "hf": 72,
    "ta": 73,
    "w": 74,
    "re": 75,
    "os": 76,
    "ir": 77,
    "pt": 78,
    "au": 79,
    "hg": 80,
    "tl": 81,
    "pb": 82,
    "bi": 83,
    "po": 84,
    "at": 85,
    "rn": 86,
    "fr": 87,
    "ra": 88,
    "ac": 89,
    "th": 90,
    "pa": 91,
    "u": 92,
    "np": 93,
    "pu": 94,
    "am": 95,
    "cm": 96,
    "bk": 97,
    "cf": 98,
    "es": 99,
    "fm": 100,
    "md": 101,
    "no": 102,
    "lr": 103,
    "q": 999,
    "x": 999,
    "cav": 1000,
    "sup": 1000000,
    "dum": 1
}

atomicNumbersReverse = {v: k for k, v in atomicNumbers.items()}
