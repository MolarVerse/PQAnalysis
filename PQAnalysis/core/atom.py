"""
A module containing the Atom class and some associated functions.

The Atom class is used to represent an atom in a molecule. It contains the
following attributes:
    
        - name: the name of the atom (e.g. 'C')
        - symbol: the symbol of the atom (e.g. 'c')
        - atomic_number: the atomic number of the atom (e.g. 6)
        - mass: the mass of the atom (e.g. 12.0107)

The atomic number and mass are automatically determined from the name or symbol
of the atom. The name and symbol are automatically determined from the atomic
number. The name is the symbol in lower case.

...

Classes
-------
Atom
    A class used to represent an atom in a molecule.

Functions
---------
guess_element
    Guesses the symbol, atomic number and mass of an atom from its name or
    atomic number.

"""

import numpy as np

from multimethod import multimethod
from beartype.typing import Any, Tuple
from numbers import Real

from PQAnalysis.utils.exceptions import ElementNotFoundError


def guess_element(id: int | str) -> Tuple[str, int, Real]:
    """
    Guesses the symbol, atomic number and mass of an atom from a string or
    integer identifier.

    Parameters
    ----------
    id : int  |  str
        The identifier of the atom. If an integer is given, it is assumed to be
        the atomic number of the atom. If a string is given, it is assumed to
        be the name of the element.

    Returns
    -------
    symbol : str
        The symbol of the atom (e.g. 'c')
    atomic_number : int
        The atomic number of the atom (e.g. 6)
    mass : Real
        The mass of the atom (e.g. 12.0107)

    Raises
    ------
    ElementNotFoundError
        If the given identifier is not a valid element identifier.
    TypeError
        If the given identifier is not an integer or string.
    """
    if isinstance(id, int):
        try:
            index = list(atomicNumbers.values()).index(id)
            symbol = list(atomicNumbers)[index]
            mass = atomicMasses[symbol]
            return symbol, id, mass
        except Exception:
            raise ElementNotFoundError(id)
    elif isinstance(id, str):
        try:
            symbol = id.lower()
            atomic_number = atomicNumbers[symbol]
            mass = atomicMasses[symbol]
            return symbol, atomic_number, mass
        except Exception:
            raise ElementNotFoundError(id)


class Atom:
    """
    A class used to represent an atom in a molecule.
    """

    @multimethod
    def __init__(self, name: str, use_guess_element: bool = True) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        If use_guess_element is True, the symbol, atomic number and mass are
        determined from the name of the atom_type. If use_guess_element is
        False, the symbol, atomic number and mass are set to None.

        Parameters
        ----------
        name : str
            The name of the atom_type (e.g. 'C1')
        use_guess_element : bool, optional
            Whether to use the guess_element function to determine the symbol,
            atomic number and mass of the atom_type. If True, the symbol,
            atomic number and mass are determined from the name of the
            atom_type. The default is True.
        """
        self._name = name
        if use_guess_element:
            self._symbol, self._atomic_number, self._mass = guess_element(name)
        else:
            self._symbol = None
            self._atomic_number = None
            self._mass = None

    @multimethod
    def __init__(self, name: str, id: int | str) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        Sets the symbol, atomic number and mass of the atom_type from the
        atomic number.

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
        self._symbol, self._atomic_number, self._mass = guess_element(id)

    @multimethod
    def __init__(self, id: int) -> None:
        """
        Constructs all the necessary attributes for the Atom object.

        Determines the symbol, atomic number and mass of the atom_type from the
        atomic number. The name of the atom_type is the symbol in lower case.

        Parameters
        ----------
        id : int
            The identifier of the atom_type. It is
            assumed to be the atomic number of the atom_type.
        """
        self._symbol, self._atomic_number, self._mass = guess_element(id)
        self._name = self._symbol

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
        is_equal &= self.atomic_number == other.atomic_number
        is_equal &= self.symbol == other.symbol
        is_equal &= self.mass == other.mass
        is_equal &= self.name == other.name
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
        The symbol of the atom_type (e.g. 'c')

        Returns
        -------
        str
            The symbol of the atom_type (e.g. 'c')
        """
        return self._symbol

    @property
    def atomic_number(self) -> int | None:
        """
        The atomic number of the atom_type (e.g. 6)

        Returns
        -------
        int
            The atomic number of the atom_type (e.g. 6)
        """
        return self._atomic_number

    @property
    def mass(self) -> Real | None:
        """
        The mass of the atom_type (e.g. 12.0107)

        Returns
        -------
        Real
            The mass of the atom_type (e.g. 12.0107)
        """
        return self._mass


atomicMasses = {"h":    1.00794,    "d":    2.014101778,   "t":    3.0160492675,
                "he":   4.002602,   "li":   6.941,         "be":   9.012182,
                "b":   10.811,      "c":   12.0107,        "n":   14.0067,
                "o":   15.9994,     "f":   18.9984032,     "ne":  20.1797,
                "na":  22.989770,   "mg":  24.3050,        "al":  26.981538,
                "si":  28.0855,     "p":   30.973761,      "s":   32.065,
                "cl":  35.453,      "ar":  39.948,         "k":   39.0983,
                "ca":  40.078,      "sc":  44.955910,      "ti":  47.880,
                "v":   50.9415,     "cr":  51.9961,        "mn":  54.938049,
                "fe":  55.845,      "co":  58.933200,      "ni":  58.6934,
                "cu":  63.546,      "zn":  65.399,         "ga":  69.723,
                "ge":  72.64,       "as":  74.92160,       "se":  78.96,
                "br":  79.904,      "kr":  83.798,         "rb":  85.4678,
                "sr":  87.62,       "y":   88.90585,       "zr":  91.224,
                "nb":  92.90638,    "mo":  95.94,          "tc":  98.9063,
                "ru": 101.07,       "rh": 102.9055,        "pd": 106.42,
                "ag": 107.8682,     "cd": 112.411,         "in": 114.818,
                "sn": 118.71,       "sb": 121.76,          "te": 127.6,
                "i":  126.90447,    "xe": 131.293,         "cs": 132.90546,
                "ba": 137.327,      "la": 138.9055,        "ce": 140.116,
                "pr": 140.90765,    "nd": 144.24,          "pm": 146.9151,
                # TODO: place lr to end of list
                "sm": 150.36,       "lr": 260.1053,        "eu": 151.964,
                "gd": 157.25,       "tb": 158.92534,       "dy": 162.5,
                "ho": 164.93032,    "er": 167.259,         "tm": 168.93421,
                "yb": 173.04,       "lu": 174.967,         "hf": 178.49,
                "ta": 180.9479,     "w":  183.84,          "re": 186.207,
                "os": 190.23,       "ir": 192.217,         "pt": 195.078,
                "au": 196.96655,    "hg": 200.59,          "tl": 204.3833,
                "pb": 207.2,        "bi": 208.98038,       "po": 208.9824,
                "at": 209.9871,     "rn": 222.0176,        "fr": 223.0197,
                "ra": 226.0254,     "ac": 227.0278,        "th": 232.0381,
                "pa": 231.03588,    "u":  238.0289,        "np": 237.0482,
                "pu": 244.0642,     "am": 243.0614,        "cm": 247.0703,
                "bk": 247.0703,     "cf": 251.0796,        "es": 252.0829,
                "fm": 257.0951,     "md": 258.0986,        "no": 259.1009,
                "q":  999.00000,    "x":  999.00000,       "cav": 1000.00000, "sup": 1000000.0, "dum": 1.0}

atomicNumbers = {"h":     1,  "d":     1,  "t":    1,
                 "he":    2,  "li":    3,  "be":   4,
                 "b":     5,  "c":     6,  "n":    7,
                 "o":     8,  "f":     9,  "ne":  10,
                 "na":   11,  "mg":   12,  "al":  13,
                 "si":   14,  "p":    15,  "s":   16,
                 "cl":   17,  "ar":   18,  "k":   19,
                 "ca":   20,  "sc":   21,  "ti":  22,
                 "v":    23,  "cr":   24,  "mn":  25,
                 "fe":   26,  "co":   27,  "ni":  28,
                 "cu":   29,  "zn":   30,  "ga":  31,
                 "ge":   32,  "as":   33,  "se":  34,
                 "br":   35,  "kr":   36,  "rb":  37,
                 "sr":   38,  "y":    39,  "zr":  40,
                 "nb":   41,  "mo":   42,  "tc":  43,
                 "ru":   44,  "rh":   45,  "pd":  46,
                 "ag":   47,  "cd":   48,  "in":  49,
                 "sn":   50,  "sb":   51,  "te":  52,
                 "i":    53,  "xe":   54,  "cs":  55,
                 "ba":   56,  "la":   57,  "ce":  58,
                 "pr":   59,  "nd":   60,  "pm":  61,
                 "sm":   62,  "lr":   103, "eu":  63,  # TODO: place lr to end of list
                 "gd":   64,  "tb":   65,  "dy":  66,
                 "ho":   67,  "er":   68,  "tm":  69,
                 "yb":   70,  "lu":   71,  "hf":  72,
                 "ta":   73,  "w":    74,  "re":  75,
                 "os":   76,  "ir":   77,  "pt":  78,
                 "au":   79,  "hg":   80,  "tl":  81,
                 "pb":   82,  "bi":   83,  "po":  84,
                 "at":   85,  "rn":   86,  "fr":  87,
                 "ra":   88,  "ac":   89,  "th":  90,
                 "pa":   91,  "u":    92,  "np":  93,
                 "pu":   94,  "am":   95,  "cm":  96,
                 "bk":   97,  "cf":   98,  "es":  99,
                 "fm":  100,  "md":  101,  "no": 102,
                 "q":   999,   "x":  999,   "cav": 1000, "sup": 1000000, "dum": 1}
