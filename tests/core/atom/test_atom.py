import pytest
import numpy as np

from multimethod import DispatchError

from PQAnalysis.core.atom import Atom, Element
from PQAnalysis.core.exceptions import ElementNotFoundError, AtomError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark
from ...conftest import assert_logging_with_exception



class TestAtom:

    def test__init__(self, caplog):
        element = Atom('C')
        assert element.symbol == 'c'
        assert element.atomic_number == 6
        assert np.isclose(element.mass, 12.0107)

        element = Atom(2)
        assert element.symbol == 'he'
        assert element.atomic_number == 2
        assert np.isclose(element.mass, 4.002602)

        element = Atom('C', 6)
        assert element.name == 'C'
        assert element.symbol == 'c'
        assert element.atomic_number == 6
        assert np.isclose(element.mass, 12.0107)

        element = Atom('C1', 'C')
        assert element.name == 'C1'
        assert element.symbol == 'c'
        assert element.atomic_number == 6
        assert np.isclose(element.mass, 12.0107)

        element = Atom('C1', use_guess_element=False)
        assert element.name == 'C1'
        assert element.symbol is None
        assert element.atomic_number is None
        assert element.mass is None

        assert_logging_with_exception(
            caplog,
            Element.__qualname__,
            "ERROR",
            "Id C1 is not a valid element identifier.",
            ElementNotFoundError,
            Atom,
            'C1'
        )

        assert_logging_with_exception(
            caplog,
            "TypeChecking",
            "ERROR",
            get_type_error_message("name", 1.2, str | int | Element),
            PQTypeError,
            Atom,
            1.2
        )

        assert_logging_with_exception(
            caplog,
            Element.__qualname__,
            "ERROR",
            "Id -1 is not a valid element identifier.",
            ElementNotFoundError,
            Atom,
            -1
        )

        assert_logging_with_exception(
            caplog,
            Atom.__qualname__,
            "ERROR",
            "The name of the atom_type cannot be an integer if the id is given.",
            AtomError,
            Atom,
            1,
            2
        )

    def test__eq__(self):
        element1 = Atom('C')
        element2 = Atom('C', 6)
        element3 = Atom('H')
        assert element1 == element2
        assert element1 != element3

        assert element1 != 1

    def test__str__(self):
        element = Atom('C')
        assert str(element) == 'Atom(C, 6, c, 12.0107)'

    def test__repr__(self):
        element = Atom('C')
        assert repr(element) == str(element)

    def test_property_element(self):
        element = Atom('C')
        assert element.element == Element('C')
