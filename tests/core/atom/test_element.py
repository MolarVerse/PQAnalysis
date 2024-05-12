import pytest
import numpy as np

from PQAnalysis.core import Element, ElementNotFoundError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark
from ...conftest import assert_logging_with_exception



class TestElement:

    def test__init__(self, caplog):
        assert_logging_with_exception(
            caplog,
            Element.__qualname__,
            "ERROR",
            "Id C1 is not a valid element identifier.",
            ElementNotFoundError,
            Element,
            'C1'
        )

        assert_logging_with_exception(
            caplog,
            "TypeChecking",
            "ERROR",
            get_type_error_message("element_id",
            1.2,
            int | str | None),
            PQTypeError,
            Element,
            1.2
        )

        element = Element()
        assert element.symbol == None
        assert element.atomic_number == None
        assert element.mass == None

        element = Element('C')
        assert element.symbol == 'c'
        assert element.atomic_number == 6
        assert np.isclose(element.mass, 12.0107)

        element = Element(2)
        assert element.symbol == 'he'
        assert element.atomic_number == 2
        assert np.isclose(element.mass, 4.002602)

    def test__eq__(self):
        element1 = Element('C')
        element2 = Element(6)
        element3 = Element('H')
        assert element1 == element2
        assert element1 != element3
        assert element2 != element3

        assert element1 != 1

    def test__str__(self):
        element = Element('C')
        assert str(element) == 'Element(c, 6, 12.0107)'

    def test__repr__(self):
        element = Element('C')
        assert repr(element) == 'Element(c, 6, 12.0107)'
