import pytest
import numpy as np

from PQAnalysis.core import Element
from PQAnalysis.core import ElementNotFoundError


class TestElement:
    def test__init__(self):
        with pytest.raises(ElementNotFoundError) as exception:
            Element('C1')
        assert str(exception.value) == "Id C1 is not a valid element identifier."

        with pytest.raises(ElementNotFoundError) as exception:
            Element(-1)
        assert str(exception.value) == "Id -1 is not a valid element identifier."

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
