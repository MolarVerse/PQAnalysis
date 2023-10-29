import pytest
import numpy as np

from PQAnalysis.atom.element import Element
from PQAnalysis.utils.exceptions import ElementNotFoundError


def test__init__():
    element = Element('C')
    assert element.name == 'c'
    assert element.atomic_number == 6
    assert np.isclose(element.mass, 12.0107)

    element = Element(2)
    assert element.name == 'he'
    assert element.atomic_number == 2
    assert np.isclose(element.mass, 4.002602)

    with pytest.raises(TypeError) as exception:
        Element(1.2)
    assert str(
        exception.value) == "Invalid type for element id - must be either int (atomic number) or str (element name)."

    with pytest.raises(ElementNotFoundError) as exception:
        Element(-1)
    assert str(
        exception.value) == "Element with id -1 is not a valid element identifier."

    with pytest.raises(ElementNotFoundError) as exception:
        Element("CH")
    assert str(
        exception.value) == "Element with id CH is not a valid element identifier."


def test__eq__():
    element1 = Element('C')
    element2 = Element(6)
    element3 = Element('H')
    assert element1 == element2
    assert element1 != element3

    assert element1 != 1
