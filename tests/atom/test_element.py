import pytest
import numpy as np

from PQAnalysis.atom.element import Element


def test__init__():
    element = Element('C')
    assert element.name == 'c'
    assert element.atomic_number == 6
    assert np.isclose(element.mass, 12.0107)

    element = Element(2)
    assert element.name == 'he'
    assert element.atomic_number == 2
    assert np.isclose(element.mass, 4.002602)

    with pytest.raises(ValueError) as exception:
        Element(1.2)
    assert str(
        exception.value) == "Invalid element id - must be either atomic number or element name."


def test__eq__():
    element1 = Element('C')
    element2 = Element(6)
    element3 = Element('H')
    assert element1 == element2
    assert element1 != element3
