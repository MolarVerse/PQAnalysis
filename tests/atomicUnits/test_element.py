import pytest
import numpy as np

from PQAnalysis.atomicUnits.element import Element, Elements
from PQAnalysis.utils.exceptions import ElementNotFoundError
from PQAnalysis.traj.selection import Selection


class TestElement:
    def test__init__(self):
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

    def test__eq__(self):
        element1 = Element('C')
        element2 = Element(6)
        element3 = Element('H')
        assert element1 == element2
        assert element1 != element3

        assert element1 != 1


class TestElements():

    # this test is representative of the elements setter method
    def test__init__(self):
        with pytest.raises(TypeError) as exception:
            Elements([1, 2.0, 3.1])
        assert str(
            "elements must be either a list of strings (elements symbols), ints (atomic numbers) or a list of Element objects.")

        with pytest.raises(ElementNotFoundError) as exception:
            Elements(["notAnElementSymbol"])
        assert str(
            "Element with id notanelementsymbol is not a valid element identifier.")

        with pytest.raises(ElementNotFoundError) as exception:
            Elements([-1])
        assert str("Element with id -1 is not a valid element identifier.")

        assert Elements(Elements(['C', 'H', 'H', 'O'])
                        ) == Elements(['C', 'H', 'H', 'O'])

        assert Elements(['C', 'H', 'H', 'O']) == Elements(['C', 'H', 'H', 'O'])

        assert Elements([Element('C'), Element('H'), Element('H'), Element('O')]
                        ) == Elements(['C', 'H', 'H', 'O'])

        assert Elements(['C', 'H', 'H', 'O']) == Elements(
            [Element('C'), Element('H'), Element('H'), Element('O')])

        assert Elements([6, 1, 1, 8]) == Elements(['C', 'H', 'H', 'O'])

    def test_names(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        assert np.array_equal(elements.names, ['c', 'h', 'h', 'o'])

    def test_atomic_numbers(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        assert np.array_equal(elements.atomic_numbers, [6, 1, 1, 8])

    def test_masses(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        assert np.allclose(elements.masses, np.array(
            [12.0107, 1.00794, 1.00794, 15.9994]))

    def test__eq__(self):
        elements1 = Elements(['C', 'H', 'H', 'O'])
        elements2 = Elements([6, 1, 1, 8])
        elements3 = Elements(['C', 'H', 'H', 'O', 'O'])
        assert elements1 == elements2
        assert elements1 != elements3
        assert elements1 != 1

    def test__getitem__(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        assert elements[0] == Element('C')
        assert elements[1] == Element('H')
        assert elements[2] == Element('H')
        assert elements[3] == Element('O')

        with pytest.raises(IndexError) as exception:
            elements[4]
        assert str(
            exception.value) == "index 4 is out of bounds for axis 0 with size 4"

        assert elements[Selection([0, 1, 2])] == Elements(
            ['C', 'H', 'H'])

    def test__len__(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        assert len(elements) == elements.n_atoms == 4

    def test__iter__(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        for i, element in enumerate(elements):
            assert element == elements[i]

    def test__contains__(self):
        elements = Elements(['C', 'H', 'H', 'O'])
        assert Element('C') in elements
        assert Element('H') in elements
        assert Element('O') in elements
        assert Element('N') not in elements
