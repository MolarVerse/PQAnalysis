import pytest
import numpy as np

from multimethod import DispatchError

from PQAnalysis.core.atom import Atom, guess_element
from PQAnalysis.utils.exceptions import ElementNotFoundError


def test_guess_element():
    with pytest.raises(ElementNotFoundError) as exception:
        guess_element("CH")
    assert str(exception.value) == "Id CH is not a valid element identifier."

    with pytest.raises(ElementNotFoundError) as exception:
        guess_element(-1)
    assert str(exception.value) == "Id -1 is not a valid element identifier."

    symbol, atomic_number, mass = guess_element('C')
    assert symbol == 'c'
    assert atomic_number == 6
    assert np.isclose(mass, 12.0107)

    symbol, atomic_number, mass = guess_element(6)
    assert symbol == 'c'
    assert atomic_number == 6
    assert np.isclose(mass, 12.0107)


class TestAtom:
    def test__init__(self):
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

        with pytest.raises(ElementNotFoundError) as exception:
            Atom('C1')
        assert str(exception.value) == "Id C1 is not a valid element identifier."

        with pytest.raises(DispatchError) as exception:
            Atom(1.2)

        with pytest.raises(ElementNotFoundError) as exception:
            Atom(-1)
        assert str(
            exception.value) == "Id -1 is not a valid element identifier."

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
