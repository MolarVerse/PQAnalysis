import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.core import Atom


def test_check_atom_number_setter():
    system = AtomicSystem()

    with pytest.raises(ValueError) as exception:
        system.pos = np.array([[0, 0, 0]])
    assert str(exception.value) == "The number of atoms in the AtomicSystem object have to be equal to the number of atoms in the new array in order to set the property."

    system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
    system.pos = np.array([[0, 0, 0], [1, 1, 1]])

    assert np.allclose(system.pos, np.array([[0, 0, 0], [1, 1, 1]]))
