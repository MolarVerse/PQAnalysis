import pytest
import numpy as np

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system.exceptions import AtomicSystemError
from PQAnalysis.core import Atom
from PQAnalysis.type_checking import get_type_error_message

from . import pytestmark  # pylint: disable=unused-import
from ..conftest import assert_logging_with_exception



def test_check_atom_number_setter(caplog):
    system = AtomicSystem()

    assert_logging_with_exception(
        caplog,
        AtomicSystem.__qualname__,
        "ERROR",
        (
        "The number of atoms in the AtomicSystem object "
        "have to be equal to the number of atoms "
        "in the new array in order to set the property."
        ),
        AtomicSystemError,
        AtomicSystem.pos.__set__,
        system,
        np.array([[0,
        0,
        0]])
    )

    system = AtomicSystem(atoms=[Atom('C'), Atom('H')])
    system.pos = np.array([[0, 0, 0], [1, 1, 1]])

    assert np.allclose(system.pos, np.array([[0, 0, 0], [1, 1, 1]]))
