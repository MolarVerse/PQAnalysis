"""
A module containing tests for the methods of the _PositionsMixin class.
"""

import numpy as np
import pytest

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system.exceptions import AtomicSystemPositionsError
from PQAnalysis.core import Atom, Cell
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.types import PositiveInt, Np1DNumberArray
from PQAnalysis.topology import SelectionCompatible
from PQAnalysis.exceptions import PQTypeError

from . import pytestmark  # pylint: disable=unused-import
from ..conftest import assert_logging_with_exception



class TestPositionsMixin:

    def test_nearest_neighbours_type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "n",
            -1,
            PositiveInt,
            ),
            exception=PQTypeError,
            function=AtomicSystem().nearest_neighbours,
            n=-1
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "selection",
            AtomicSystem(),
            SelectionCompatible,
            ),
            exception=PQTypeError,
            function=AtomicSystem().nearest_neighbours,
            selection=AtomicSystem()
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "use_full_atom_info",
            "2",
            bool,
            ),
            exception=PQTypeError,
            function=AtomicSystem().nearest_neighbours,
            use_full_atom_info="2"
        )

    def test_nearest_neighbours(self, caplog):
        positions = np.array([[0, 0, 0], [10, 0, 0], [2, 1, 0], [1, 0, 0]])
        atoms = [Atom('C'), Atom('H1', 1), Atom(1), Atom('H1', 1)]

        system = AtomicSystem(pos=positions, atoms=atoms)

        indices, distances = system.nearest_neighbours()
        assert np.allclose(indices, [[3], [2], [3], [0]])
        assert np.allclose(
            distances,
            [[1.0],
            [np.sqrt(8 * 8 + 1 * 1)],
            [np.sqrt(2)],
            [1.0]]
        )

        indices, distances = system._nearest_neighbours()
        assert np.allclose(indices, [[3], [2], [3], [0]])
        assert np.allclose(
            distances,
            [[1.0],
            [np.sqrt(8 * 8 + 1 * 1)],
            [np.sqrt(2)],
            [1.0]]
        )

        indices, distances = system.nearest_neighbours(n=2)
        assert np.allclose(indices, [[3, 2], [2, 3], [3, 0], [0, 2]])
        assert np.allclose(
            distances,
            [
            [1.0,
            np.sqrt(5)],
            [np.sqrt(8 * 8 + 1),
            np.sqrt(9 * 9)],
            [np.sqrt(2),
            np.sqrt(5)],
            [1.0,
            np.sqrt(2)]
            ]
        )

        indices, distances = system.nearest_neighbours(
            n=1,
            selection=[Atom('H1', 1)],
            use_full_atom_info=True
        )
        assert np.allclose(indices, [[2], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [1.0]])

        indices, distances = system.nearest_neighbours(n=1, selection=['H1'])
        assert np.allclose(indices, [[2], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1,
            selection=[Atom('H1', 1)],
            use_full_atom_info=False
        )
        assert np.allclose(indices, [[2], [3], [0]])
        assert np.allclose(distances, [[np.sqrt(65)], [np.sqrt(2.0)], [1.0]])

        indices, distances = system.nearest_neighbours(
            n=1,
            selection=np.array([0, 2])
        )
        assert np.allclose(indices, [[3], [3]])
        assert np.allclose(distances, [[1.0], [np.sqrt(2)]])

        system = AtomicSystem(atoms=atoms)
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=AtomicSystem.__qualname__,
            logging_level="ERROR",
            message_to_test=AtomicSystemPositionsError.message,
            exception=AtomicSystemPositionsError,
            function=system.nearest_neighbours,
            n=1
        )

    def test_image(self):
        system = AtomicSystem(
            atoms=[Atom('C'),
            Atom('H1',
            1)],
            pos=np.array([[0,
            0,
            0],
            [10,
            0,
            0]]),
            cell=Cell(8,
            8,
            8)
        )
        system.image()

        assert np.allclose(system.pos, np.array([[0, 0, 0], [2, 0, 0]]))

    def test_center_type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "position",
            "2",
            Np1DNumberArray,
            ),
            exception=PQTypeError,
            function=AtomicSystem().center,
            position="2"
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "image",
            "2",
            bool,
            ),
            exception=PQTypeError,
            function=AtomicSystem().center,
            position=np.array([1,
            0,
            0]),
            image="2"
        )

    def test_center(self):
        system = AtomicSystem(
            atoms=[Atom('C'),
            Atom('H1',
            1)],
            pos=np.array([[0,
            0,
            0],
            [10,
            0,
            0]]),
            cell=Cell(8,
            8,
            8)
        )
        system.center(np.array([1, 0, 0]), image=False)

        assert np.allclose(system.pos, np.array([[-1, 0, 0], [9, 0, 0]]))

        system.center(np.array([1, 0, 0]), image=True)

        assert np.allclose(system.pos, np.array([[-2, 0, 0], [0, 0, 0]]))
