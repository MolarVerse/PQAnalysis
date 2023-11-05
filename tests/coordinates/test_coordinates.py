import numpy as np
import pytest

from PQAnalysis.coordinates.coordinates import Coordinates, image
from PQAnalysis.core.cell import Cell


def test_image():
    coords = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
    coords = image(coords)
    assert np.allclose(coords.xyz, np.array([[1, 2, 3], [4, 5, 6]]))
    assert coords.cell is None

    coords = Coordinates(
        np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.5, 1.5, 1.5))
    coords = image(coords)
    assert np.allclose(coords.xyz, np.array([[-0.5, 0.5, 0], [-0.5, 0.5, 0]]))
    assert coords.cell == Cell(1.5, 1.5, 1.5)


class TestCoordinates:
    def test__init__(self):

        coords = Coordinates()
        assert np.array_equal(coords.xyz, np.zeros((0, 3)))
        assert coords.cell is None

        coords = Coordinates(np.array([1, 2, 3]))
        assert np.allclose(coords.xyz, np.array([[1, 2, 3]]))
        assert coords.cell is None

        coords = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        assert np.allclose(coords.xyz, np.array([[1, 2, 3], [4, 5, 6]]))
        assert coords.cell is None

        coords = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        assert np.allclose(coords.xyz, np.array([[1, 2, 3], [4, 5, 6]]))
        assert coords.cell == Cell(1.0, 1.0, 1.0)

    def test_n_atoms(self):
        coords = Coordinates()
        assert coords.n_atoms == 0

        coords = Coordinates(np.array([1, 2, 3]))
        assert coords.n_atoms == 1

        coords = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        assert coords.n_atoms == 2

    def test_PBC(self):
        coords = Coordinates()
        assert coords.PBC is False

        coords = Coordinates(np.array([1, 2, 3]))
        assert coords.PBC is False

        coords = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        assert coords.PBC is False

        coords = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        assert coords.PBC is True

    def test__add__(self):
        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords3 = coords1 + coords2
        assert np.allclose(coords3.xyz, np.array([[2, 4, 6], [8, 10, 12]]))
        assert coords3.cell is None

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords3 = coords1 + coords2
        assert np.allclose(coords3.xyz, np.array([[2, 4, 6], [8, 10, 12]]))
        assert coords3.cell == Cell(1.0, 1.0, 1.0)

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords3 = coords1 + coords2
        assert np.allclose(coords3.xyz, np.array([[2, 4, 6], [8, 10, 12]]))
        assert coords3.cell == Cell(1.0, 1.0, 1.0)

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords3 = coords1 + coords2
        assert np.allclose(coords3.xyz, np.array([[2, 4, 6], [8, 10, 12]]))
        assert coords3.cell == Cell(1.0, 1.0, 1.0)

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(2.0, 2.0, 2.0))
        with pytest.raises(ValueError) as exception:
            coords1 + coords2
        assert str(
            exception.value) == 'The cells of the two Coordinates objects must be the same.'

    def test__sub__(self):
        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords3 = coords1 - coords2
        assert np.allclose(coords3.xyz, np.array([[0, 0, 0], [0, 0, 0]]))
        assert coords3.cell is None

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords3 = coords1 - coords2
        assert np.allclose(coords3.xyz, np.array([[0, 0, 0], [0, 0, 0]]))
        assert coords3.cell == Cell(1.0, 1.0, 1.0)

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords3 = coords1 - coords2
        assert np.allclose(coords3.xyz, np.array([[[0, 0, 0], [0, 0, 0]]]))
        assert coords3.cell == Cell(1.0, 1.0, 1.0)

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords3 = coords1 - coords2
        assert np.allclose(coords3.xyz, np.array([[[0, 0, 0], [0, 0, 0]]]))
        assert coords3.cell == Cell(1.0, 1.0, 1.0)

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(2.0, 2.0, 2.0))
        with pytest.raises(ValueError) as exception:
            coords1 - coords2
        assert str(
            exception.value) == 'The cells of the two Coordinates objects must be the same.'

    def test__eq__(self):
        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        assert coords1 == coords2

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 7]]))
        assert coords1 != coords2

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        assert coords1 != coords2

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        assert coords1 == coords2

    def test__getitem__(self):
        coords = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        assert np.allclose(coords[0].xyz, np.array([1, 2, 3]))
        assert np.allclose(coords[1].xyz, np.array([4, 5, 6]))
        assert np.allclose(coords[-1].xyz, np.array([4, 5, 6]))
        assert np.allclose(coords[-2].xyz, np.array([1, 2, 3]))
        with pytest.raises(IndexError) as exception:
            coords[2]
        assert str(
            exception.value) == 'index 2 is out of bounds for axis 0 with size 2'

    def test__cell_compatible__(self):
        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        assert coords1.__cell_compatible__(coords2) is True

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        assert coords1.__cell_compatible__(coords2) is True

        coords1 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(1.0, 1.0, 1.0))
        coords2 = Coordinates(
            np.array([[1, 2, 3], [4, 5, 6]]), Cell(2.0, 2.0, 2.0))
        assert coords1.__cell_compatible__(coords2) is False

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]), Cell(1, 1, 1))
        coords2 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]), Cell(1, 1, 1))
        assert coords1.__cell_compatible__(coords2) is True

    def test_append(self):
        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(np.array([[7, 8, 9], [10, 11, 12]]))
        coords1.append(coords2)
        assert np.allclose(coords1.xyz, np.array(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]))
        assert coords1.cell is None

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]), Cell(1, 1, 1))
        coords2 = Coordinates(np.array([[7, 8, 9], [10, 11, 12]]))
        coords1.append(coords2)
        assert np.allclose(coords1.xyz, np.array(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]))
        assert coords1.cell == Cell(1, 1, 1)

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]))
        coords2 = Coordinates(
            np.array([[7, 8, 9], [10, 11, 12]]), Cell(1, 1, 1))
        with pytest.raises(ValueError) as exception:
            coords1.append(coords2)
        assert str(
            exception.value) == 'The cells of the two Coordinates objects must be the same or the one of the appended None.'

        coords1 = Coordinates(np.array([[1, 2, 3], [4, 5, 6]]), Cell(1, 1, 1))
        coords2 = Coordinates(
            np.array([[7, 8, 9], [10, 11, 12]]), Cell(1, 1, 1))
        coords1.append(coords2)
        assert np.allclose(coords1.xyz, np.array(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]))
        assert coords1.cell == Cell(1, 1, 1)
