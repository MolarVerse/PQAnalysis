import numpy as np
import pytest

from . import pytestmark

from PQAnalysis.core import Cell



class TestCell:

    def test__init__(self):
        cell = Cell(1, 2, 3)
        assert cell.x == 1
        assert cell.y == 2
        assert cell.z == 3
        assert cell.alpha == 90
        assert cell.beta == 90
        assert cell.gamma == 90
        assert np.allclose(
            cell.box_matrix, np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        )

        cell = Cell(1, 2, 3, 60, 90, 120)
        assert cell.x == 1
        assert cell.y == 2
        assert cell.z == 3
        assert cell.alpha == 60
        assert cell.beta == 90
        assert cell.gamma == 120
        assert np.allclose(
            cell.box_matrix,
            np.array(
                [[1, -1, 0], [0, 1.73205081, 1.73205081], [0, 0, 2.44948974]]
            )
        )

    def test_box_lengths(self):
        cell = Cell(1, 2, 3)
        assert np.allclose(cell.box_lengths, np.array([1, 2, 3]))

    def test_set_box_lengths(self):
        cell = Cell(1, 2, 3)
        cell.box_lengths = np.array([2, 3, 4])
        assert np.allclose(cell.box_lengths, np.array([2, 3, 4]))
        assert np.allclose(
            cell.box_matrix, np.array([[2, 0, 0], [0, 3, 0], [0, 0, 4]])
        )

    def test_box_angles(self):
        cell = Cell(1, 2, 3)
        assert np.allclose(cell.box_angles, np.array([90, 90, 90]))

    def test_set_box_angles(self):
        cell = Cell(1, 2, 3)
        cell.box_angles = np.array([60, 90, 120])
        assert np.allclose(cell.box_angles, np.array([60, 90, 120]))
        assert np.allclose(
            cell.box_matrix,
            np.array(
                [[1, -1, 0], [0, 1.73205081, 1.73205081], [0, 0, 2.44948974]]
            )
        )

    def test_volume(self):
        cell = Cell(1, 2, 3, 60, 90, 120)
        box_angles = np.deg2rad(cell.box_angles)
        assert np.isclose(
            cell.volume,
            np.prod(cell.box_lengths) * np.sqrt(
                1 - sum(np.cos(box_angles)**2) +
                2 * np.prod(np.cos(box_angles))
            )
        )

    def test_bounding_edges(self):
        cell = Cell(1, 2, 3, 60, 90, 120)
        edges = cell.bounding_edges
        assert np.allclose(edges[0], np.array([0, -1.73205081, -1.22474487]))
        assert np.allclose(edges[1], np.array([0, 0, 1.22474487]))
        assert np.allclose(edges[2], np.array([-1, 0, -1.22474487]))
        assert np.allclose(edges[3], np.array([-1, 1.73205081, 1.22474487]))
        assert np.allclose(edges[4], np.array([1, -1.73205081, -1.22474487]))
        assert np.allclose(edges[5], np.array([1, 0, 1.22474487]))
        assert np.allclose(edges[6], np.array([0, 0, -1.22474487]))
        assert np.allclose(edges[7], np.array([0, 1.73205081, 1.22474487]))

    def test__eq__(self):

        cell1 = Cell(1, 2, 3, 60, 90, 120)
        cell2 = Cell(1, 2, 3, 60, 90, 120)
        assert cell1 == cell2

        cell1 = Cell(1, 2, 3, 60, 90, 120)
        cell2 = Cell(1, 2, 3, 60, 90, 90)
        assert cell1 != cell2

        cell1 = Cell(1, 2, 3, 60, 90, 120)
        assert cell1 != 1

        cell1 = Cell(1, 2, 3, 60, 90, 120)
        cell2 = Cell(1, 3, 3, 60, 90, 120)

        assert cell1 != cell2

    def test_isclose(self):
        cell1 = Cell(1, 2, 3, 60, 90, 120.1)
        cell2 = Cell(1, 2, 3, 60, 90, 120.0)
        assert cell1.isclose(cell2, atol=1e-1)
        assert not cell1.isclose(cell2, atol=1e-2)
        assert cell1.isclose(cell2, rtol=1e-3)
        assert not cell1.isclose(cell2, rtol=1e-4)

    def test_image(self):
        cell = Cell(1, 2, 3, 60, 90, 120)
        assert np.allclose(
            cell.image(np.array([0, 0, 0])), np.array([0, 0, 0])
        )
        assert np.allclose(
            cell.image(np.array([0.75, 0.5, 0.5])),
            np.array([-0.25, 0.5, 0.5])
        )
        assert np.allclose(
            cell.image(np.array([1, 2, 3])),
            np.array([0., 0.267949192, 0.550510257])
        )
        assert np.allclose(
            cell.image(np.array([-1, -2, -3])),
            np.array([0., -0.267949192, -0.550510257])
        )

    def test__str__(self):
        cell = Cell(1, 2, 3, 60, 90, 120)
        assert str(
            cell
        ) == "Cell(x=1.0, y=2.0, z=3.0, alpha=60.0, beta=90.0, gamma=120.0)"
        assert str(cell) == repr(cell)

        cell = Cell()
        assert str(cell) == "Cell()"
        assert str(cell) == repr(cell)

    def test_init_from_box_matrix(self):
        box_matrix = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        cell = Cell.init_from_box_matrix(box_matrix)
        assert cell.x == 1
        assert cell.y == 2
        assert cell.z == 3
        assert cell.alpha == 90
        assert cell.beta == 90
        assert cell.gamma == 90
        assert np.allclose(cell.box_matrix, box_matrix)
        assert np.allclose(cell.box_lengths, np.array([1, 2, 3]))
        assert np.allclose(cell.box_angles, np.array([90, 90, 90]))
        assert np.isclose(cell.volume, 6)

        box_matrix = np.array(
            [[1, -1, 0], [0, 1.73205081, 1.73205081], [0, 0, 2.44948974]]
        )
        cell = Cell.init_from_box_matrix(box_matrix)
        assert np.isclose(cell.x, 1)
        assert np.isclose(cell.y, 2)
        assert np.isclose(cell.z, 3)
        assert np.isclose(cell.alpha, 60)
        assert np.isclose(cell.beta, 90)
        assert np.isclose(cell.gamma, 120)
