import pytest

from . import pytestmark

from PQAnalysis.core import Cell
from PQAnalysis.traj import check_trajectory_pbc, check_trajectory_vacuum



def test_check_trajectory_PBC():
    cells = [Cell(10, 10, 10), Cell(10, 10, 10)]

    assert check_trajectory_pbc(cells) == True

    cells = [Cell(10, 10, 10), Cell(10, 10, 10), Cell()]

    assert check_trajectory_pbc(cells) == False

    cells = []

    assert check_trajectory_pbc(cells) == False



def test_check_trajectory_vacuum():
    cells = [Cell(), Cell(), Cell()]

    assert check_trajectory_vacuum(cells) == True

    cells = [Cell(), Cell(), Cell(10, 10, 10)]

    assert check_trajectory_vacuum(cells) == False
