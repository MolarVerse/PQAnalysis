import numpy as np
import pytest

from PQAnalysis.traj.frame import Frame
from PQAnalysis.pbc.cell import Cell
from PQAnalysis.selection.selection import Selection


def test__init__():
    frame = Frame(2, [[0, 0, 0], [1, 0, 0]], ['C', 'H'])

    assert frame.n_atoms == 2
    assert np.allclose(frame.xyz, np.array([[0, 0, 0], [1, 0, 0]]))
    assert frame.atoms[0] == 'C'
    assert frame.atoms[1] == 'H'
    assert frame.cell is None

    frame = Frame(2, [[0, 0, 0], [1, 0, 0]], ['C', 'H'], Cell(10, 10, 10))

    assert frame.n_atoms == 2
    assert np.allclose(frame.xyz, np.array([[0, 0, 0], [1, 0, 0]]))
    assert frame.atoms[0] == 'C'
    assert frame.atoms[1] == 'H'
    assert frame.cell == Cell(10, 10, 10)


def test_PBC():
    frame = Frame(2, [[0, 0, 0], [1, 0, 0]], ['C', 'H'], Cell(10, 10, 10))
    assert frame.PBC

    frame = Frame(2, [[0, 0, 0], [1, 0, 0]], ['C', 'H'])
    assert not frame.PBC


def test__getindex__():
    frame = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                  ['C', 'H', 'H'], Cell(10, 10, 10))

    assert frame[0].n_atoms == 1
    assert np.allclose(frame[0].xyz, np.array([[0, 0, 0]]))
    assert frame[0].atoms == 'C'
    assert frame[0].cell == Cell(10, 10, 10)

    sel = Selection(['H'], frame)
    assert frame[sel].n_atoms == 2
    assert np.allclose(frame[sel].xyz, np.array([[1, 0, 0], [2, 0, 0]]))
    assert frame[sel].atoms[0] == 'H'
    assert frame[sel].atoms[1] == 'H'
    assert frame[sel].cell == Cell(10, 10, 10)

    sel = Selection([1, 2], frame)
    assert frame[sel].n_atoms == 2
    assert np.allclose(frame[sel].xyz, np.array([[1, 0, 0], [2, 0, 0]]))
    assert frame[sel].atoms[0] == 'H'
    assert frame[sel].atoms[1] == 'H'
    assert frame[sel].cell == Cell(10, 10, 10)

    sel = Selection(['O'], frame)
    with pytest.raises(ValueError) as exception:
        frame[sel]
    assert str(exception.value) == 'Selection is empty.'


def test_compute_com():
    frame = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                  ['C', 'H', 'H'], Cell(10, 10, 10))

    com_frame = frame.compute_com()
    assert com_frame.n_atoms == 1
    assert np.allclose(com_frame.xyz, np.array([[0.21557785, 0, 0]]))
    assert com_frame.atoms == "chh"

    sel = Selection(['H'], frame)
    com_frame = frame[sel].compute_com()
    assert com_frame.n_atoms == 1
    assert np.allclose(com_frame.xyz, np.array([[1.5, 0, 0]]))
    assert com_frame.atoms == "hh"

    with pytest.raises(ValueError) as exception:
        frame.compute_com(group=2)
    assert str(
        exception.value) == 'Number of atoms in selection is not a multiple of group.'


def test__eq__():
    frame1 = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                   ['C', 'H', 'H'], Cell(10, 10, 10))
    frame2 = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                   ['C', 'H', 'H'], Cell(10, 10, 10))
    assert frame1 == frame2
    frame2 = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                   ['C', 'H', 'H'])
    assert frame1 != frame2

    frame2 = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 1]],
                   ['C', 'H', 'H'], Cell(10, 10, 10))
    assert frame1 != frame2

    frame2 = Frame(3, [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                   ['C', 'H', 'O'], Cell(10, 10, 10))
    assert frame1 != frame2

    assert frame1 != 1
