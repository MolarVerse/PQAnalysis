import pytest
import numpy as np

from PQAnalysis.traj.frame import Frame
from PQAnalysis.selection.selection import Selection


def test__init__():

    frame = Frame(2, [[0, 0, 0], [1, 0, 0]], ['C', 'H'])

    with pytest.raises(ValueError) as exception:
        Selection(['C', 'H'])
    assert str(
        exception.value) == "Frame must be provided when selection is a string."

    selection = Selection(['C', 'H'], frame)
    assert np.allclose(selection.selection, [0, 1])

    selection = Selection([0, 1], frame)
    assert np.allclose(selection.selection, [0, 1])

    selection = Selection(0, frame)
    assert np.allclose(selection.selection, [0])

    with pytest.raises(ValueError) as exception:
        Selection('C')
    assert str(
        exception.value) == "Frame must be provided when selection is a string."

    selection = Selection('C', frame)
    assert np.allclose(selection.selection, [0])

    with pytest.raises(TypeError) as exception:
        Selection(1.2)
    assert str(exception.value) == "Invalid selection type."
