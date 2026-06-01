import numpy as np
import pytest

from PQAnalysis.io.traj_file import _process_lines_py

from . import pytestmark



def test_process_lines_with_atoms():
    atoms, xyz = _process_lines_py.process_lines_with_atoms(
        ["h 1.0 2.0 3.0", "o 4.0 5.0 6.0"], 2
    )

    assert atoms == ["h", "o"]
    assert np.allclose(xyz, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


def test_process_lines():
    xyz = _process_lines_py.process_lines(
        ["h 1.0 2.0 3.0", "o 4.0 5.0 6.0"], 2
    )

    assert np.allclose(xyz, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


@pytest.mark.parametrize("line", ["h 1.0 2.0", "h 1.0 2.0 bad"])
def test_process_lines_invalid_line(line):
    with pytest.raises(ValueError) as exception:
        _process_lines_py.process_lines_with_atoms([line], 1)

    assert str(exception.value) == "Could not parse line"
