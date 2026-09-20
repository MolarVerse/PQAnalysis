import numpy as np
import pytest

from PQAnalysis.io.traj_file import _process_lines_py
from PQAnalysis.io.traj_file.exceptions import FrameReaderError

from . import pytestmark

try:
    from PQAnalysis.io.traj_file import process_lines as _process_lines_ext
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _process_lines_ext = None

#: Both line parser implementations.
PARSER_MODULES = [
    pytest.param(_process_lines_py, id="python-fallback"),
    pytest.param(
        _process_lines_ext,
        id="cython",
        marks=pytest.mark.skipif(
            _process_lines_ext is None,
            reason="compiled line parser not available",
        ),
    ),
]



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


@pytest.mark.parametrize("module", PARSER_MODULES)
@pytest.mark.parametrize("length", [255, 300])
def test_process_lines_with_atoms_overlong_label(module, length):
    line = "A" * length + " 1.0 2.0 3.0"

    with pytest.raises(FrameReaderError) as exception:
        module.process_lines_with_atoms([line], 1)

    assert str(exception.value) == (
        "Atom type is too long: "
        "the maximum supported length is 254 characters"
    )


@pytest.mark.parametrize("module", PARSER_MODULES)
def test_process_lines_with_atoms_long_custom_label(module):
    atoms, xyz = module.process_lines_with_atoms(
        ["custom_label1 1.0 2.0 3.0"], 1
    )

    assert atoms == ["custom_label1"]
    assert np.allclose(xyz, [[1.0, 2.0, 3.0]])
