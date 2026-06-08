import builtins
import importlib

import numpy as np
import pytest

from PQAnalysis.core import Cell, Atom
from PQAnalysis.io import _FrameReader
from PQAnalysis.io.traj_file import _process_lines_py
import PQAnalysis.io.traj_file.frame_reader as frame_reader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError
from PQAnalysis.topology import Topology
from PQAnalysis.traj import MDEngineFormat, TrajectoryFormat
from PQAnalysis.traj.exceptions import TrajectoryFormatError

from . import pytestmark


def test_frame_reader_base_class_is_abstract():
    with pytest.raises(TypeError):
        frame_reader.BaseFrameReader()


def test_frame_reader_compatibility_alias():
    reader = frame_reader._FrameReader()

    assert isinstance(reader, frame_reader.BaseFrameReader)
    assert isinstance(reader, frame_reader.XYZFrameReader)


@pytest.mark.parametrize(
    "traj_format",
    [
        TrajectoryFormat.XYZ,
        TrajectoryFormat.VEL,
        TrajectoryFormat.FORCE,
        TrajectoryFormat.CHARGE,
    ],
)
def test_get_frame_reader(traj_format):
    reader = frame_reader.get_frame_reader(traj_format)

    assert isinstance(reader, frame_reader.XYZFrameReader)


def test_frame_reader_uses_python_fallback(monkeypatch):
    real_import = builtins.__import__

    def fail_process_lines_import(
        name, globals=None, locals=None, fromlist=(), level=0
    ):
        if level == 1 and name == "process_lines":
            raise ModuleNotFoundError(
                "No module named 'PQAnalysis.io.traj_file.process_lines'"
            )

        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fail_process_lines_import)

    reloaded_frame_reader = importlib.reload(frame_reader)

    assert (
        reloaded_frame_reader.process_lines_with_atoms
        is _process_lines_py.process_lines_with_atoms
    )

    monkeypatch.undo()
    importlib.reload(frame_reader)


def test_get_frame_reader_keeps_md_format():
    reader = frame_reader.get_frame_reader(
        TrajectoryFormat.XYZ, md_format="qmcfc"
    )

    assert isinstance(reader, frame_reader.XYZFrameReader)
    assert reader.md_format is MDEngineFormat.QMCFC


def test_get_frame_reader_rejects_unsupported_format(monkeypatch):
    monkeypatch.setattr(
        frame_reader.BaseFrameReader.logger,
        "error",
        lambda *args, **kwargs: None,
    )

    with pytest.raises(FrameReaderError) as exception:
        frame_reader.get_frame_reader(TrajectoryFormat.AUTO)

    assert str(exception.value) == (
        "Invalid TrajectoryFormat given. "
        "traj_format=<TrajectoryFormat.AUTO: 'auto'>"
    )


@pytest.mark.parametrize(
    ("traj_format", "attribute", "expected"),
    [
        (TrajectoryFormat.XYZ, "pos", [[1.0, 2.0, 3.0]]),
        (TrajectoryFormat.VEL, "vel", [[1.0, 2.0, 3.0]]),
        (TrajectoryFormat.FORCE, "forces", [[1.0, 2.0, 3.0]]),
        (TrajectoryFormat.CHARGE, "charges", [1.0]),
    ],
)
def test_xyz_frame_reader_preserves_qmcfc_handling(
    traj_format, attribute, expected
):
    reader = frame_reader.XYZFrameReader(md_format="qmcfc")
    if traj_format is TrajectoryFormat.CHARGE:
        frame_string = "2\n\nX 0.0\nh 1.0"
    else:
        frame_string = "2\n\nX 0.0 0.0 0.0\nh 1.0 2.0 3.0"

    frame = reader.read(frame_string, traj_format=traj_format)

    assert frame.n_atoms == 1
    assert frame.atoms == [Atom("h")]
    assert np.allclose(getattr(frame, attribute), expected)


def test_xyz_frame_reader_rejects_invalid_qmcfc_first_atom():
    reader = frame_reader.XYZFrameReader(md_format="qmcfc")

    with pytest.raises(FrameReaderError) as exception:
        reader.read("1\n\nh 0.0 0.0 0.0")

    assert str(exception.value) == (
        "The first atom in one of the frames is not X. "
        "Please use PQ (default) md engine instead"
    )


def test_xyz_frame_reader_resets_topology_between_reads():
    reader = frame_reader.XYZFrameReader()
    topology = Topology(atoms=[Atom("h")])

    frame = reader.read("1\n\nh 0.0 0.0 0.0", topology=topology)
    assert frame.topology is topology

    frame = reader.read("1\n\no 1.0 2.0 3.0")

    assert frame.topology is not topology
    assert frame.atoms == [Atom("o")]


def test_xyz_frame_reader_defensive_invalid_format_branch(monkeypatch):
    class UnknownTrajectoryFormat:
        XYZ = object()
        VEL = object()
        FORCE = object()
        CHARGE = object()

        def __new__(cls, value):
            return object()

    monkeypatch.setattr(frame_reader, "TrajectoryFormat", UnknownTrajectoryFormat)

    reader = frame_reader.XYZFrameReader()
    monkeypatch.setattr(reader.logger, "error", lambda *args, **kwargs: None)

    with pytest.raises(FrameReaderError) as exception:
        reader.read("", traj_format="unknown")

    assert str(exception.value) == "Invalid TrajectoryFormat given. traj_format='unknown'"


class TestFrameReader:

    def test__read_header_line(self):
        reader = _FrameReader()

        with pytest.raises(FrameReaderError) as exception:
            reader._read_header_line("1 2.0 3.0")
        assert str(
            exception.value
        ) == "Invalid file format in header line of Frame."

        n_atoms, cell = reader._read_header_line("1 2.0 3.0 4.0 5.0 6.0 7.0")
        assert n_atoms == 1
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [5.0, 6.0, 7.0])

        n_atoms, cell = reader._read_header_line("2 2.0 3.0 4.0")
        assert n_atoms == 2
        assert np.allclose(cell.box_lengths, [2.0, 3.0, 4.0])
        assert np.allclose(cell.box_angles, [90.0, 90.0, 90.0])

        n_atoms, cell = reader._read_header_line("3")
        assert n_atoms == 3
        assert cell == Cell()

    def test__read_xyz(self):
        reader = _FrameReader()

        with pytest.raises(FrameReaderError) as exception:
            reader._read_xyz(["", "", "h 1.0 2.0 3.0", "o 2.0 2.0"], n_atoms=2)
        assert str(
            exception.value
        ) == "Invalid file format in xyz coordinates of Frame."

        xyz, atoms = reader._read_xyz(
            ["", "", "h 1.0 2.0 3.0", "o 2.0 2.0 2.0"], n_atoms=2)
        assert np.allclose(xyz, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert atoms == ["h", "o"]

    def test__read_scalar(self):
        reader = _FrameReader()

        with pytest.raises(FrameReaderError) as exception:
            reader._read_scalar(["", "", "h 1.0 2.0 3.0"], n_atoms=1)
        assert str(
            exception.value
        ) == "Invalid file format in scalar values of Frame."

        scalar, atoms = reader._read_scalar(["", "", "h 1.0"], n_atoms=1)
        assert np.allclose(scalar, [1.0])
        assert atoms == ["h"]

    def test_read(self):
        reader = _FrameReader()

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no 2.0 2.0 2.0"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [Atom(atom) for atom in ["h", "o"]]
        assert np.allclose(frame.pos, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.pos, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0",
            traj_format="vel"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.vel, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0 2.0 3.0\no1 2.0 2.0 2.0",
            traj_format="force"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.forces, [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

        frame = reader.read(
            "2 2.0 3.0 4.0 5.0 6.0 7.0\n\nh 1.0\no1 2.0",
            traj_format="charge"
        )
        assert frame.n_atoms == 2
        assert frame.atoms == [
            Atom(atom,
            use_guess_element=False) for atom in ["h",
            "o1"]
        ]
        assert np.allclose(frame.charges, [1.0, 2.0])
        assert frame.cell == Cell(2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

    def test_read_invalid_format(self):
        reader = _FrameReader()

        with pytest.raises(TrajectoryFormatError) as exception:
            reader.read("", traj_format="invalid")
        assert str(exception.value) == (
            "\n"
            "'invalid' is not a valid TrajectoryFormat.\n"
            f"Possible values are: {TrajectoryFormat.member_repr()} "
            "or their case insensitive string representation: "
            f"{TrajectoryFormat.value_repr()}"
        )

    def test__get_topology(self):
        reader = _FrameReader()

        topology = reader._get_topology(["h", "o"], None)
        assert topology.atoms == [Atom(atom) for atom in ["h", "o"]]

        topology = Topology()
        topology = reader._get_topology(["h", "o"], topology)
        assert topology == Topology()
