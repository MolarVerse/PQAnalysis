import builtins
import importlib

import numpy as np
import pytest

from PQAnalysis.core import Cell, Atom
from PQAnalysis.io import ExtXYZFrameReader, _FrameReader
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


def test_get_frame_reader_extxyz():
    reader = frame_reader.get_frame_reader(TrajectoryFormat.EXTXYZ)

    assert isinstance(reader, frame_reader.ExtXYZFrameReader)
    assert isinstance(reader, ExtXYZFrameReader)


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


def test_extxyz_frame_reader_reads_properties_and_metadata():
    reader = frame_reader.ExtXYZFrameReader()
    frame_string = (
        '2\n'
        'Lattice="10 0 0 0 20 0 0 0 30" '
        'Properties=species:S:1:pos:R:3:forces:R:3:charge:R:1 '
        'energy=-1.5 virial="1 0 0 0 2 0 0 0 3" '
        'stress="4 0 0 0 5 0 0 0 6"\n'
        'H 0 1 2 0.1 0.2 0.3 -0.1\n'
        'O 3 4 5 0.4 0.5 0.6 -0.2'
    )

    frame = reader.read(frame_string)

    assert frame.n_atoms == 2
    assert frame.atoms == [Atom("h"), Atom("o")]
    assert np.allclose(frame.pos, [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    assert np.allclose(frame.forces, [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    assert np.allclose(frame.charges, [-0.1, -0.2])
    assert frame.energy == -1.5
    assert np.allclose(frame.virial, np.diag([1.0, 2.0, 3.0]))
    assert np.allclose(frame.stress, np.diag([4.0, 5.0, 6.0]))
    assert frame.cell == Cell(10.0, 20.0, 30.0)


def test_extxyz_frame_reader_reads_nep_writer_style_header():
    reader = frame_reader.ExtXYZFrameReader()
    frame_string = (
        '1\n'
        'energy=-2.0 config_type=nep2xyz '
        'lattice="2 0 0 0 3 0 0 0 4 " '
        'properties=species:S:1:pos:R:3:forces:R:3\n'
        'C 0.0 0.5 1.0 1.0 2.0 3.0'
    )

    frame = reader.read(frame_string)

    assert frame.atoms == [Atom("c")]
    assert np.allclose(frame.pos, [[0.0, 0.5, 1.0]])
    assert np.allclose(frame.forces, [[1.0, 2.0, 3.0]])
    assert frame.energy == -2.0
    assert frame.cell == Cell(2.0, 3.0, 4.0)


def test_extxyz_frame_reader_reads_spaced_metadata_assignment():
    reader = frame_reader.ExtXYZFrameReader()
    frame_string = (
        '1\n'
        'Lattice = "4 0 0 0 5 0 0 0 6" '
        'Properties=species:S:1:pos:R:3:force:R:3 Energy= -1.0\n'
        'H 0.0 0.5 1.0 1.0 2.0 3.0'
    )

    frame = reader.read(frame_string)

    assert frame.atoms == [Atom("h")]
    assert np.allclose(frame.pos, [[0.0, 0.5, 1.0]])
    assert np.allclose(frame.forces, [[1.0, 2.0, 3.0]])
    assert frame.energy == -1.0
    assert frame.cell == Cell(4.0, 5.0, 6.0)


def test_extxyz_frame_reader_reads_noncanonical_lattice_vectors():
    reader = frame_reader.ExtXYZFrameReader()
    lattice = np.array(
        [
            [0.0, 0.0, 5.785795211791992],
            [5.928857803344727, 0.0, 0.0],
            [-2.9644289016723633, 8.167890548706055, -2.892897605895996],
        ]
    )
    frame_string = (
        '1\n'
        'Lattice="0.0 0.0 5.785795211791992 '
        '5.928857803344727 0.0 0.0 '
        '-2.9644289016723633 8.167890548706055 -2.892897605895996" '
        'Properties=species:S:1:pos:R:3\n'
        'H 0.0 0.0 0.0'
    )

    frame = reader.read(frame_string)

    assert np.isclose(frame.cell.volume, abs(np.linalg.det(lattice)))


def test_extxyz_frame_reader_rejects_degenerate_lattice_vectors():
    reader = frame_reader.ExtXYZFrameReader()
    frame_string = (
        '1\n'
        'Lattice="0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0" '
        'Properties=species:S:1:pos:R:3\n'
        'H 0.0 0.0 0.0'
    )

    with pytest.raises(FrameReaderError) as exception:
        reader.read(frame_string)

    assert str(exception.value) == "Invalid Lattice metadata in extended xyz frame."


def test_extxyz_frame_reader_preserves_qmcfc_handling():
    reader = frame_reader.ExtXYZFrameReader(md_format="qmcfc")
    frame_string = (
        '2\n'
        'Properties=species:S:1:pos:R:3:forces:R:3:charge:R:1\n'
        'X 0 0 0 0 0 0 0.0\n'
        'H 1 2 3 0.1 0.2 0.3 -0.1'
    )

    frame = reader.read(frame_string)

    assert frame.n_atoms == 1
    assert frame.atoms == [Atom("h")]
    assert np.allclose(frame.pos, [[1.0, 2.0, 3.0]])
    assert np.allclose(frame.forces, [[0.1, 0.2, 0.3]])
    assert np.allclose(frame.charges, [-0.1])


def test_extxyz_frame_reader_rejects_invalid_format():
    reader = frame_reader.ExtXYZFrameReader()

    with pytest.raises(FrameReaderError) as exception:
        reader.read("", traj_format=TrajectoryFormat.XYZ)

    assert str(exception.value) == (
        "Invalid TrajectoryFormat given. "
        "traj_format=<TrajectoryFormat.XYZ: 'XYZ'>"
    )


@pytest.mark.parametrize(
    ("frame_string", "expected", "md_format"),
    [
        (
            "not-an-int\nProperties=species:S:1:pos:R:3\nH 0.0 0.0 0.0",
            "Invalid number of atoms in extended xyz frame.",
            "pq",
        ),
        (
            "1\ncomment\nH 0.0 0.0 0.0",
            "Extended xyz frame does not define Properties metadata.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R\nH 0.0 0.0 0.0",
            "Invalid Properties metadata in extended xyz frame.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R:three\nH 0.0 0.0 0.0",
            "Invalid Properties metadata in extended xyz frame.",
            "pq",
        ),
        (
            "1\nProperties=pos:R:3\n0.0 0.0 0.0",
            "Extended xyz frame requires species and pos Properties.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R:2\nH 0.0 0.0",
            "Invalid Properties metadata in extended xyz frame.",
            "pq",
        ),
        (
            (
                "1\nProperties=species:S:1:pos:R:3:forces:R:2\n"
                "H 0.0 0.0 0.0 1.0 2.0"
            ),
            "Invalid Properties metadata in extended xyz frame.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R:3",
            "Invalid atom line in extended xyz frame.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R:3\nH 0.0 0.0",
            "Invalid atom line in extended xyz frame.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R:3\nH 0.0 0.0 0.0",
            (
                "The first atom in one of the frames is not X. "
                "Please use PQ (default) md engine instead"
            ),
            "qmcfc",
        ),
        (
            "1\nProperties=species:S:1:pos:R:3\nH 0.0 nope 0.0",
            "Invalid numeric value in extended xyz frame.",
            "pq",
        ),
        (
            (
                '1\nLattice="1 2" Properties=species:S:1:pos:R:3\n'
                "H 0.0 0.0 0.0"
            ),
            "Invalid Lattice metadata in extended xyz frame.",
            "pq",
        ),
        (
            "1\nProperties=species:S:1:pos:R:3 energy=nope\nH 0.0 0.0 0.0",
            "Invalid energy metadata in extended xyz frame.",
            "pq",
        ),
        (
            (
                '1\nProperties=species:S:1:pos:R:3 virial="1 2"\n'
                "H 0.0 0.0 0.0"
            ),
            "Invalid virial metadata in extended xyz frame.",
            "pq",
        ),
    ],
)
def test_extxyz_frame_reader_rejects_invalid_frames(
    frame_string,
    expected,
    md_format,
):
    reader = frame_reader.ExtXYZFrameReader(md_format=md_format)

    with pytest.raises(FrameReaderError) as exception:
        reader.read(frame_string)

    assert str(exception.value) == expected


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
