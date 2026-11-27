"""
End-to-end tests of the structure conversion commands rst2xyz, xyz2rst,
xyz2gen and gen2xyz through the real command line, on the water example
(one H2O), the 2000 molecule TIP3P water box of examples/h2o_mm-mshake and
the triclinic acof example.

The expected data comes from parsing the input files in the tests with plain
Python and NumPy. Round trips check the commands against each other; the
2000 molecule box also checks that the moldescriptor file reproduces the
molecule columns of a real restart file for every atom.

The standard reader parses xyz values as float32, so everything that goes
through an xyz file is compared to about seven significant digits.
"""

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.io.exceptions import FileWritingModeError

pytestmark = pytest.mark.integration

FLOAT32 = {"rtol": 2e-7, "atol": 2e-7}
N_WATERS = 2000


def run_cli(monkeypatch, *arguments):
    monkeypatch.setattr(
        sys, "argv", ["pqanalysis", "--log-file", "off", *arguments]
    )
    pqanalysis_main()


def read_restart(filename):
    """The box and the rows (name, integer columns, floats) of a restart."""

    box, rows = None, []
    for line in Path(filename).read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or fields[0] == "Step":
            continue
        if fields[0] == "Box":
            box = np.array([float(value) for value in fields[1:7]])
            continue
        rows.append((
            fields[0],
            [int(value) for value in fields[1:3]],
            np.array([float(value) for value in fields[3:]]),
        ))

    return box, rows


def read_xyz(filename):
    """The frames (box or None, names, positions) of an xyz file."""

    lines = Path(filename).read_text(encoding="utf-8").splitlines()
    frames, start = [], 0
    while start < len(lines):
        header = lines[start].split()
        n_atoms = int(header[0])
        box = np.array([float(v) for v in header[1:7]]) if len(header) > 1 else None
        atoms = [line.split() for line in lines[start + 2:start + 2 + n_atoms]]
        frames.append((
            box,
            [atom[0] for atom in atoms],
            np.array([[float(x) for x in atom[1:4]] for atom in atoms]),
        ))
        start += n_atoms + 2

    return frames


def read_gen(filename):
    """
    A gen file: the periodicity flag, the species, the (species index,
    positions) of the atoms, and the origin and lattice vectors if periodic.
    """

    lines = Path(filename).read_text(encoding="utf-8").splitlines()
    n_atoms, flag = int(lines[0].split()[0]), lines[0].split()[1]
    species = lines[1].split()
    atoms = [lines[2 + i].split() for i in range(n_atoms)]
    lattice = None
    if flag == "S":
        lattice = np.array([
            [float(v) for v in lines[2 + n_atoms + i].split()] for i in range(4)
        ])

    return flag, species, atoms, lattice


def lattice_vectors(box):
    a, b, c = box[:3]
    alpha, beta, gamma = np.radians(box[3:])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)

    return np.array([
        [a, 0.0, 0.0],
        [b * np.cos(gamma), b * np.sin(gamma), 0.0],
        [cx, cy, np.sqrt(c**2 - cx**2 - cy**2)],
    ])


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestRst2XyzSmall:

    def test_the_positions_and_the_box(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "rst2xyz", "restart.rst", "-o", "out.xyz")

        box, rows = read_restart("restart.rst")
        (xyz_box, names, positions), = read_xyz("out.xyz")

        assert np.allclose(xyz_box, box)
        assert names == [row[0] for row in rows]
        assert np.allclose(positions, [row[2][:3] for row in rows], atol=1e-9, rtol=0)

    def test_nobox(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "rst2xyz", "restart.rst", "--nobox", "-o", "out.xyz")

        lines = Path("out.xyz").read_text(encoding="utf-8").splitlines()

        assert lines[0].split() == ["3"]
        assert len(read_xyz("out.xyz")[0][1]) == 3

    def test_stdout(self, test_integration_folder, monkeypatch, capsys):
        run_cli(monkeypatch, "rst2xyz", "restart.rst")

        out = capsys.readouterr().out.splitlines()
        atoms = [line.split() for line in out if line.split()[:1] in (["O"], ["H"])]

        assert [atom[0] for atom in atoms] == ["O", "H", "H"]
        assert np.allclose([float(v) for v in atoms[1][1:]], [5.96, 5.0, 5.0])

    def test_the_writing_modes(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "rst2xyz", "restart.rst", "-o", "out.xyz")
        first = Path("out.xyz").read_text(encoding="utf-8")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "rst2xyz", "restart.rst", "-o", "out.xyz")

        run_cli(monkeypatch, "rst2xyz", "restart.rst", "-o", "out.xyz", "--mode", "a")
        assert len(read_xyz("out.xyz")) == 2

        run_cli(monkeypatch, "rst2xyz", "restart.rst", "-o", "out.xyz", "--mode", "o")
        assert Path("out.xyz").read_text(encoding="utf-8") == first


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestWaterBox:

    def test_rst2xyz(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "rst2xyz", "input.rst", "-o", "box.xyz")

        box, rows = read_restart("input.rst")
        (xyz_box, names, positions), = read_xyz("box.xyz")

        assert np.allclose(xyz_box, box)
        assert names == ["O", "H", "H"] * N_WATERS
        assert np.allclose(positions, [row[2][:3] for row in rows], atol=1e-9, rtol=0)

    def test_the_round_trip_restores_the_restart_columns(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(monkeypatch, "rst2xyz", "input.rst", "-o", "box.xyz")
        run_cli(
            monkeypatch, "xyz2rst", "box.xyz", "--moldescriptor_file",
            "moldescriptor.dat", "-o", "back.rst",
        )

        box, rows = read_restart("input.rst")
        back_box, back_rows = read_restart("back.rst")

        assert np.allclose(back_box, box)
        assert len(back_rows) == len(rows) == 3 * N_WATERS
        # the names, the index in the molecule and the molecule type of every
        # atom are inferred from the moldescriptor file
        assert [r[:2] for r in back_rows] == [r[:2] for r in rows]
        assert np.allclose(
            [r[2][:3] for r in back_rows], [r[2][:3] for r in rows], **FLOAT32
        )

    def test_xyz2rst_without_a_moldescriptor_has_no_molecules(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(monkeypatch, "rst2xyz", "input.rst", "-o", "box.xyz")
        run_cli(monkeypatch, "xyz2rst", "box.xyz", "-o", "plain.rst")

        _, rows = read_restart("plain.rst")

        assert {row[1][1] for row in rows} == {0}

    def test_randomize_displaces_the_positions(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(monkeypatch, "rst2xyz", "input.rst", "-o", "box.xyz")
        run_cli(monkeypatch, "xyz2rst", "box.xyz", "-o", "plain.rst")
        run_cli(
            monkeypatch, "xyz2rst", "box.xyz", "--randomize", "0.1", "-o", "rnd.rst"
        )
        run_cli(
            monkeypatch, "xyz2rst", "box.xyz", "--randomize", "0.1", "-o", "rnd2.rst"
        )

        plain = np.array([r[2][:3] for r in read_restart("plain.rst")[1]])
        random = np.array([r[2][:3] for r in read_restart("rnd.rst")[1]])
        again = np.array([r[2][:3] for r in read_restart("rnd2.rst")[1]])
        displacement = random - plain

        # normally distributed with a standard deviation of 0.1 A
        assert displacement.std() == pytest.approx(0.1, rel=0.03)
        assert abs(displacement.mean()) < 0.005
        # the default seed of the conversion is fixed, so a rerun is the same
        assert np.array_equal(random, again)
        # without the option nothing moves
        assert np.abs(plain).max() > 1.0


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestXyz2Rst:

    def test_the_last_frame_with_velocities(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(
            monkeypatch, "xyz2rst", "trajectory.xyz", "--velocity_file",
            "trajectory.vel", "-o", "out.rst",
        )

        box, rows = read_restart("out.rst")
        positions = read_xyz("trajectory.xyz")[-1][2]
        velocities = read_xyz("trajectory.vel")[-1][2]

        assert len(read_xyz("trajectory.xyz")) == 25
        assert np.allclose(box, [10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
        assert [row[0] for row in rows] == ["O", "H", "H"]
        assert np.allclose([r[2][:3] for r in rows], positions, **FLOAT32)
        assert np.allclose([r[2][3:6] for r in rows], velocities, **FLOAT32)
        # no forces were given
        assert np.all(np.array([r[2][6:9] for r in rows]) == 0.0)

    def test_the_forces(self, test_integration_folder, monkeypatch):
        forces = np.arange(9.0).reshape(3, 3) * 1.5
        lines = ["3 10.0 10.0 10.0 90.0 90.0 90.0", ""]
        lines += [
            f"{name} {x} {y} {z}"
            for name, (x, y, z) in zip("OHH", forces)
        ]
        Path("trajectory.force").write_text("\n".join(lines) + "\n", encoding="utf-8")

        run_cli(
            monkeypatch, "xyz2rst", "trajectory.xyz", "--force_file",
            "trajectory.force", "-o", "out.rst",
        )

        _, rows = read_restart("out.rst")

        assert np.allclose([r[2][6:9] for r in rows], forces, **FLOAT32)

    def test_the_moldescriptor_defines_the_molecule(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(
            monkeypatch, "xyz2rst", "trajectory.xyz", "--moldescriptor_file",
            "moldescriptor.dat", "-o", "out.rst",
        )

        _, rows = read_restart("out.rst")

        # the index within the molecule and the molecule type
        assert [r[1] for r in rows] == [[1, 1], [2, 1], [3, 1]]

    def test_the_qmcfc_engine_round_trip(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "rst2xyz", "restart.rst", "-o", "r.xyz")
        run_cli(monkeypatch, "traj2qmcfc", "r.xyz", "-o", "r.qmcfc.xyz")
        run_cli(
            monkeypatch, "xyz2rst", "r.qmcfc.xyz", "--engine", "QMCFC", "-o", "q.rst"
        )
        run_cli(
            monkeypatch, "rst2xyz", "q.rst", "--engine", "QMCFC", "-o", "back.xyz"
        )

        original = read_xyz("r.qmcfc.xyz")[0]
        back = read_xyz("back.xyz")[0]

        # the QMCFC restart file carries the positions twice (the second
        # block is the previous step), the dummy atom X comes first
        _, rows = read_restart("q.rst")
        assert len(rows[0][2]) == 18
        assert back[1] == original[1] == ["X", "O", "H", "H"]
        assert np.allclose(back[2], original[2], **FLOAT32)

    def test_the_writing_modes(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "xyz2rst", "trajectory.xyz", "-o", "out.rst")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "xyz2rst", "trajectory.xyz", "-o", "out.rst")


@pytest.mark.parametrize("example_dir", ["traj2box"], indirect=False)
class TestGen:

    @pytest.fixture(autouse=True)
    def single_frame(self, test_integration_folder):
        lines = Path("acof_triclinic.xyz").read_text(encoding="utf-8").splitlines()
        Path("ac.xyz").write_text("\n".join(lines[:218]) + "\n", encoding="utf-8")

    def test_xyz2gen_of_a_periodic_structure(self, monkeypatch):
        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "-o", "ac.gen")

        box, names, positions = read_xyz("ac.xyz")[0]
        flag, species, atoms, lattice = read_gen("ac.gen")

        assert flag == "S"
        # the species in the order of their first appearance
        assert species == list(dict.fromkeys(names))
        assert [species[int(atom[1]) - 1] for atom in atoms] == names
        assert [int(atom[0]) for atom in atoms] == list(range(1, len(names) + 1))
        assert np.allclose(
            [[float(v) for v in atom[2:5]] for atom in atoms], positions, **FLOAT32
        )
        # the origin and the three lattice vectors of the triclinic box
        assert np.allclose(lattice[0], 0.0)
        assert np.allclose(lattice[1:], lattice_vectors(box), atol=1e-6)

    def test_a_cluster_has_no_lattice(self, monkeypatch):
        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "--periodic", "False", "-o", "c.gen")

        flag, _, atoms, lattice = read_gen("c.gen")

        assert flag == "C"
        assert lattice is None
        assert len(atoms) == 216
        assert len(Path("c.gen").read_text(encoding="utf-8").splitlines()) == 2 + 216

    def test_the_round_trip(self, monkeypatch):
        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "-o", "ac.gen")
        run_cli(monkeypatch, "gen2xyz", "ac.gen", "-o", "back.xyz")

        box, names, positions = read_xyz("ac.xyz")[0]
        back_box, back_names, back_positions = read_xyz("back.xyz")[0]

        assert back_names == names
        assert np.allclose(back_positions, positions, **FLOAT32)
        assert np.allclose(back_box, box, rtol=1e-9, atol=1e-8)

    def test_the_round_trip_of_a_cluster_has_no_box(self, monkeypatch):
        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "--periodic", "False", "-o", "c.gen")
        run_cli(monkeypatch, "gen2xyz", "c.gen", "-o", "back.xyz")

        _, names, positions = read_xyz("ac.xyz")[0]
        back_box, back_names, back_positions = read_xyz("back.xyz")[0]

        assert back_box is None
        assert back_names == names
        assert np.allclose(back_positions, positions, **FLOAT32)

    def test_gen2xyz_without_the_box(self, monkeypatch):
        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "-o", "ac.gen")
        run_cli(monkeypatch, "gen2xyz", "ac.gen", "--nobox", "-o", "nobox.xyz")

        assert Path("nobox.xyz").read_text(encoding="utf-8").splitlines()[0] == "216"

    def test_the_writing_modes(self, monkeypatch):
        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "-o", "ac.gen")
        first = Path("ac.gen").read_text(encoding="utf-8")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "xyz2gen", "ac.xyz", "-o", "ac.gen")

        run_cli(monkeypatch, "xyz2gen", "ac.xyz", "-o", "ac.gen", "--mode", "o")
        assert Path("ac.gen").read_text(encoding="utf-8") == first

        run_cli(monkeypatch, "gen2xyz", "ac.gen", "-o", "back.xyz")
        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "gen2xyz", "ac.gen", "-o", "back.xyz")
