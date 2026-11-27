"""
End-to-end tests of the trajectory conversion commands traj2box, traj2qmcfc
and traj2extxyz on the triclinic examples of examples/traj2box (216 atoms,
100 frames each, a box that changes from frame to frame).

The expected data is derived in the tests from the input files with plain
Python and NumPy, and the extended xyz output is additionally read back with
ASE, which does not share any code with PQAnalysis.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.io.exceptions import FileWritingModeError

pytestmark = pytest.mark.integration

FILES = ["acof_triclinic.xyz", "acof_triclinic_2.xyz"]
N_ATOMS = 216
N_FRAMES = 100


#: The standard reader parses xyz values as float32, so converted positions
#: agree with the input to about seven significant digits.
FLOAT32 = {"rtol": 2e-7, "atol": 2e-7}


def run_cli(monkeypatch, *arguments):
    monkeypatch.setattr(
        sys, "argv", ["pqanalysis", "--log-file", "off", *arguments]
    )
    pqanalysis_main()


def read_xyz(filename):
    """
    The frames of an xyz file: the box parameters (a, b, c, alpha, beta,
    gamma), the atom names and the positions.
    """

    lines = Path(filename).read_text(encoding="utf-8").splitlines()
    size = int(lines[0].split()[0]) + 2
    frames = []
    for start in range(0, len(lines) - size + 1, size):
        n_atoms = int(lines[start].split()[0])
        box = [float(value) for value in lines[start].split()[1:7]]
        atoms = [line.split() for line in lines[start + 2:start + 2 + n_atoms]]
        frames.append((
            np.array(box),
            [atom[0] for atom in atoms],
            np.array([[float(x) for x in atom[1:4]] for atom in atoms]),
        ))

    return frames


def lattice(box):
    """The lattice vectors (rows) of a, b, c, alpha, beta, gamma."""

    a, b, c = box[:3]
    alpha, beta, gamma = np.radians(box[3:])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)

    return np.array([
        [a, 0.0, 0.0],
        [b * np.cos(gamma), b * np.sin(gamma), 0.0],
        [cx, cy, np.sqrt(c**2 - cx**2 - cy**2)],
    ])


@pytest.fixture(name="frames")
def frames_fixture(test_integration_folder):
    return [frame for name in FILES for frame in read_xyz(name)]


@pytest.mark.parametrize("example_dir", ["traj2box"], indirect=False)
class TestTraj2Box:

    def test_the_box_of_every_frame(
        self, frames, monkeypatch
    ):
        run_cli(monkeypatch, "traj2box", *FILES, "-o", "box.dat")

        lines = Path("box.dat").read_text(encoding="utf-8").splitlines()

        assert len(lines) == 2 * N_FRAMES
        for index, (line, (box, _, _)) in enumerate(zip(lines, frames), 1):
            fields = line.split()
            assert int(fields[0]) == index
            assert np.allclose([float(v) for v in fields[1:]], box)
        # the box changes from frame to frame
        assert len({line.split()[1] for line in lines}) > 50

    def test_the_vmd_box_vertices(self, frames, monkeypatch):
        run_cli(monkeypatch, "traj2box", *FILES, "--vmd", "-o", "box.vmd.xyz")

        lines = Path("box.vmd.xyz").read_text(encoding="utf-8").splitlines()

        # every frame: the number of atoms (8), the box, and 8 vertices
        assert len(lines) == 2 * N_FRAMES * 10
        for index, (box, _, _) in enumerate(frames):
            block = lines[10 * index:10 * index + 10]
            vectors = lattice(box)
            expected = sorted(
                tuple(np.round(
                    i * vectors[0] + j * vectors[1] + k * vectors[2]
                    - 0.5 * vectors.sum(axis=0), 6
                ))
                for i in (0, 1) for j in (0, 1) for k in (0, 1)
            )
            vertices = sorted(
                tuple(np.round([float(v) for v in row.split()[1:4]], 6))
                for row in block[2:]
            )

            assert block[0].split()[0] == "8"
            assert all(row.split()[0] == "X" for row in block[2:])
            assert np.allclose(expected, vertices, atol=1e-5)

    def test_the_output_goes_to_stdout_without_a_file(
        self, frames, monkeypatch, capsys
    ):
        run_cli(monkeypatch, "traj2box", FILES[0])

        rows = [
            line.split() for line in capsys.readouterr().out.splitlines()
            if line.split() and line.split()[0].isdigit()
        ]

        assert len(rows) == N_FRAMES
        assert np.allclose([float(v) for v in rows[0][1:]], frames[0][0])

    def test_the_writing_modes(self, frames, monkeypatch):
        run_cli(monkeypatch, "traj2box", FILES[0], "-o", "box.dat")
        first = Path("box.dat").read_text(encoding="utf-8")

        # the default refuses to replace an existing file
        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "traj2box", FILES[0], "-o", "box.dat")
        assert Path("box.dat").read_text(encoding="utf-8") == first

        run_cli(monkeypatch, "traj2box", FILES[0], "-o", "box.dat", "--mode", "a")
        assert len(Path("box.dat").read_text().splitlines()) == 2 * N_FRAMES

        run_cli(monkeypatch, "traj2box", FILES[0], "-o", "box.dat", "--mode", "o")
        assert Path("box.dat").read_text(encoding="utf-8") == first


@pytest.mark.parametrize("example_dir", ["traj2box"], indirect=False)
class TestTraj2QMCFC:

    def test_a_dummy_atom_is_added_to_every_frame(self, frames, monkeypatch):
        run_cli(monkeypatch, "traj2qmcfc", *FILES, "-o", "traj.qmcfc.xyz")

        converted = read_xyz("traj.qmcfc.xyz")

        assert len(converted) == len(frames) == 2 * N_FRAMES
        for (box, names, positions), (new_box, new_names, new_positions) in zip(
            frames, converted
        ):
            assert np.allclose(new_box, box)
            assert new_names == ["X", *names]
            assert np.allclose(new_positions[0], 0.0)
            assert np.allclose(new_positions[1:], positions, **FLOAT32)

    def test_the_round_trip_through_the_extended_xyz(self, frames, monkeypatch):
        run_cli(monkeypatch, "traj2qmcfc", FILES[0], "-o", "traj.qmcfc.xyz")
        run_cli(
            monkeypatch, "traj2extxyz", "traj.qmcfc.xyz", "--engine", "QMCFC",
            "-o", "back.extxyz",
        )

        ase_io = pytest.importorskip("ase.io")
        atoms = ase_io.read("back.extxyz", index=":")

        # the QMCFC engine drops the dummy atom again
        assert len(atoms) == N_FRAMES
        assert all(len(frame) == N_ATOMS for frame in atoms)
        assert np.allclose(atoms[10].positions, frames[10][2], **FLOAT32)

    def test_the_writing_modes(self, frames, monkeypatch):
        run_cli(monkeypatch, "traj2qmcfc", FILES[0], "-o", "q.xyz")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "traj2qmcfc", FILES[0], "-o", "q.xyz")

        run_cli(monkeypatch, "traj2qmcfc", FILES[0], "-o", "q.xyz", "--mode", "a")
        assert len(read_xyz("q.xyz")) == 2 * N_FRAMES


@pytest.mark.parametrize("example_dir", ["traj2box"], indirect=False)
class TestTraj2ExtXYZ:

    def test_read_back_with_ase(self, frames, monkeypatch):
        ase_io = pytest.importorskip("ase.io")

        run_cli(monkeypatch, "traj2extxyz", *FILES, "-o", "traj.extxyz")

        atoms = ase_io.read("traj.extxyz", index=":")

        assert len(atoms) == 2 * N_FRAMES
        for frame, (box, names, positions) in zip(atoms, frames):
            assert frame.get_chemical_symbols() == names
            assert np.allclose(frame.positions, positions, **FLOAT32)
            assert np.allclose(frame.cell.array, lattice(box), rtol=1e-9, atol=1e-8)
            assert frame.pbc.all()

    @pytest.mark.parametrize("profile", ["ase", "pq", "nep"])
    def test_the_profiles_keep_the_geometry(self, frames, monkeypatch, profile):
        ase_io = pytest.importorskip("ase.io")

        run_cli(
            monkeypatch, "traj2extxyz", FILES[0], "--extxyz-profile", profile,
            "-o", f"{profile}.extxyz",
        )

        atoms = ase_io.read(f"{profile}.extxyz", index=":")

        assert len(atoms) == N_FRAMES
        assert np.allclose(
            atoms[-1].positions, frames[N_FRAMES - 1][2], **FLOAT32
        )

        # the nep profile writes the lower case keys that GPUMD expects, which
        # ASE does not read as a cell; the lattice is parsed from the line
        header = Path(f"{profile}.extxyz").read_text(
            encoding="utf-8").splitlines()[1]
        values = header.split('attice="')[1].split('"')[0].split()
        assert np.allclose(
            np.array(values, dtype=float).reshape(3, 3),
            lattice(frames[0][0]), rtol=1e-9, atol=1e-8,
        )
        if profile != "nep":
            assert np.allclose(
                atoms[0].cell.array, lattice(frames[0][0]), rtol=1e-9, atol=1e-8
            )

    def test_the_nep_profile_writes_lower_case_keys(self, frames, monkeypatch):
        run_cli(
            monkeypatch, "traj2extxyz", FILES[0], "--extxyz-profile", "nep",
            "-o", "nep.extxyz",
        )
        run_cli(
            monkeypatch, "traj2extxyz", FILES[0], "--extxyz-profile", "ase",
            "-o", "ase.extxyz",
        )

        nep = Path("nep.extxyz").read_text(encoding="utf-8").splitlines()[1]
        ase = Path("ase.extxyz").read_text(encoding="utf-8").splitlines()[1]

        # GPUMD does not strip quotes around the properties value
        assert nep.startswith('lattice="') and " properties=species" in nep
        assert ase.startswith('Lattice="') and 'Properties="' in ase

    def test_the_writing_modes(self, frames, monkeypatch):
        run_cli(monkeypatch, "traj2extxyz", FILES[0], "-o", "e.extxyz")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "traj2extxyz", FILES[0], "-o", "e.extxyz")
        first = Path("e.extxyz").read_text(encoding="utf-8")

        run_cli(monkeypatch, "traj2extxyz", FILES[0], "-o", "e.extxyz", "--mode", "o")
        assert Path("e.extxyz").read_text(encoding="utf-8") == first
