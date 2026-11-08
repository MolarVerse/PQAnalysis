"""
Integration tests of the traj2com command: the pqanalysis command line tool
on the TIP3P water box (examples/h2o_mm-mshake) and on UMCM-9
(examples/traj2comtraj), compared with the Python functions.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.io import (
    TrajectoryReader,
    read_restart_file,
    rst2xyz,
)
from PQAnalysis.tools import traj2com, traj_to_com_traj
from PQAnalysis.traj import Trajectory

pytestmark = pytest.mark.integration

N_MOLECULES = 2000


def run_cli(monkeypatch, *arguments):
    monkeypatch.setattr(
        sys, "argv", ["pqanalysis", "--log-file", "off", "traj2com", *arguments]
    )
    pqanalysis_main()


def read_com(filename):
    return TrajectoryReader(filename).read()


def moldescriptor_with_water_type(tmp_path):
    """The TIP3P moldescriptor with the WATER_TYPE and AMMONIA_TYPE lines."""

    text = open("moldescriptor.dat", encoding="utf-8").read()
    path = tmp_path / "moldescriptor_water_type.dat"
    path.write_text(
        "WATER_TYPE 1\nAMMONIA_TYPE 2\n\n" + text + "\n"
        "# Molecule Type 2\n"
        "Ammonia 4 0.0\n"
        "N 0 -0.8022\nH 1 0.2674\nH 1 0.2674\nH 1 0.2674\n",
        encoding="utf-8",
    )

    return str(path)


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestWaterBox:

    @pytest.fixture
    def water(self, test_integration_folder):
        rst2xyz("input.rst", output="water.xyz")
        system = read_restart_file("input.rst", "moldescriptor.dat")

        return Trajectory([system]), system

    def test_by_residue_equals_the_python_call(self, water, monkeypatch):
        traj, _ = water

        run_cli(
            monkeypatch, "water.xyz", "--by-residue",
            "--restart-file", "input.rst",
            "--moldescriptor-file", "moldescriptor.dat",
            "-o", "com.xyz",
        )

        expected = traj_to_com_traj(traj, by_residue=True)[0]
        written = read_com("com.xyz")[0]

        assert written.n_atoms == N_MOLECULES
        assert all(atom.name == "TIP3P-Water" for atom in written.atoms)
        assert written.cell == expected.cell
        difference = written.cell.image(written.pos - expected.pos)
        assert np.abs(difference).max() < 1e-5

    def test_every_center_of_mass_is_next_to_its_oxygen(self, water, monkeypatch):
        _, system = water

        run_cli(
            monkeypatch, "water.xyz", "--by-residue",
            "--restart-file", "input.rst",
            "--moldescriptor-file", "moldescriptor.dat",
            "-o", "com.xyz",
        )
        written = read_com("com.xyz")[0]

        difference = written.pos - system.pos[::3]
        distance = np.linalg.norm(written.cell.image(difference), axis=1)
        assert distance.min() > 0.05
        assert distance.max() < 0.08

    def test_groups_of_three_atoms_give_the_same_centers(self, water, monkeypatch):
        run_cli(
            monkeypatch, "water.xyz", "--group", "3", "-o", "groups.xyz"
        )
        run_cli(
            monkeypatch, "water.xyz", "--by-residue",
            "--restart-file", "input.rst",
            "--moldescriptor-file", "moldescriptor.dat",
            "-o", "residues.xyz",
        )

        groups = read_com("groups.xyz")[0]
        residues = read_com("residues.xyz")[0]

        assert groups.n_atoms == residues.n_atoms == N_MOLECULES
        assert np.abs(groups.cell.image(groups.pos - residues.pos)).max() < 1e-5

    def test_a_moldescriptor_with_water_type_lines(
        self, water, monkeypatch, tmp_path
    ):
        # WATER_TYPE and AMMONIA_TYPE only point at molecule types of the
        # file; the molecules are still defined by the entries
        with_types = moldescriptor_with_water_type(tmp_path)

        run_cli(
            monkeypatch, "water.xyz", "--by-residue",
            "--restart-file", "input.rst",
            "--moldescriptor-file", with_types,
            "-o", "com_types.xyz",
        )
        run_cli(
            monkeypatch, "water.xyz", "--by-residue",
            "--restart-file", "input.rst",
            "--moldescriptor-file", "moldescriptor.dat",
            "-o", "com.xyz",
        )

        with_types = read_com("com_types.xyz")[0]
        plain = read_com("com.xyz")[0]

        assert with_types.n_atoms == N_MOLECULES
        assert np.allclose(with_types.pos, plain.pos)
        assert [a.name for a in with_types.atoms] == [a.name for a in plain.atoms]

    def test_water_type_without_a_molecule_entry_defines_no_molecules(
        self, water, tmp_path
    ):
        only_type = tmp_path / "only_water_type.dat"
        only_type.write_text("WATER_TYPE 1\n", encoding="utf-8")

        with pytest.raises(PQValueError, match="no residues"):
            traj2com(
                ["water.xyz"],
                output="com.xyz",
                by_residue=True,
                restart_file="input.rst",
                moldescriptor_file=str(only_type),
            )

    def test_by_residue_needs_the_restart_file(self, water, monkeypatch):
        with pytest.raises(PQValueError, match="no residues"):
            run_cli(monkeypatch, "water.xyz", "--by-residue", "-o", "com.xyz")


@pytest.mark.parametrize("example_dir", ["traj2comtraj"], indirect=False)
class TestUMCM9:

    N_FRAMES = 3

    @pytest.fixture
    def trajectory_file(self, test_integration_folder):
        # the first frames only: groups of atoms are slow for 100 frames
        lines = open("umcm-9-md-01.xyz", encoding="utf-8").read().splitlines()
        n_lines = self.N_FRAMES * (int(lines[0].split()[0]) + 2)
        with open("umcm-9.xyz", "w", encoding="utf-8") as file:
            file.write("\n".join(lines[:n_lines]) + "\n")

        return "umcm-9.xyz"

    def test_selection_and_groups_equal_the_python_call(
        self, trajectory_file, monkeypatch
    ):
        # "*|X" selects everything but the dummy atom X of mass 999
        run_cli(
            monkeypatch, trajectory_file, "--selection", "*|X",
            "--group", "8", "-o", "com.xyz",
        )

        traj = TrajectoryReader(trajectory_file).read()
        expected = traj_to_com_traj(traj, selection=slice(1, None), group=8)
        written = read_com("com.xyz")

        assert len(written) == len(expected) == self.N_FRAMES
        for frame, reference in zip(written, expected):
            assert frame.n_atoms == reference.n_atoms == 101
            assert frame.combined_name == reference.combined_name
            difference = frame.cell.image(frame.pos - reference.pos)
            assert np.abs(difference).max() < 1e-5

    def test_one_center_of_mass_of_the_selection(self, trajectory_file, monkeypatch):
        run_cli(
            monkeypatch, trajectory_file, "--selection", "Zn", "-o", "zn.xyz"
        )

        written = read_com("zn.xyz")

        assert len(written) == self.N_FRAMES
        assert all(frame.n_atoms == 1 for frame in written)

    def test_output_to_stdout(self, trajectory_file, monkeypatch, capsys):
        run_cli(
            monkeypatch, trajectory_file, "--selection", "Zn",
        )

        assert capsys.readouterr().out.count("Zn") >= self.N_FRAMES
