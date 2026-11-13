"""
Integration tests of the traj2com command: the pqanalysis command line tool
on the TIP3P water box (examples/h2o_mm-mshake) and on UMCM-9
(examples/traj2comtraj), compared with the Python functions.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.topology import Topology
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.io import (
    TrajectoryReader,
    TrajectoryWriter,
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


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestWaterBoxVelocities:

    """
    A vel file of the water box (the restart file has the velocities): the
    result is the velocity of the center of mass of every molecule.
    """

    @pytest.fixture
    def velocities(self, test_integration_folder):
        system = read_restart_file("input.rst", "moldescriptor.dat")
        TrajectoryWriter("water.vel").write(
            Trajectory([system, system]), traj_type="vel"
        )
        masses = system.atomic_masses.reshape(-1, 3)
        expected = (
            system.vel.reshape(-1, 3, 3) * masses[:, :, None]
        ).sum(1) / masses.sum(1)[:, None]

        return expected

    def test_by_residue_gives_the_center_of_mass_velocities(
        self, velocities, monkeypatch
    ):
        run_cli(
            monkeypatch, "water.vel", "--by-residue",
            "--restart-file", "input.rst",
            "--moldescriptor-file", "moldescriptor.dat",
            "-o", "com.vel",
        )

        written = read_com("com.vel")

        assert len(written) == 2
        for frame in written:
            assert frame.n_atoms == N_MOLECULES
            assert frame.has_vel and not frame.has_pos
            assert np.allclose(
                frame.vel, velocities, rtol=1e-5, atol=1e-5 * np.abs(velocities).max()
            )

    def test_groups_of_three_atoms_give_the_same_velocities(
        self, velocities, monkeypatch
    ):
        run_cli(monkeypatch, "water.vel", "--group", "3", "-o", "groups.vel")

        written = read_com("groups.vel")

        assert np.allclose(
            written[0].vel, velocities, rtol=1e-5,
            atol=1e-5 * np.abs(velocities).max(),
        )

    def test_the_velocity_of_the_center_of_all_atoms_is_the_total_momentum(
        self, test_integration_folder, monkeypatch
    ):
        system = read_restart_file("input.rst", "moldescriptor.dat")
        TrajectoryWriter("water.vel").write(Trajectory([system]), traj_type="vel")

        run_cli(monkeypatch, "water.vel", "-o", "total.vel")

        # the single center is named after all 6000 atoms, which is too long
        # for the reader, so the numbers are taken from the line itself
        line = open("total.vel", encoding="utf-8").read().splitlines()[2]
        velocity = np.array([float(value) for value in line.split()[-3:]])
        masses = system.atomic_masses
        momentum = (system.vel * masses[:, None]).sum(0)

        assert np.allclose(
            velocity * masses.sum(), momentum,
            rtol=1e-5, atol=1e-5 * np.abs(momentum).max(),
        )

    def test_continuous_needs_positions(self, velocities, monkeypatch):
        with pytest.raises(PQValueError, match="needs the positions"):
            run_cli(
                monkeypatch, "water.vel", "--group", "3", "--continuous",
                "-o", "x.vel",
            )


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestContinuousCli:

    def test_the_option_follows_the_centers_over_the_boundary(self, test_integration_folder, monkeypatch):
        # one water that moves 0.3 of the box per frame along x
        system = read_restart_file("input.rst", "moldescriptor.dat")
        box = system.cell.box_lengths[0]
        frames = []
        for frame in range(6):
            moved = system.pos.copy()
            moved[:, 0] += 0.3 * box * frame
            frames.append(
                AtomicSystem(
                    topology=system.topology, pos=system.cell.image(moved),
                    cell=system.cell,
                )
            )
        TrajectoryWriter("moving.xyz").write(Trajectory(frames))
        arguments = [
            "moving.xyz", "--by-residue", "--restart-file", "input.rst",
            "--moldescriptor-file", "moldescriptor.dat",
        ]

        run_cli(monkeypatch, *arguments, "-o", "folded.xyz")
        run_cli(monkeypatch, *arguments, "--continuous", "-o", "continuous.xyz")

        folded = np.array([f.pos[0] for f in read_com("folded.xyz")])
        continuous = np.array([f.pos[0] for f in read_com("continuous.xyz")])

        assert np.diff(continuous[:, 0]).min() == pytest.approx(0.3 * box, abs=1e-4)
        assert np.abs(np.diff(folded[:, 0])).max() > 0.5 * box
        assert np.allclose(continuous[0], folded[0])


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


class TestMoleculeTypes:

    """
    12 waters, 5 ammonias and 3 ions in a 15 A box, with a moldescriptor file
    that names the water and the ammonia type: the centers of one kind of
    molecules through the command line.
    """

    MOLDESCRIPTOR = (
        "WATER_TYPE 2\nAMMONIA_TYPE 3\n\n"
        "Sodium 1 1.0\nNa 0 1.0\n\n"
        "H2O 3 0.0\nO 0 -0.8\nH 1 0.4\nH 1 0.4\n\n"
        "NH3 4 0.0\nN 0 -0.8\nH 1 0.27\nH 1 0.27\nH 1 0.27\n"
    )
    KINDS = [("H2O", 2, "OHH", 12), ("NH3", 3, "NHHH", 5), ("Sodium", 1, ("Na",), 3)]

    @pytest.fixture
    def mixed(self, tmp_path, monkeypatch):
        from PQAnalysis.core import Cell
        from PQAnalysis.io import MoldescriptorReader
        from PQAnalysis.io.restart_file.api import write_restart_file

        monkeypatch.chdir(tmp_path)
        (tmp_path / "moldescriptor.dat").write_text(self.MOLDESCRIPTOR, encoding="utf-8")
        residues = MoldescriptorReader("moldescriptor.dat").read()
        by_name = {residue.name: residue for residue in residues}

        rng = np.random.default_rng(4)
        atoms, residue_ids, centers = [], [], []
        for name, number, symbols, count in self.KINDS:
            for _ in range(count):
                atoms += [Atom(symbol) for symbol in symbols]
                residue_ids += [number] * len(symbols)
                centers.append(rng.uniform(0.0, 15.0, 3))
        positions = np.vstack([
            c + 0.5 * rng.normal(size=(n, 3))
            for c, n in zip(centers, [len(s) for _, _, s, k in self.KINDS for _ in range(k)])
        ])
        topology = Topology(
            atoms=atoms, residue_ids=np.array(residue_ids),
            reference_residues=residues,
        )
        system = AtomicSystem(
            topology=topology, pos=Cell(15, 15, 15).image(positions),
            cell=Cell(15, 15, 15),
        )
        TrajectoryWriter("mixed.xyz").write(Trajectory([system]))
        write_restart_file(system, "mixed.rst")

        return system, by_name

    def run(self, monkeypatch, *types, output="out.xyz"):
        arguments = [
            "mixed.xyz", "--by-residue", "--restart-file", "mixed.rst",
            "--moldescriptor-file", "moldescriptor.dat", "-o", output,
        ]
        for molecule_type in types:
            arguments += ["--molecule-type", molecule_type]
        run_cli(monkeypatch, *arguments)

        return read_com(output)[0]

    @pytest.mark.parametrize(
        "types, names",
        [
            (["water"], ["H2O"] * 12),
            (["Water"], ["H2O"] * 12),
            (["ammonia"], ["NH3"] * 5),
            (["2"], ["H2O"] * 12),
            (["NH3"], ["NH3"] * 5),
            (["water", "Sodium"], ["H2O"] * 12 + ["Sodium"] * 3),
            ([], ["H2O"] * 12 + ["NH3"] * 5 + ["Sodium"] * 3),
        ],
    )
    def test_the_selected_molecules(self, mixed, monkeypatch, types, names):
        com = self.run(monkeypatch, *types)

        assert [atom.name for atom in com.atoms] == names

    def test_the_water_centers_are_those_of_all_molecules(self, mixed, monkeypatch):
        everything = self.run(monkeypatch, output="all.xyz")
        waters = self.run(monkeypatch, "water")

        # the molecules are in the order of the file: the waters come first
        assert np.allclose(waters.pos, everything.pos[:12], atol=1e-5)

    def test_water_without_a_water_type_line_raises(self, mixed, monkeypatch, tmp_path):
        (tmp_path / "moldescriptor.dat").write_text(
            self.MOLDESCRIPTOR.replace("WATER_TYPE 2\n", ""), encoding="utf-8"
        )

        with pytest.raises(PQValueError, match="has no WATER_TYPE line"):
            self.run(monkeypatch, "water")
