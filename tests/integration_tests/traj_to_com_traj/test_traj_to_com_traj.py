"""
Integration test of traj_to_com_traj on the UMCM-9 example: read a real
trajectory file, compute the center of mass trajectory for the atoms of
the whole selection or for groups of N consecutive atoms, write it and read
it back.

The expected centers of mass are computed independently from the raw text of
the xyz file. They are compared modulo the box, because the tool folds the
center of mass into [-L/2, L/2) while the reference does not need to.
"""

import numpy as np
import pytest

from PQAnalysis.atomic_system.exceptions import AtomicSystemError
from PQAnalysis.core import Atom
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.io import (
    TrajectoryReader,
    TrajectoryWriter,
    read_restart_file,
)
from PQAnalysis.tools import traj_to_com_traj
from PQAnalysis.traj import Trajectory

pytestmark = pytest.mark.integration

FILENAME = "umcm-9-md-01.xyz"
N_FRAMES = 3  # the tool takes about 0.2 s per frame for 8-atom groups

# The first atom of the file is a dummy "X" of mass 999. It is excluded so that
# the remaining 808 = 8 * 101 atoms can be grouped.
SELECTION = slice(1, None)


def read_raw_xyz(filename, n_frames):
    """Reads names, positions and box of every frame from the xyz text."""

    lines = open(filename, encoding="utf-8").read().splitlines()
    frames = []
    position = 0

    for _ in range(n_frames):
        header = lines[position].split()
        n_atoms = int(header[0])
        box = np.array([float(value) for value in header[1:4]])
        atoms = [line.split() for line in lines[position + 2:position + 2 + n_atoms]]
        frames.append((
            [atom[0] for atom in atoms],
            np.array([[float(value) for value in atom[1:4]] for atom in atoms]),
            box,
        ))
        position += 2 + n_atoms

    return frames


def read_raw_restart(filename):
    """Reads names, positions and box of a PQ restart file from its text."""

    names, positions, box = [], [], None

    for line in open(filename, encoding="utf-8").read().splitlines():
        tokens = line.split()

        if not tokens or tokens[0] == "Step":
            continue
        if tokens[0] == "Box":
            box = np.array([float(value) for value in tokens[1:4]])
            continue

        names.append(tokens[0])
        positions.append([float(value) for value in tokens[3:6]])

    return names, np.array(positions), box


def reference_coms(names, positions, box, group):
    """Mass weighted centers of mass of consecutive groups of atoms."""

    masses = np.array([Atom(name).mass for name in names])
    group = len(names) if group is None else group
    coms = []

    for start in range(0, len(names), group):
        pos = positions[start:start + group]
        mass = masses[start:start + group]

        relative = pos - pos[0]
        relative -= box * np.round(relative / box)
        coms.append(pos[0] + (relative * mass[:, None]).sum(axis=0) / mass.sum())

    return np.array(coms)


def assert_equal_modulo_box(actual, expected, box, atol=1e-4):
    difference = actual - expected
    difference -= box * np.round(difference / box)

    assert np.abs(difference).max() < atol


@pytest.mark.parametrize("example_dir", ["traj2comtraj"], indirect=False)
class TestTrajToComTraj:

    @pytest.mark.parametrize("group", [None, 101, 8, 4])
    def test_groups_of_consecutive_atoms(self, test_integration_folder, group):
        traj = TrajectoryReader(FILENAME).read()[:N_FRAMES]
        raw_frames = read_raw_xyz(FILENAME, N_FRAMES)

        com_traj = traj_to_com_traj(traj, selection=SELECTION, group=group)

        assert len(com_traj) == N_FRAMES

        for frame, (names, positions, box) in zip(com_traj, raw_frames):
            names, positions = names[1:], positions[1:]
            size = len(names) if group is None else group

            assert frame.n_atoms == len(names) // size
            assert [atom.name for atom in frame.atoms] == [
                "".join(names[start:start + size])
                for start in range(0, len(names), size)
            ]
            assert np.allclose(frame.cell.box_lengths, box)
            assert_equal_modulo_box(
                frame.pos, reference_coms(names, positions, box, group), box
            )

    def test_written_and_read_back_trajectory(self, test_integration_folder):
        traj = TrajectoryReader(FILENAME).read()[:N_FRAMES]
        com_traj = traj_to_com_traj(traj, selection=SELECTION, group=8)

        TrajectoryWriter("com.xyz").write(com_traj)
        read_back = TrajectoryReader("com.xyz").read()

        assert len(read_back) == N_FRAMES

        for expected, actual in zip(com_traj, read_back):
            assert actual.n_atoms == expected.n_atoms == 101
            assert actual.combined_name == expected.combined_name
            assert actual.cell == expected.cell
            assert np.allclose(actual.pos, expected.pos, atol=1e-5)

    def test_group_must_divide_the_number_of_atoms(
        self, test_integration_folder
    ):
        traj = TrajectoryReader(FILENAME).read()[:1]

        with pytest.raises(AtomicSystemError, match="multiple of group"):
            traj_to_com_traj(traj, selection=SELECTION, group=3)

    def test_the_dummy_atom_pulls_the_center_of_mass_towards_itself(
        self, test_integration_folder
    ):
        # documents why the tests above exclude the dummy atom X (mass 999)
        traj = TrajectoryReader(FILENAME).read()[:1]
        names, positions, box = read_raw_xyz(FILENAME, 1)[0]

        with_dummy = traj_to_com_traj(traj, selection=slice(0, 9), group=9)
        without_dummy = traj_to_com_traj(traj, selection=slice(1, 9), group=8)

        assert_equal_modulo_box(
            with_dummy[0].pos,
            reference_coms(names[:9], positions[:9], box, 9),
            box,
        )

        def distance_to_dummy(com):
            difference = com - positions[0]
            difference -= box * np.round(difference / box)

            return np.linalg.norm(difference)

        # 5.7 A away without the dummy atom, 1.7 A with it
        assert distance_to_dummy(without_dummy[0].pos[0]) > 5.0
        assert distance_to_dummy(with_dummy[0].pos[0]) < 2.0


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestWaterBox:

    """
    2000 rigid TIP3P waters in a cubic box of 39.4491 A (the water box of the
    PQ example used for the O-O radial distribution function). 128 molecules
    are split by the periodic boundary.
    """

    N_MOLECULES = 2000

    @pytest.fixture
    def water_box(self, test_integration_folder):
        system = read_restart_file("input.rst", "moldescriptor.dat")
        raw = read_raw_restart("input.rst")

        return system, raw

    def test_the_example_has_molecules_split_by_the_boundary(self, water_box):
        _, (_, positions, _) = water_box

        molecules = positions.reshape(-1, 3, 3)
        longest_oh = np.maximum(
            np.linalg.norm(molecules[:, 1] - molecules[:, 0], axis=1),
            np.linalg.norm(molecules[:, 2] - molecules[:, 0], axis=1),
        )

        assert len(molecules) == self.N_MOLECULES
        assert (longest_oh > 2.0).sum() == 128

    def test_groups_of_three_atoms_are_the_residues_of_the_moldescriptor(
        self, water_box
    ):
        system, _ = water_box

        residues = system.topology.residue_atom_indices

        assert len(residues) == self.N_MOLECULES
        assert all(
            list(atoms) == [3 * index, 3 * index + 1, 3 * index + 2]
            for index, atoms in enumerate(residues)
        )

    def test_group_of_three_gives_the_center_of_mass_of_every_water(
        self, water_box
    ):
        system, (names, positions, box) = water_box

        com_traj = traj_to_com_traj(Trajectory([system]), group=3)
        frame = com_traj[0]

        assert frame.n_atoms == self.N_MOLECULES
        assert all(atom.name == "OHH" for atom in frame.atoms)
        assert np.allclose(frame.cell.box_lengths, box)
        assert_equal_modulo_box(
            frame.pos, reference_coms(names, positions, box, 3), box
        )

    def test_every_center_of_mass_lies_next_to_its_own_oxygen(
        self, water_box
    ):
        # independent of the reference above: the center of mass of a rigid
        # TIP3P water is about 0.065 A from the oxygen atom, also for the
        # molecules that are split by the periodic boundary
        system, (_, positions, box) = water_box

        frame = traj_to_com_traj(Trajectory([system]), group=3)[0]

        difference = frame.pos - positions[::3]
        difference -= box * np.round(difference / box)
        distance = np.linalg.norm(difference, axis=1)

        assert distance.min() > 0.05
        assert distance.max() < 0.08

    def test_written_and_read_back_trajectory(self, water_box):
        system, _ = water_box
        com_traj = traj_to_com_traj(Trajectory([system]), group=3)

        TrajectoryWriter("com.xyz").write(com_traj)
        read_back = TrajectoryReader("com.xyz").read()

        assert len(read_back) == 1
        assert read_back[0].n_atoms == self.N_MOLECULES
        assert read_back[0].cell == com_traj[0].cell
        assert np.allclose(read_back[0].pos, com_traj[0].pos, atol=1e-5)

    def test_by_residue_gives_the_same_as_groups_of_three(self, water_box):
        system, (names, positions, box) = water_box
        traj = Trajectory([system])

        by_residue = traj_to_com_traj(traj, by_residue=True)[0]
        by_group = traj_to_com_traj(traj, group=3)[0]

        assert by_residue.n_atoms == self.N_MOLECULES
        assert by_residue.atoms == by_group.atoms
        assert np.allclose(by_residue.pos, by_group.pos, atol=1e-12)
        assert_equal_modulo_box(
            by_residue.pos, reference_coms(names, positions, box, 3), box
        )

    def test_by_residue_for_an_xyz_trajectory_with_the_restart_topology(
        self, water_box
    ):
        # the usual workflow: a trajectory without molecule information (xyz)
        # and the topology from the restart file and the moldescriptor file
        system, _ = water_box
        TrajectoryWriter("water.xyz").write(Trajectory([system]))
        traj = TrajectoryReader("water.xyz").read()

        with pytest.raises(PQValueError, match="no residues"):
            traj_to_com_traj(traj, by_residue=True)

        com_traj = traj_to_com_traj(
            traj, by_residue=True, topology=system.topology
        )
        expected = traj_to_com_traj(Trajectory([system]), by_residue=True)

        assert com_traj[0].n_atoms == self.N_MOLECULES
        assert_equal_modulo_box(
            com_traj[0].pos, expected[0].pos, expected[0].cell.box_lengths, atol=1e-4
        )
