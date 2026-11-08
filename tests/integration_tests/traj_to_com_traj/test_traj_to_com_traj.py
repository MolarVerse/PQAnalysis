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
from PQAnalysis.io import TrajectoryReader, TrajectoryWriter
from PQAnalysis.tools import traj_to_com_traj

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
