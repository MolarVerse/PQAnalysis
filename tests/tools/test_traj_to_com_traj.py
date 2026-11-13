import logging

import numpy as np
import pytest

from PQAnalysis.tools import traj_to_com_traj

from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system.exceptions import (
    AtomicSystemError,
    AtomicSystemMassError,
)
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.topology import Topology



def test_traj_to_com_traj():
    traj = Trajectory()

    assert traj_to_com_traj(traj) == Trajectory()

    atoms = [Atom(atom) for atom in ['C', 'C', 'C']]
    coordinates = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]

    frame = AtomicSystem(atoms=atoms, pos=np.array(coordinates))
    traj.append(frame)

    print(traj[0].n_atoms)

    traj_output = traj_to_com_traj(traj)

    assert np.allclose(traj_output[0].pos, [[1, 1, 1]])
    assert traj_output[0].combined_name == 'CCC'

    coordinates = [[0, 0, 1], [1, 1, 2], [2, 2, 3]]
    frame = AtomicSystem(atoms=atoms, pos=np.array(coordinates))

    traj.append(frame)

    traj_output = traj_to_com_traj(traj)

    assert np.allclose(traj_output[0].pos, [[1, 1, 1]])
    assert traj_output[0].combined_name == 'CCC'
    assert np.allclose(traj_output[1].pos, [[1, 1, 2]])
    assert traj_output[1].combined_name == 'CCC'

    # selection=None must mean "all atoms of each frame", not a
    # slice frozen from the first frame's atom count.
    traj_output = traj_to_com_traj(traj, selection=None)
    assert np.allclose(traj_output[0].pos, [[1, 1, 1]])
    assert np.allclose(traj_output[1].pos, [[1, 1, 2]])

    traj_output = traj_to_com_traj(traj, selection=slice(0, 2))
    assert np.allclose(traj_output[0].pos, [[0.5, 0.5, 0.5]])
    assert traj_output[0].combined_name == 'CC'
    assert np.allclose(traj_output[1].pos, [[0.5, 0.5, 1.5]])
    assert traj_output[1].combined_name == 'CC'

    traj = Trajectory()
    atoms = [Atom(atom) for atom in ['C', 'C', 'H', 'H']]
    coordinates1 = [[0, 0, 1], [1, 1, 2], [2, 2, 3], [3, 3, 4]]
    coordinates2 = [[0, 1, 1], [1, 2, 2], [2, 3, 3], [3, 4, 4]]
    frame1 = AtomicSystem(atoms=atoms, pos=np.array(coordinates1))
    frame2 = AtomicSystem(atoms=atoms, pos=np.array(coordinates2))
    traj.append(frame1)
    traj.append(frame2)

    traj_output = traj_to_com_traj(traj, group=2)

    assert np.allclose(traj_output[0].pos, [[0.5, 0.5, 1.5], [2.5, 2.5, 3.5]])
    assert traj_output[0].atoms == [
        Atom('CC',
        use_guess_element=False),
        Atom('HH',
        use_guess_element=False)
    ]
    assert traj_output[0].combined_name == 'CCHH'
    assert np.allclose(traj_output[1].pos, [[0.5, 1.5, 1.5], [2.5, 3.5, 3.5]])
    assert traj_output[1].atoms == [
        Atom('CC',
        use_guess_element=False),
        Atom('HH',
        use_guess_element=False)
    ]
    assert traj_output[1].combined_name == 'CCHH'



class TestByResidue:

    """
    Tests of the by_residue option: water, methane, water in a periodic box.
    Molecules of different sizes cannot be grouped by a number of atoms.
    """

    CELL = Cell(10, 10, 10)

    @staticmethod
    def topology():
        residues = [
            Residue(
                name="WAT",
                elements=[Element("O"), Element("H"), Element("H")],
                residue_id=1,
                total_charge=0,
                atom_types=np.array([0, 1, 1]),
                partial_charges=np.array([0, 0, 0]),
            ),
            Residue(
                name="MET",
                elements=[Element("C")] + [Element("H")] * 4,
                residue_id=2,
                total_charge=0,
                atom_types=np.array([0, 1, 1, 1, 1]),
                partial_charges=np.zeros(5),
            ),
        ]
        names = "OHH" + "CHHHH" + "OHH"

        return Topology(
            atoms=[Atom(name) for name in names],
            residue_ids=np.array([1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1]),
            reference_residues=residues,
        )

    def positions(self, shift=0.0):
        water = np.array([[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
        methane = np.array(
            [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
             [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]]
        )
        # the first water is split by the x boundary, the methane by y
        pos = np.vstack([
            water + [4.7, 1.0, 1.0],
            methane + [1.0, 4.8, -2.0],
            water + [-2.0, -3.0, 3.0],
        ]) + shift

        return self.CELL.image(pos)

    def trajectory(self, n_frames=2, with_topology=True):
        topology = self.topology()
        frames = [
            AtomicSystem(
                topology=topology if with_topology else None,
                atoms=None if with_topology else topology.atoms,
                pos=self.positions(shift=0.1 * frame),
                cell=self.CELL,
            ) for frame in range(n_frames)
        ]

        return Trajectory(frames)

    @staticmethod
    def expected(positions, atoms, cell):
        masses = np.array([atom.mass for atom in atoms])
        relative = cell.image(positions - positions[0])
        center = positions[0] + (
            relative * masses[:, None]
        ).sum(axis=0) / masses.sum()

        return center

    def assert_equal_modulo_box(self, actual, expected):
        assert np.abs(self.CELL.image(actual - expected)).max() < 1e-10

    def test_every_molecule_gets_its_center_of_mass(self):
        traj = self.trajectory()
        molecules = [slice(0, 3), slice(3, 8), slice(8, 11)]

        com_traj = traj_to_com_traj(traj, by_residue=True)

        assert len(com_traj) == 2
        for com_frame, frame in zip(com_traj, traj):
            assert com_frame.n_atoms == 3
            assert [atom.name for atom in com_frame.atoms] == [
                "WAT", "MET", "WAT"
            ]
            for index, atoms in enumerate(molecules):
                self.assert_equal_modulo_box(
                    com_frame.pos[index],
                    self.expected(frame.pos[atoms], frame.atoms[atoms], self.CELL),
                )

    def test_molecules_split_by_the_boundary_are_handled(self):
        frame = self.trajectory(1)[0]

        # without imaging the water at x = 4.7 would have its H at -4.34
        assert frame.pos[1, 0] < 0 < frame.pos[0, 0]

        com = traj_to_com_traj(Trajectory([frame]), by_residue=True)[0].pos[0]

        assert abs(com[0]) > 4.5  # next to the boundary, not near the middle

    def test_the_topology_argument_replaces_the_one_of_the_trajectory(self):
        with_topology = traj_to_com_traj(self.trajectory(), by_residue=True)
        bare = self.trajectory(with_topology=False)

        explicit = traj_to_com_traj(
            bare, by_residue=True, topology=self.topology()
        )

        for a, b in zip(with_topology, explicit):
            assert np.allclose(a.pos, b.pos)

    def test_the_selection_picks_the_atoms_of_every_molecule(self):
        traj = self.trajectory(1)
        frame = traj[0]
        hydrogens = np.array([1, 2, 4, 5, 6, 7, 9, 10])

        by_array = traj_to_com_traj(traj, selection=hydrogens, by_residue=True)
        by_atom = traj_to_com_traj(traj, selection=Atom("H"), by_residue=True)

        # named by the molecule, not by the selected atoms
        assert [atom.name for atom in by_array[0].atoms] == [
            "WAT", "MET", "WAT"
        ]
        assert np.allclose(by_array[0].pos, by_atom[0].pos)
        for index, atoms in enumerate([[1, 2], [4, 5, 6, 7], [9, 10]]):
            self.assert_equal_modulo_box(
                by_array[0].pos[index],
                self.expected(
                    frame.pos[atoms],
                    [frame.atoms[i] for i in atoms],
                    self.CELL,
                ),
            )

    def test_molecules_without_a_selected_atom_are_left_out(self):
        traj = self.trajectory(1)

        com_traj = traj_to_com_traj(traj, selection=Atom("C"), by_residue=True)

        assert com_traj[0].n_atoms == 1
        assert com_traj[0].atoms[0].name == "MET"
        self.assert_equal_modulo_box(com_traj[0].pos[0], traj[0].pos[3])

    def test_selection_by_slice_and_index(self):
        traj = self.trajectory(1)

        sliced = traj_to_com_traj(traj, selection=slice(0, 8), by_residue=True)
        single = traj_to_com_traj(traj, selection=3, by_residue=True)

        assert sliced[0].n_atoms == 2
        assert single[0].n_atoms == 1
        self.assert_equal_modulo_box(single[0].pos[0], traj[0].pos[3])

    def test_an_empty_trajectory_gives_an_empty_trajectory(self):
        assert traj_to_com_traj(Trajectory(), by_residue=True) == Trajectory()

    def test_group_and_by_residue_cannot_be_combined(self):
        with pytest.raises(PQValueError, match="cannot be combined"):
            traj_to_com_traj(self.trajectory(), group=3, by_residue=True)

    def test_a_topology_without_residues_raises(self):
        bare = self.trajectory(with_topology=False)

        with pytest.raises(PQValueError, match="no residues"):
            traj_to_com_traj(bare, by_residue=True)

    def test_a_topology_with_another_number_of_atoms_raises(self):
        traj = self.trajectory()[0:1]
        short = Trajectory([traj[0][0:8]])

        with pytest.raises(PQValueError, match="atoms per frame"):
            traj_to_com_traj(short, by_residue=True, topology=self.topology())


class TestByResidueValidation:

    """
    Wrong results that must fail loudly: molecules larger than half of the
    box, placeholder masses and atoms without a mass.
    """

    @staticmethod
    def system(
        elements, positions, cell, residue_name="MOL", atom_names=None
    ):
        residue = Residue(
            name=residue_name,
            elements=[Element(name) for name in elements],
            residue_id=1,
            total_charge=0,
            atom_types=np.arange(len(elements)),
            partial_charges=np.zeros(len(elements)),
        )
        topology = Topology(
            atoms=[
                Atom(name) if atom_names is None else Atom(
                    name, use_guess_element=False
                ) for name in (elements if atom_names is None else atom_names)
            ],
            residue_ids=np.ones(len(elements), dtype=int),
            reference_residues=[residue],
        )

        return AtomicSystem(
            topology=topology, pos=np.array(positions, float), cell=cell
        )

    def test_a_molecule_larger_than_half_the_box_raises(self):
        # 6.4 A long molecule in a 10 A box: its atoms cannot be told apart
        # from their periodic images (without the check the center of mass
        # is -0.13 A instead of 3.2 A)
        frame = self.system("CCC", [[0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]], Cell(10, 10, 10))

        with pytest.raises(AtomicSystemError, match="spans half of the box"):
            traj_to_com_traj(Trajectory([frame]), by_residue=True)

    def test_a_molecule_smaller_than_half_the_box_is_accepted(self):
        frame = self.system("CCC", [[0, 0, 0], [2.4, 0, 0], [4.8, 0, 0]], Cell(10, 10, 10))

        com_traj = traj_to_com_traj(Trajectory([frame]), by_residue=True)

        assert np.allclose(com_traj[0].pos, [[2.4, 0, 0]])

    def test_a_large_molecule_is_fine_in_a_vacuum_cell(self):
        frame = self.system("CCC", [[0, 0, 0], [30, 0, 0], [60, 0, 0]], Cell())

        com_traj = traj_to_com_traj(Trajectory([frame]), by_residue=True)

        assert np.allclose(com_traj[0].pos, [[30, 0, 0]])

    def test_group_mode_does_not_check_the_extent(self):
        frame = self.system("CCC", [[0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]], Cell(10, 10, 10))

        com_traj = traj_to_com_traj(Trajectory([frame]), group=3)

        assert com_traj[0].n_atoms == 1

    def test_atoms_with_a_placeholder_mass_log_a_warning(self, caplog):
        frame = self.system("OHX", [[0, 0, 0], [1, 0, 0], [0, 1, 0]], Cell(10, 10, 10))

        with caplog.at_level(logging.WARNING):
            traj_to_com_traj(Trajectory([frame]), by_residue=True)

        assert "placeholder mass (x)" in caplog.text.lower()

    def test_the_placeholder_warning_is_not_logged_when_they_are_excluded(
        self, caplog
    ):
        frame = self.system("OHX", [[0, 0, 0], [1, 0, 0], [0, 1, 0]], Cell(10, 10, 10))

        with caplog.at_level(logging.WARNING):
            traj_to_com_traj(
                Trajectory([frame]), selection=slice(0, 2), by_residue=True
            )

        assert "placeholder" not in caplog.text

    def test_atoms_without_a_mass_raise(self):
        frame = self.system(
            "OHH",
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            Cell(10, 10, 10),
            atom_names=["O1", "H1", "H2"],
        )

        with pytest.raises(AtomicSystemMassError):
            traj_to_com_traj(Trajectory([frame]), by_residue=True)

    def test_the_result_carries_no_residues(self):
        frame = self.system("OHH", [[0, 0, 0], [1, 0, 0], [0, 1, 0]], Cell(10, 10, 10), "WAT")

        com_traj = traj_to_com_traj(Trajectory([frame]), by_residue=True)

        assert com_traj[0].atoms[0].name == "WAT"
        assert not com_traj[0].topology.residues


class TestTraj2COM:

    """Tests of the file based traj2com function."""

    @staticmethod
    def write_xyz(filename, frames=2):
        with open(filename, "w", encoding="utf-8") as file:
            for frame in range(frames):
                shift = 0.1 * frame
                print("4 10.0 10.0 10.0 90.0 90.0 90.0\n", file=file)
                print(f"C {shift} 0.0 0.0", file=file)
                print(f"C {2.0 + shift} 0.0 0.0", file=file)
                print(f"H {shift} 2.0 0.0", file=file)
                print(f"H {shift} 4.0 0.0", file=file)

    @staticmethod
    def read_positions(filename):
        from PQAnalysis.io import TrajectoryReader

        return [frame.pos for frame in TrajectoryReader(filename).read()]

    @pytest.mark.usefixtures("tmpdir")
    def test_groups_of_atoms(self):
        from PQAnalysis.tools import traj2com

        self.write_xyz("traj.xyz")

        traj2com(["traj.xyz"], output="com.xyz", group=2)

        first, second = self.read_positions("com.xyz")
        assert np.allclose(first, [[1.0, 0, 0], [0, 3.0, 0]], atol=1e-6)
        assert np.allclose(second, [[1.1, 0, 0], [0.1, 3.0, 0]], atol=1e-6)

    @pytest.mark.usefixtures("tmpdir")
    def test_selection_string(self):
        from PQAnalysis.tools import traj2com

        self.write_xyz("traj.xyz", frames=1)

        traj2com(["traj.xyz"], output="com.xyz", selection="C")

        (positions,) = self.read_positions("com.xyz")
        assert np.allclose(positions, [[1.0, 0, 0]], atol=1e-6)

    @pytest.mark.usefixtures("tmpdir")
    def test_stdout_without_output_file(self, capsys):
        from PQAnalysis.tools import traj2com

        self.write_xyz("traj.xyz", frames=1)

        traj2com(["traj.xyz"], selection="H")

        assert capsys.readouterr().out.splitlines()[0].split()[0] == "1"

    @pytest.mark.usefixtures("tmpdir")
    def test_an_empty_selection_raises(self):
        from PQAnalysis.tools import traj2com

        self.write_xyz("traj.xyz", frames=1)

        with pytest.raises(PQValueError, match="does not select any atoms"):
            traj2com(["traj.xyz"], output="com.xyz", selection="N")

    @pytest.mark.usefixtures("tmpdir")
    def test_a_moldescriptor_file_needs_a_restart_file(self):
        from PQAnalysis.tools import traj2com

        self.write_xyz("traj.xyz", frames=1)

        with pytest.raises(PQValueError, match="needs a restart file"):
            traj2com(
                ["traj.xyz"], output="com.xyz", moldescriptor_file="mol.dat"
            )

    @pytest.mark.usefixtures("tmpdir")
    def test_by_residue_without_molecules_raises(self):
        from PQAnalysis.tools import traj2com

        self.write_xyz("traj.xyz", frames=1)

        with pytest.raises(PQValueError, match="no residues"):
            traj2com(["traj.xyz"], output="com.xyz", by_residue=True)


def moving_molecules(cell, n_frames=8, seed=3):
    """
    Water, methane and ion that move with a constant velocity, so that they
    cross the boundary of the box several times, with every atom wrapped on its
    own. Returns the trajectory, the topology and the true (unwrapped) centers
    of mass of every frame.
    """

    topology = TestByResidue.topology()
    masses = np.array([atom.mass for atom in topology.atoms])
    rng = np.random.default_rng(seed)
    start = TestByResidue().positions()
    # the largest step is 0.4 of the shortest box length, below one half
    steps = rng.uniform(-0.4, 0.4, (3, 3)) @ cell.box_matrix.T
    molecule = topology.residue_numbers
    names = np.unique(molecule)

    frames, centers = [], []
    for frame in range(n_frames):
        unwrapped = start.copy()
        for index, molecule_id in enumerate(names):
            unwrapped[molecule == molecule_id] += frame * steps[index]
        frames.append(
            AtomicSystem(topology=topology, pos=cell.image(unwrapped), cell=cell)
        )
        centers.append(
            np.array([
                (unwrapped[molecule == i] * masses[molecule == i, None]).sum(0)
                / masses[molecule == i].sum() for i in names
            ])
        )

    return Trajectory(frames), topology, np.array(centers)


class TestContinuous:

    CELLS = {
        "orthorhombic": Cell(10, 12, 9),
        "triclinic": Cell(10, 10, 10, 80, 95, 100),
    }

    @pytest.mark.parametrize("cell", CELLS.values(), ids=CELLS.keys())
    def test_the_centers_follow_their_way_over_the_boundary(self, cell):
        traj, _, true_centers = moving_molecules(cell)

        folded = traj_to_com_traj(traj, by_residue=True)
        continuous = traj_to_com_traj(traj, by_residue=True, continuous=True)

        folded_pos = np.array([frame.pos for frame in folded])
        pos = np.array([frame.pos for frame in continuous])

        # the displacement from the first frame is the true displacement
        assert np.allclose(
            pos - pos[0], true_centers - true_centers[0], atol=1e-10
        )
        # the first frame is folded, like without the option
        assert np.allclose(pos[0], folded_pos[0])
        # the folded centers jump over the boundary, the continuous ones not
        box = np.linalg.norm(cell.box_matrix, axis=0).min()
        assert np.abs(np.diff(folded_pos, axis=0)).max() > 0.5 * box
        assert np.abs(np.diff(pos, axis=0)).max() < 0.5 * box
        # and the folded centers are the continuous ones folded back
        assert np.allclose(
            cell.image(pos - folded_pos), 0.0, atol=1e-10
        )

    def test_it_gives_the_mean_square_displacement_of_the_true_motion(self):
        cell = self.CELLS["orthorhombic"]
        traj, _, true_centers = moving_molecules(cell)

        continuous = traj_to_com_traj(traj, by_residue=True, continuous=True)

        pos = np.array([frame.pos for frame in continuous])
        msd = ((pos - pos[0])**2).sum(-1).mean(-1)
        expected = ((true_centers - true_centers[0])**2).sum(-1).mean(-1)
        assert np.allclose(msd, expected)
        assert msd[-1] > 100.0  # many box lengths squared

    def test_groups_of_atoms_can_be_continuous_as_well(self):
        cell = self.CELLS["orthorhombic"]
        traj, _, true_centers = moving_molecules(cell)

        # the first three atoms are the water
        continuous = traj_to_com_traj(
            traj, selection=slice(0, 3), group=3, continuous=True
        )

        pos = np.array([frame.pos[0] for frame in continuous])
        assert np.allclose(
            pos - pos[0], true_centers[:, 0] - true_centers[0, 0], atol=1e-10
        )

    def test_a_cell_that_changes_in_size_is_followed(self):
        # NPT: the box grows by 5 % per frame while the centers travel 0.4 of
        # a box per frame, so after a few frames they are several boxes away
        # from the start and the growth alone moves them by more than half a
        # box. Following the Cartesian positions would pick the wrong image;
        # the fractional coordinates change smoothly.
        topology = TestByResidue.topology()
        start = np.array([0.45, 0.2, 0.3])
        speeds = np.array([0.4, -0.4, 0.3])
        molecule = topology.residue_numbers
        names = np.unique(molecule)
        frames, true_fractional = [], []
        for frame in range(8):
            box = 10.0 * (1.0 + 0.05 * frame)
            cell = Cell(box, box, box)
            pos = np.empty((11, 3))
            for index, molecule_id in enumerate(names):
                fractional = start + speeds[index] * frame
                pos[molecule == molecule_id] = fractional * box + 1.5 * index
            frames.append(
                AtomicSystem(topology=topology, pos=cell.image(pos), cell=cell)
            )
            true_fractional.append(
                np.array([
                    start + speeds[i] * frame + 1.5 * i / box for i in range(3)
                ])
            )

        continuous = traj_to_com_traj(
            Trajectory(frames), by_residue=True, continuous=True
        )
        fractional = np.array([
            frame.pos @ frame.cell.inverse_box_matrix.T for frame in continuous
        ])

        # the continuous fractional coordinates differ from the true ones by
        # the same lattice vector in every frame
        offset = fractional - np.array(true_fractional)
        assert np.allclose(offset, np.round(offset), atol=1e-9)
        assert np.allclose(offset, offset[0], atol=1e-9)
        assert np.abs(np.diff(fractional, axis=0)).max() < 0.5

    def test_a_vacuum_cell_is_left_alone(self):
        topology = TestByResidue.topology()
        frames = [
            AtomicSystem(
                topology=topology,
                pos=TestByResidue().positions() + 30.0 * frame,
                cell=Cell(),
            ) for frame in range(3)
        ]

        continuous = traj_to_com_traj(
            Trajectory(frames), by_residue=True, continuous=True
        )
        plain = traj_to_com_traj(Trajectory(frames), by_residue=True)

        for a, b in zip(continuous, plain):
            assert np.array_equal(a.pos, b.pos)

    def test_a_single_frame_is_unchanged(self):
        traj, _, _ = moving_molecules(self.CELLS["orthorhombic"], n_frames=1)

        assert np.allclose(
            traj_to_com_traj(traj, by_residue=True, continuous=True)[0].pos,
            traj_to_com_traj(traj, by_residue=True)[0].pos,
        )

    def test_the_default_is_unchanged(self):
        traj, _, _ = moving_molecules(self.CELLS["orthorhombic"])

        default = traj_to_com_traj(traj, by_residue=True)
        explicit = traj_to_com_traj(traj, by_residue=True, continuous=False)

        for a, b in zip(default, explicit):
            assert np.array_equal(a.pos, b.pos)

    def test_positions_are_needed(self):
        topology = TestByResidue.topology()
        frame = AtomicSystem(
            topology=topology, vel=np.ones((11, 3)), cell=Cell(10, 10, 10)
        )

        with pytest.raises(PQValueError, match="needs the positions"):
            traj_to_com_traj(Trajectory([frame]), by_residue=True, continuous=True)


def with_properties(seed=5, positions=True, **present):
    """
    The molecules of TestByResidue with random velocities, forces and charges
    (as far as requested) and the independent per molecule references.
    """

    rng = np.random.default_rng(seed)
    topology = TestByResidue.topology()
    masses = np.array([atom.mass for atom in topology.atoms])
    molecule = topology.residue_numbers
    cell = TestByResidue.CELL
    data = {}
    if present.get("vel"):
        data["vel"] = rng.normal(size=(11, 3))
    if present.get("forces"):
        data["forces"] = rng.normal(size=(11, 3))
    if present.get("charges"):
        data["charges"] = rng.normal(size=11)
    if positions:
        data["pos"] = TestByResidue().positions()

    frame = AtomicSystem(topology=topology, cell=cell, **data)

    expected = {"vel": [], "forces": [], "charges": []}
    for molecule_id in np.unique(molecule):
        atoms = molecule == molecule_id
        if "vel" in data:
            expected["vel"].append(
                (data["vel"][atoms] * masses[atoms, None]).sum(0) / masses[atoms].sum()
            )
        if "forces" in data:
            expected["forces"].append(data["forces"][atoms].sum(0))
        if "charges" in data:
            expected["charges"].append(data["charges"][atoms].sum())

    return frame, {key: np.array(value) for key, value in expected.items()}


class TestVelocitiesForcesAndCharges:

    def test_by_residue_combines_all_the_data(self):
        frame, expected = with_properties(vel=True, forces=True, charges=True)

        com = traj_to_com_traj(Trajectory([frame]), by_residue=True)[0]

        assert np.allclose(com.vel, expected["vel"])
        assert np.allclose(com.forces, expected["forces"])
        assert np.allclose(com.charges, expected["charges"])
        # the positions do not depend on the other data
        plain, _ = with_properties()
        assert np.allclose(
            com.pos, traj_to_com_traj(Trajectory([plain]), by_residue=True)[0].pos
        )

    @pytest.mark.parametrize("key", ["vel", "forces", "charges"])
    def test_data_that_is_not_there_is_not_made_up(self, key):
        frame, _ = with_properties(**{key: True})

        com = traj_to_com_traj(Trajectory([frame]), by_residue=True)[0]

        assert {
            "vel": com.has_vel, "forces": com.has_forces, "charges": com.has_charges
        } == {
            name: name == key for name in ("vel", "forces", "charges")
        }

    @pytest.mark.parametrize("key", ["vel", "forces", "charges"])
    def test_frames_without_positions(self, key):
        # a vel, force or charge trajectory has no positions
        frame, expected = with_properties(positions=False, **{key: True})

        com = traj_to_com_traj(Trajectory([frame]), by_residue=True)[0]

        assert not com.has_pos
        assert com.n_atoms == 3
        assert [atom.name for atom in com.atoms] == ["WAT", "MET", "WAT"]
        assert np.allclose(getattr(com, key), expected[key])

    def test_a_selection_picks_the_atoms_of_the_data_too(self):
        frame, _ = with_properties(vel=True, charges=True)
        hydrogens = np.array([1, 2, 4, 5, 6, 7, 9, 10])
        masses = np.array([atom.mass for atom in frame.atoms])
        molecule = frame.topology.residue_numbers

        com = traj_to_com_traj(
            Trajectory([frame]), selection=hydrogens, by_residue=True
        )[0]

        for index, molecule_id in enumerate([0, 1, 2]):
            atoms = hydrogens[molecule[hydrogens] == molecule_id]
            assert np.allclose(
                com.vel[index],
                (frame.vel[atoms] * masses[atoms, None]).sum(0) / masses[atoms].sum(),
            )
            assert com.charges[index] == pytest.approx(frame.charges[atoms].sum())

    def test_groups_of_atoms_combine_the_data(self):
        frame, _ = with_properties(vel=True, forces=True, charges=True)
        masses = np.array([atom.mass for atom in frame.atoms])

        com = traj_to_com_traj(
            Trajectory([frame]), selection=slice(0, 8), group=4
        )[0]

        for index in range(2):
            atoms = slice(4 * index, 4 * index + 4)
            assert np.allclose(
                com.vel[index],
                (frame.vel[atoms] * masses[atoms, None]).sum(0) / masses[atoms].sum(),
            )
            assert np.allclose(com.forces[index], frame.forces[atoms].sum(0))
            assert com.charges[index] == pytest.approx(frame.charges[atoms].sum())

    def test_the_center_of_mass_of_all_atoms_has_the_total_momentum(self):
        frame, _ = with_properties(vel=True, forces=True, charges=True)
        masses = np.array([atom.mass for atom in frame.atoms])

        com = traj_to_com_traj(Trajectory([frame]))[0]

        assert np.allclose(
            com.vel[0] * masses.sum(), (frame.vel * masses[:, None]).sum(0)
        )
        assert np.allclose(com.forces[0], frame.forces.sum(0))

    def test_atoms_without_a_mass_raise_also_without_positions(self):
        topology = TestByResidue.topology()
        topology.atoms[0].element = Element()  # no mass
        frame = AtomicSystem(
            topology=topology, vel=np.ones((11, 3)), cell=Cell(10, 10, 10)
        )

        with pytest.raises(AtomicSystemMassError):
            traj_to_com_traj(Trajectory([frame]), by_residue=True)


class TestMoleculeTypes:

    """Only the molecules of some types, by name or number (WAT, MET, WAT)."""

    @staticmethod
    def traj():
        return TestByResidue().trajectory(1)

    @pytest.mark.parametrize(
        "types, names",
        [(["WAT"], ["WAT", "WAT"]), ([2], ["MET"]), (["WAT", 2], ["WAT", "MET", "WAT"]), ([1], ["WAT", "WAT"])],
    )
    def test_only_the_molecules_of_the_types(self, types, names):
        com = traj_to_com_traj(self.traj(), by_residue=True, molecule_types=types)[0]

        assert [atom.name for atom in com.atoms] == names

    def test_the_centers_equal_those_of_all_molecules(self):
        traj = self.traj()
        everything = traj_to_com_traj(traj, by_residue=True)[0]
        waters = traj_to_com_traj(traj, by_residue=True, molecule_types=["WAT"])[0]
        methane = traj_to_com_traj(traj, by_residue=True, molecule_types=["MET"])[0]

        assert np.array_equal(waters.pos, everything.pos[[0, 2]])
        assert np.array_equal(methane.pos, everything.pos[[1]])

    def test_the_selection_is_combined_with_the_types(self):
        traj = self.traj()
        hydrogens = np.array([1, 2, 4, 5, 6, 7, 9, 10])

        com = traj_to_com_traj(
            traj, selection=hydrogens, by_residue=True, molecule_types=["MET"]
        )[0]
        reference = traj_to_com_traj(
            traj, selection=hydrogens[2:6], by_residue=True
        )[0]

        assert com.n_atoms == 1
        assert np.array_equal(com.pos, reference.pos)

    def test_velocities_of_the_waters_only(self):
        frame, _ = with_properties(vel=True)
        masses = np.array([atom.mass for atom in frame.atoms])
        waters = np.r_[0:3, 8:11]

        com = traj_to_com_traj(
            Trajectory([frame]), by_residue=True, molecule_types=["WAT"]
        )[0]

        for row, atoms in enumerate((waters[:3], waters[3:])):
            assert np.allclose(
                com.vel[row],
                (frame.vel[atoms] * masses[atoms, None]).sum(0) / masses[atoms].sum(),
            )

    def test_an_unknown_type_raises_and_lists_the_known_ones(self):
        with pytest.raises(PQValueError, match=r"\['ION'\] are not in the topology.*\['MET', 'WAT'\].*\[1, 2\]"):
            traj_to_com_traj(self.traj(), by_residue=True, molecule_types=["ION"])

    def test_types_need_by_residue(self):
        with pytest.raises(PQValueError, match="needs by_residue"):
            traj_to_com_traj(self.traj(), group=3, molecule_types=["WAT"])


class TestTypeKeywords:

    TEXT = (
        "WATER_TYPE 2\nAMMONIA_TYPE 3\n"
        "Na 1 1.0\nNa 0 1.0\n"
        "H2O 3 0.0\nO 0 -0.8\nH 1 0.4\nH 1 0.4\n"
        "NH3 4 0.0\nN 0 -0.8\nH 1 0.27\nH 1 0.27\nH 1 0.27\n"
    )

    @staticmethod
    def resolve(types, text=None, tmp_path=None):
        from PQAnalysis.tools.traj_to_com_traj import _resolve_type_keywords

        path = None
        if text is not None:
            path = tmp_path / "moldescriptor.dat"
            path.write_text(text, encoding="utf-8")
            path = str(path)

        return _resolve_type_keywords(types, path)

    def test_keywords_become_the_numbers_of_the_file(self, tmp_path):
        assert self.resolve(["Water", "AMMONIA", "Na", 1], self.TEXT, tmp_path) == [2, 3, "Na", 1]

    def test_nothing_to_resolve(self):
        assert self.resolve(None) is None
        assert self.resolve(["WAT", 2]) == ["WAT", 2]

    def test_a_keyword_needs_the_moldescriptor_file(self):
        with pytest.raises(PQValueError, match="needs the moldescriptor file"):
            self.resolve(["water"])

    def test_a_keyword_needs_its_line(self, tmp_path):
        with pytest.raises(PQValueError, match="has no AMMONIA_TYPE line"):
            self.resolve(["ammonia"], "WATER_TYPE 1\nH2O 3 0.0\nO 0 -0.8\nH 1 0.4\nH 1 0.4\n", tmp_path)
