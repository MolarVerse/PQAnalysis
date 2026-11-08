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
