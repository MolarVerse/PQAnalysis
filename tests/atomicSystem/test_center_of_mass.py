"""
Tests for the center_of_mass_groups function.
"""

import numpy as np
import pytest

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system._center_of_mass import (
    center_of_mass_groups,
    center_of_mass_properties,
    group_mass_weighted_mean,
    group_sums,
)
from PQAnalysis.atomic_system.exceptions import (
    AtomicSystemError,
    AtomicSystemPositionsError,
)
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.topology import Topology

from . import pytestmark  # pylint: disable=unused-import


def reference_center_of_mass(positions, masses, cell):
    """The center of mass of one group as computed before the function."""

    relative = cell.image(positions - positions[0]) + positions[0]
    average = np.sum(relative * masses[:, None], axis=0) / np.sum(masses)

    return cell.image(average)


def difference_modulo_box(a, b, cell):
    return np.abs(cell.image(np.asarray(a) - np.asarray(b))).max()


class TestCenterOfMassGroups:

    def test_a_single_group_is_the_mass_weighted_average(self):
        positions = np.array([[0.0, 0, 0], [3.0, 0, 0]])
        masses = np.array([1.0, 3.0])

        com = center_of_mass_groups(positions, masses, np.zeros(2, int), Cell())

        assert com.shape == (1, 3)
        assert np.allclose(com, [[2.25, 0, 0]])

    def test_groups_of_different_sizes(self):
        positions = np.array(
            [[0.0, 0, 0], [2, 0, 0], [10, 0, 0], [10, 3, 0], [10, 0, 6]]
        )
        masses = np.array([1.0, 1.0, 2.0, 1.0, 3.0])
        group_ids = np.array([0, 0, 1, 1, 1])

        com = center_of_mass_groups(positions, masses, group_ids, Cell())

        assert np.allclose(com, [[1, 0, 0], [10, 0.5, 3.0]])

    def test_ids_do_not_need_to_be_sorted_or_contiguous(self):
        rng = np.random.default_rng(1)
        positions = rng.uniform(0, 5, (12, 3))
        masses = rng.uniform(1, 16, 12)
        group_ids = np.repeat([0, 1, 2, 3], 3)

        expected = center_of_mass_groups(positions, masses, group_ids, Cell())

        order = rng.permutation(12)
        shuffled = center_of_mass_groups(
            positions[order], masses[order], (group_ids * 7 + 3)[order], Cell()
        )

        # every group keeps its atoms, only the order of the atoms (and with
        # it the first atom of a group) changes; the rows follow the ids
        assert np.allclose(shuffled, expected)

    def test_groups_are_returned_in_ascending_order_of_their_ids(self):
        positions = np.array([[1.0, 0, 0], [2, 0, 0], [3, 0, 0]])
        masses = np.ones(3)

        com = center_of_mass_groups(
            positions, masses, np.array([9, 2, 5]), Cell()
        )

        assert np.allclose(com[:, 0], [2.0, 3.0, 1.0])

    @pytest.mark.parametrize(
        "cell",
        [
            Cell(10, 12, 9),
            Cell(10, 10, 10, 80, 95, 100),
            Cell(),
        ],
        ids=["orthorhombic", "triclinic", "vacuum"],
    )
    def test_agrees_with_the_formula_for_every_group(self, cell):
        rng = np.random.default_rng(7)
        n_groups, size = 20, 4
        masses = rng.uniform(1, 16, n_groups * size)
        positions = rng.uniform(-3, 3, (n_groups, 1, 3)) + rng.uniform(
            -1, 1, (n_groups, size, 3)
        )
        positions = cell.image(positions.reshape(-1, 3))

        com = center_of_mass_groups(
            positions, masses, np.repeat(np.arange(n_groups), size), cell
        )

        for index in range(n_groups):
            atoms = slice(index * size, (index + 1) * size)
            expected = reference_center_of_mass(
                positions[atoms], masses[atoms], cell
            )
            assert difference_modulo_box(com[index], expected, cell) < 1e-12

    def test_a_group_split_by_the_boundary_gets_the_center_of_the_group(self):
        cell = Cell(10, 10, 10)
        positions = np.array([[4.8, 0, 0], [-4.8, 0, 0]])  # 0.4 A apart
        masses = np.array([1.0, 1.0])

        com = center_of_mass_groups(positions, masses, np.zeros(2, int), cell)

        assert difference_modulo_box(com[0], [5.0, 0, 0], cell) < 1e-12

    def test_a_group_split_by_a_triclinic_lattice_vector(self):
        cell = Cell(10, 10, 10, 90, 90, 60)
        positions = np.array([[1.0, 1, 1], [1.5, 1, 1]])
        shifted = positions.copy()
        shifted[1] += cell.box_matrix @ np.array([1, 0, 0])
        masses = np.array([16.0, 1.0])

        com = center_of_mass_groups(positions, masses, np.zeros(2, int), cell)
        com_shifted = center_of_mass_groups(
            shifted, masses, np.zeros(2, int), cell
        )

        assert difference_modulo_box(com, com_shifted, cell) < 1e-12

    def test_the_result_lies_inside_the_cell(self):
        cell = Cell(10, 10, 10)
        positions = np.array([[30.0, -20, 14]])

        com = center_of_mass_groups(positions, np.ones(1), np.zeros(1, int), cell)

        assert np.all(np.abs(com) <= 5.0 + 1e-12)

    def test_integer_arrays_are_accepted(self):
        com = center_of_mass_groups(
            np.array([[0, 0, 0], [2, 0, 0]]),
            np.array([1, 1]),
            np.array([4, 4], dtype=np.int32),
            Cell(),
        )

        assert np.allclose(com, [[1, 0, 0]])

    def test_no_atoms_give_no_groups(self):
        com = center_of_mass_groups(
            np.zeros((0, 3)), np.zeros(0), np.zeros(0, dtype=int), Cell()
        )

        assert com.shape == (0, 3)

    def test_a_group_without_mass_raises(self):
        with pytest.raises(AtomicSystemError, match="total mass"):
            center_of_mass_groups(
                np.zeros((2, 3)), np.array([1.0, 0.0]), np.array([0, 1]), Cell()
            )

    def test_unequal_lengths_raise(self):
        with pytest.raises(AtomicSystemPositionsError, match="must be equal"):
            center_of_mass_groups(
                np.zeros((3, 3)), np.ones(2), np.zeros(3, dtype=int), Cell()
            )

    def test_positions_need_three_columns(self):
        with pytest.raises(AtomicSystemError, match="shape"):
            center_of_mass_groups(
                np.zeros((3, 2)), np.ones(3), np.zeros(3, dtype=int), Cell()
            )


class TestCheckExtent:

    """
    The optional check for groups that span half of the box or more.

    The positions of a group are imaged around its first atom, so the check
    sees groups with atoms on both sides of the first atom.
    """

    CELL = Cell(10, 10, 10)

    @staticmethod
    def run(positions, cell, group_ids=None, check_extent=True):
        positions = np.array(positions, float)

        return center_of_mass_groups(
            positions,
            np.ones(len(positions)),
            np.zeros(len(positions), int) if group_ids is None else group_ids,
            cell,
            check_extent=check_extent,
        )

    def test_the_example_of_a_molecule_larger_than_half_the_box(self):
        positions = [[0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]]

        with pytest.raises(AtomicSystemError, match="spans half of the box"):
            self.run(positions, self.CELL)

        # without the check the center of mass is silently wrong
        wrong = self.run(positions, self.CELL, check_extent=False)
        assert abs(wrong[0, 0] - 3.2) > 3

    def test_the_extent_limit_is_half_the_box(self):
        inside = [[0, 0, 0], [-2.4, 0, 0], [2.4, 0, 0]]  # spans 4.8 A
        outside = [[0, 0, 0], [-2.6, 0, 0], [2.6, 0, 0]]  # spans 5.2 A

        assert self.run(inside, self.CELL).shape == (1, 3)
        with pytest.raises(AtomicSystemError):
            self.run(outside, self.CELL)

    def test_every_direction_is_checked(self):
        for direction in np.eye(3):
            with pytest.raises(AtomicSystemError):
                self.run(
                    [[0, 0, 0], -2.6 * direction, 2.6 * direction], self.CELL
                )

    def test_the_check_uses_the_cell_vectors_of_a_triclinic_cell(self):
        cell = Cell(10, 10, 10, 90, 90, 60)
        a_vector = cell.box_matrix[:, 0]

        assert self.run(
            [[0, 0, 0], -0.24 * a_vector, 0.24 * a_vector], cell
        ).shape == (1, 3)
        with pytest.raises(AtomicSystemError):
            self.run([[0, 0, 0], -0.26 * a_vector, 0.26 * a_vector], cell)

    def test_a_vacuum_cell_has_no_limit(self):
        positions = [[0, 0, 0], [-1e4, 0, 0], [1e4, 0, 0]]

        assert self.run(positions, Cell()).shape == (1, 3)

    def test_only_the_large_groups_are_reported(self):
        positions = [
            [0, 0, 0], [-1, 0, 0], [1, 0, 0],
            [0, 0, 0], [-2.6, 0, 0], [2.6, 0, 0],
            [0, 0, 0], [0, -1, 0], [0, 1, 0],
        ]
        group_ids = np.repeat([5, 7, 9], 3)

        with pytest.raises(AtomicSystemError, match=r"group 7 .*\(1 of 3 groups"):
            self.run(positions, self.CELL, group_ids)

    def test_a_pair_farther_apart_than_half_the_box_cannot_be_detected(self):
        # 5.1 A apart cannot be told from the image 4.9 A apart (documented
        # limit of the check): the center of mass is taken of the image
        com = self.run([[0, 0, 0], [5.1, 0, 0]], self.CELL)

        assert np.isclose(com[0, 0], -2.45)

    def test_it_is_off_by_default(self):
        positions = np.array([[0.0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]])

        com = center_of_mass_groups(
            positions, np.ones(3), np.zeros(3, int), self.CELL
        )

        assert com.shape == (1, 3)

    def test_no_atoms_pass_the_check(self):
        com = self.run(np.zeros((0, 3)), self.CELL, np.zeros(0, int))

        assert com.shape == (0, 3)


class TestAtomicSystemUsesTheFunction:

    """
    The methods of AtomicSystem keep their results (they are wrappers of
    center_of_mass_groups).
    """

    @staticmethod
    def water_box(n_molecules=6, box=8.0):
        rng = np.random.default_rng(3)
        atoms = [Atom(name) for name in "OHH" * n_molecules]
        centers = rng.uniform(-box / 2, box / 2, (n_molecules, 1, 3))
        offsets = np.array([[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
        positions = (centers + offsets).reshape(-1, 3)
        cell = Cell(box, box, box)

        return AtomicSystem(atoms=atoms, pos=cell.image(positions), cell=cell)

    def test_center_of_mass_property(self):
        system = self.water_box(n_molecules=1)

        expected = reference_center_of_mass(
            system.pos, system.atomic_masses, system.cell
        )

        assert difference_modulo_box(system.center_of_mass, expected, system.cell) < 1e-12

    def test_compute_com_atomic_system_groups_consecutive_atoms(self):
        system = self.water_box()

        com_system = system.compute_com_atomic_system(group=3)

        assert com_system.n_atoms == 6
        assert [atom.name for atom in com_system.atoms] == ["OHH"] * 6
        for index in range(6):
            atoms = slice(3 * index, 3 * index + 3)
            expected = reference_center_of_mass(
                system.pos[atoms], system.atomic_masses[atoms], system.cell
            )
            assert difference_modulo_box(
                com_system.pos[index], expected, system.cell
            ) < 1e-12

    def test_compute_com_atomic_system_without_group_uses_all_atoms(self):
        system = self.water_box(n_molecules=1)

        com_system = system.compute_com_atomic_system()

        assert com_system.n_atoms == 1
        assert np.allclose(com_system.pos[0], system.center_of_mass)

    def test_compute_com_atomic_system_needs_a_divisor_of_the_atoms(self):
        with pytest.raises(AtomicSystemError, match="multiple of group"):
            self.water_box().compute_com_atomic_system(group=4)

    def test_center_of_mass_residues_of_molecules_split_by_the_boundary(self):
        reference_residues = [
            Residue(
                name="WAT",
                elements=[Element("O"), Element("H"), Element("H")],
                residue_id=1,
                total_charge=0,
                atom_types=np.array([0, 1, 1]),
                partial_charges=np.array([0, 0, 0]),
            )
        ]
        cell = Cell(10, 10, 10)
        topology = Topology(
            atoms=[Atom(name) for name in "OHH" * 2],
            residue_ids=np.array([1, 1, 1, 1, 1, 1]),
            reference_residues=reference_residues,
        )
        positions = np.array([
            [4.8, 0, 0], [-4.8, 0, 0], [4.9, 0.9, 0],  # split along x
            [0.0, 4.9, 0], [0.9, -4.8, 0], [-0.9, 4.9, 0],  # split along y
        ])
        system = AtomicSystem(topology=topology, pos=positions, cell=cell)

        com_system = system.center_of_mass_residues

        assert [atom.name for atom in com_system.atoms] == ["WAT", "WAT"]
        assert np.allclose(
            com_system.pos, system.compute_com_atomic_system(group=3).pos
        )
        for index in range(2):
            atoms = slice(3 * index, 3 * index + 3)
            expected = reference_center_of_mass(
                positions[atoms], system.atomic_masses[atoms], cell
            )
            assert difference_modulo_box(
                com_system.pos[index], expected, cell
            ) < 1e-12


class TestGroupReductions:

    GROUPS = np.array([4, 4, 1, 4, 1])
    MASSES = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    VALUES = np.arange(15, dtype=float).reshape(5, 3)

    def test_mass_weighted_mean(self):
        mean = group_mass_weighted_mean(self.VALUES, self.MASSES, self.GROUPS)

        # ascending ids: group 1 (atoms 2 and 4), then group 4 (0, 1, 3)
        assert mean.shape == (2, 3)
        assert np.allclose(
            mean[0],
            (3 * self.VALUES[2] + 5 * self.VALUES[4]) / 8,
        )
        assert np.allclose(
            mean[1],
            (self.VALUES[0] + 2 * self.VALUES[1] + 4 * self.VALUES[3]) / 7,
        )

    def test_sums_of_vectors_and_scalars(self):
        vectors = group_sums(self.VALUES, self.GROUPS)
        scalars = group_sums(self.MASSES, self.GROUPS)

        assert np.allclose(vectors[0], self.VALUES[[2, 4]].sum(0))
        assert np.allclose(vectors[1], self.VALUES[[0, 1, 3]].sum(0))
        assert np.allclose(scalars, [8.0, 7.0])

    def test_a_single_group_of_one_atom(self):
        assert np.allclose(
            group_mass_weighted_mean(
                self.VALUES[:1], self.MASSES[:1], np.array([0])
            ),
            self.VALUES[:1],
        )


class TestCenterOfMassProperties:

    CELL = Cell(10, 10, 10)

    def arrays(self):
        rng = np.random.default_rng(2)

        return {
            "positions": rng.uniform(0, 10, (6, 3)),
            "velocities": rng.normal(size=(6, 3)),
            "forces": rng.normal(size=(6, 3)),
            "charges": rng.normal(size=6),
        }

    def test_everything_that_is_given_is_returned(self):
        data = self.arrays()
        masses = np.arange(1.0, 7.0)
        groups = np.array([0, 0, 0, 1, 1, 1])

        result = center_of_mass_properties(masses, groups, self.CELL, **data)

        assert set(result) == {"pos", "vel", "forces", "charges"}
        assert np.allclose(
            result["pos"],
            center_of_mass_groups(data["positions"], masses, groups, self.CELL),
        )
        assert np.allclose(
            result["vel"],
            group_mass_weighted_mean(data["velocities"], masses, groups),
        )
        assert np.allclose(result["forces"], group_sums(data["forces"], groups))
        assert np.allclose(result["charges"], group_sums(data["charges"], groups))

    def test_none_and_empty_data_are_left_out(self):
        data = self.arrays()

        result = center_of_mass_properties(
            np.ones(6), np.zeros(6, int), self.CELL,
            positions=np.zeros((0, 3)),
            velocities=data["velocities"],
            forces=None,
            charges=np.zeros(0),
        )

        assert set(result) == {"vel"}

    def test_nothing_given_returns_nothing(self):
        assert center_of_mass_properties(np.ones(3), np.zeros(3, int), self.CELL) == {}

    def test_data_with_another_number_of_atoms_raises(self):
        with pytest.raises(AtomicSystemPositionsError, match="velocities of 3"):
            center_of_mass_properties(
                np.ones(6), np.zeros(6, int), self.CELL,
                velocities=np.zeros((3, 3)),
            )

    def test_a_group_without_mass_raises_even_without_positions(self):
        with pytest.raises(AtomicSystemError, match="total mass"):
            center_of_mass_properties(
                np.array([1.0, 0.0]), np.array([0, 1]), self.CELL,
                velocities=np.ones((2, 3)),
            )

    def test_the_extent_check_is_passed_on(self):
        with pytest.raises(AtomicSystemError, match="spans half of the box"):
            center_of_mass_properties(
                np.ones(3), np.zeros(3, int), self.CELL,
                positions=np.array([[0.0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]]),
                check_extent=True,
            )


class TestAtomicSystemWithVelocitiesForcesAndCharges:

    @staticmethod
    def system(positions=True):
        water = TestAtomicSystemUsesTheFunction.water_box(n_molecules=4)
        rng = np.random.default_rng(9)
        kwargs = {
            "vel": rng.normal(size=(12, 3)),
            "forces": rng.normal(size=(12, 3)),
            "charges": rng.normal(size=12),
        }
        reference_residues = [
            Residue(
                name="WAT",
                elements=[Element("O"), Element("H"), Element("H")],
                residue_id=1,
                total_charge=0,
                atom_types=np.array([0, 1, 1]),
                partial_charges=np.zeros(3),
            )
        ]
        topology = Topology(
            atoms=water.atoms,
            residue_ids=np.ones(12, dtype=int),
            reference_residues=reference_residues,
        )

        return AtomicSystem(
            topology=topology,
            pos=water.pos if positions else None,
            cell=water.cell,
            **kwargs,
        )

    @staticmethod
    def expected(system, molecule):
        atoms = slice(3 * molecule, 3 * molecule + 3)
        masses = system.atomic_masses[atoms]

        return (
            (system.vel[atoms] * masses[:, None]).sum(0) / masses.sum(),
            system.forces[atoms].sum(0),
            system.charges[atoms].sum(),
        )

    @pytest.mark.parametrize("positions", [True, False])
    def test_center_of_mass_residues(self, positions):
        system = self.system(positions)

        com = system.center_of_mass_residues

        assert com.n_atoms == 4
        assert com.has_pos == positions
        for molecule in range(4):
            vel, forces, charge = self.expected(system, molecule)
            assert np.allclose(com.vel[molecule], vel)
            assert np.allclose(com.forces[molecule], forces)
            assert com.charges[molecule] == pytest.approx(charge)

    def test_compute_com_atomic_system(self):
        system = self.system()

        com = system.compute_com_atomic_system(group=3)

        assert com.n_atoms == 4
        for molecule in range(4):
            vel, forces, charge = self.expected(system, molecule)
            assert np.allclose(com.vel[molecule], vel)
            assert np.allclose(com.forces[molecule], forces)
            assert com.charges[molecule] == pytest.approx(charge)

    def test_the_results_agree_with_each_other(self):
        system = self.system()

        by_residue = system.center_of_mass_residues
        by_group = system.compute_com_atomic_system(group=3)

        assert np.allclose(by_residue.vel, by_group.vel)
        assert np.allclose(by_residue.forces, by_group.forces)
        assert np.allclose(by_residue.charges, by_group.charges)
