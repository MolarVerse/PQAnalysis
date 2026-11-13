import numpy as np
import pytest

from scipy.spatial.transform import Rotation

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system._center_of_mass import center_of_mass_groups
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.tools import traj_to_com_traj
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory

#: Volume per water molecule in Angstrom^3 (about 1 g/cm^3).
VOLUME_PER_MOLECULE = 30.0

WATER = np.array([
    [0.0, 0.0, 0.0],
    [0.9572, 0.0, 0.0],
    [-0.2399, 0.9266, 0.0],
])

#: (number of molecules, number of frames) of the water boxes.
SIZES = [(500, 1), (2000, 1), (8000, 1), (2000, 5)]
SIZE_IDS = [f"{n}mol-{f}frames" for n, f in SIZES]


def _water_topology(n_molecules):
    residue = Residue(
        name="WAT",
        elements=[Element("O"), Element("H"), Element("H")],
        residue_id=1,
        total_charge=0,
        atom_types=np.array([0, 1, 1]),
        partial_charges=np.zeros(3),
    )

    return Topology(
        atoms=[Atom(name) for name in "OHH" * n_molecules],
        residue_ids=np.ones(3 * n_molecules, dtype=int),
        reference_residues=[residue],
    )


def _water_positions(n_molecules, cell, rng):
    centers = rng.uniform(0.0, 1.0, (n_molecules, 1, 3)) * cell.box_lengths
    matrices = Rotation.random(n_molecules, random_state=rng).as_matrix()
    molecules = np.einsum("nij,kj->nki", matrices, WATER) + centers

    return cell.image(molecules.reshape(-1, 3))


def _water_box(n_molecules, n_frames, seed=42):
    """A seeded trajectory of random water molecules, built in memory."""

    rng = np.random.default_rng(seed)
    box = (n_molecules * VOLUME_PER_MOLECULE)**(1.0 / 3.0)
    cell = Cell(box, box, box)
    topology = _water_topology(n_molecules)

    return Trajectory([
        AtomicSystem(
            topology=topology,
            pos=_water_positions(n_molecules, cell, rng),
            cell=cell,
        ) for _ in range(n_frames)
    ])


@pytest.mark.benchmark(group="CenterOfMass")
class BenchmarkCenterOfMassGroups:

    @pytest.mark.parametrize("n_molecules", [500, 2000, 8000])
    def benchmark_groups(self, benchmark, n_molecules):
        frame = _water_box(n_molecules, 1)[0]
        masses = frame.atomic_masses
        group_ids = np.repeat(np.arange(n_molecules), 3)

        benchmark(center_of_mass_groups, frame.pos, masses, group_ids, frame.cell)

    def benchmark_groups_with_extent_check(self, benchmark):
        frame = _water_box(2000, 1)[0]
        masses = frame.atomic_masses
        group_ids = np.repeat(np.arange(2000), 3)

        benchmark(
            center_of_mass_groups,
            frame.pos,
            masses,
            group_ids,
            frame.cell,
            check_extent=True,
        )

    def benchmark_residues_of_an_atomic_system(self, benchmark):
        frame = _water_box(2000, 1)[0]

        benchmark(lambda: frame.center_of_mass_residues)


@pytest.mark.benchmark(group="TrajToComTraj")
class BenchmarkTrajToComTraj:

    @pytest.mark.parametrize("size", SIZES, ids=SIZE_IDS)
    def benchmark_group_of_three_atoms(self, benchmark, size):
        traj = _water_box(*size)

        benchmark(lambda: traj_to_com_traj(traj, group=3))

    @pytest.mark.parametrize("size", SIZES, ids=SIZE_IDS)
    def benchmark_by_residue(self, benchmark, size):
        traj = _water_box(*size)

        benchmark(lambda: traj_to_com_traj(traj, by_residue=True))

    def benchmark_by_residue_with_selection(self, benchmark):
        # the hydrogens of every molecule
        traj = _water_box(2000, 1)

        benchmark(
            lambda: traj_to_com_traj(
                traj, selection=Atom("H"), by_residue=True
            )
        )

    def benchmark_all_atoms(self, benchmark):
        traj = _water_box(2000, 5)

        benchmark(lambda: traj_to_com_traj(traj))
