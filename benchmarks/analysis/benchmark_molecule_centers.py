import numpy as np
import pytest

from PQAnalysis.analysis._molecule_centers import MoleculeCenters
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.topology import Topology

#: Volume per water molecule in Angstrom^3 (about 1 g/cm^3).
VOLUME_PER_MOLECULE = 30.0
WATER = np.array([
    [0.0, 0.0, 0.0],
    [0.9572, 0.0, 0.0],
    [-0.2399, 0.9266, 0.0],
])


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


def _from_topology(n_molecules, topology):
    return MoleculeCenters.from_topology(
        topology, np.arange(3 * n_molecules), None, ValueError
    )


def _water_centers(n_molecules, seed=42):
    centers = _from_topology(n_molecules, _water_topology(n_molecules))

    rng = np.random.default_rng(seed)
    box = (n_molecules * VOLUME_PER_MOLECULE)**(1.0 / 3.0)
    cell = Cell(box, box, box)
    positions = (
        rng.uniform(0.0, box, (n_molecules, 1, 3)) + WATER[None]
    ).reshape(-1, 3)

    return centers, cell.image(positions), cell, rng


@pytest.mark.benchmark(group="MoleculeCenters")
class BenchmarkMoleculeCenters:

    @pytest.mark.parametrize("n_molecules", [500, 2000, 8000])
    def benchmark_from_topology(self, benchmark, n_molecules):
        topology = _water_topology(n_molecules)

        benchmark(_from_topology, n_molecules, topology)

    @pytest.mark.parametrize("n_molecules", [500, 2000, 8000])
    def benchmark_extend(self, benchmark, n_molecules):
        centers, positions, cell, _ = _water_centers(n_molecules)

        benchmark(centers.extend, positions, cell)

    @pytest.mark.parametrize("n_frames", [1, 50])
    def benchmark_mean_velocities(self, benchmark, n_frames):
        centers, _, _, rng = _water_centers(2000)
        velocities = rng.uniform(-1.0, 1.0, (n_frames, 6000, 3))

        benchmark(centers.mean_velocities, velocities)
