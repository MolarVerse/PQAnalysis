import numpy as np
import pytest

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.topology import Topology
from PQAnalysis.topology.shake_topology import ShakeTopologyGenerator
from PQAnalysis.traj import Trajectory

#: Volume per water molecule in Angstrom^3 (about 1 g/cm^3).
VOLUME_PER_MOLECULE = 30.0
WATER = np.array([
    [0.0, 0.0, 0.0],
    [0.9572, 0.0, 0.0],
    [-0.2399, 0.9266, 0.0],
])


def _residue():
    return Residue(
        name="WAT",
        elements=[Element("O"), Element("H"), Element("H")],
        residue_id=1,
        total_charge=0,
        atom_types=np.array([0, 1, 1]),
        partial_charges=np.zeros(3),
    )


def _topology(n_molecules):
    return Topology(
        atoms=[Atom(name) for name in "OHH" * n_molecules],
        residue_ids=np.ones(3 * n_molecules, dtype=int),
        reference_residues=[_residue()],
    )


def _trajectory(n_molecules, n_frames, seed=42):
    rng = np.random.default_rng(seed)
    box = (n_molecules * VOLUME_PER_MOLECULE)**(1.0 / 3.0)
    cell = Cell(box, box, box)
    topology = _topology(n_molecules)
    centers = rng.uniform(0.0, box, (n_molecules, 1, 3))

    frames = []
    for _ in range(n_frames):
        jitter = rng.normal(0.0, 0.02, (n_molecules, 3, 3))
        positions = (centers + WATER[None] + jitter).reshape(-1, 3)
        frames.append(
            AtomicSystem(
                topology=topology, pos=cell.image(positions), cell=cell
            )
        )

    return Trajectory(frames)


@pytest.mark.benchmark(group="Topology")
class BenchmarkTopology:

    @pytest.mark.parametrize("n_molecules", [500, 2000, 8000])
    def benchmark_build(self, benchmark, n_molecules):
        atoms = [Atom(name) for name in "OHH" * n_molecules]
        residue_ids = np.ones(3 * n_molecules, dtype=int)
        reference_residues = [_residue()]

        benchmark(
            Topology,
            atoms=atoms,
            residue_ids=residue_ids,
            reference_residues=reference_residues,
        )

    def benchmark_residue_atom_indices(self, benchmark):
        topology = _topology(2000)

        benchmark(lambda: topology.residue_atom_indices)

    def benchmark_slice(self, benchmark):
        topology = _topology(2000)
        indices = np.arange(0, topology.n_atoms, 2)

        benchmark(topology.__getitem__, indices)


@pytest.mark.benchmark(group="ShakeTopology")
class BenchmarkShakeTopology:

    @pytest.mark.parametrize("n_molecules, n_frames", [(500, 20), (2000, 5)])
    def benchmark_generate(self, benchmark, n_molecules, n_frames):
        trajectory = _trajectory(n_molecules, n_frames)

        def generate():
            ShakeTopologyGenerator("H").generate_topology(trajectory)

        benchmark(generate)
