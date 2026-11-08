import numpy as np
import pytest

from scipy.spatial.transform import Rotation

from PQAnalysis.analysis import RDF
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory


def _random_positions(n_frames, n_atoms, box=20.0, seed=42):
    rng = np.random.default_rng(seed)

    return rng.uniform(0.0, box, (n_frames, n_atoms, 3))


def _random_trajectory(n_frames, n_atoms, box=20.0, seed=42):
    cell = Cell(box, box, box)
    atoms = [Atom("O") if i % 2 == 0 else Atom("H") for i in range(n_atoms)]

    positions = _random_positions(n_frames, n_atoms, box=box, seed=seed)

    return Trajectory([
        AtomicSystem(atoms=atoms, pos=pos, cell=cell) for pos in positions
    ])


def _write_random_file(path, n_frames, n_atoms, box=20.0, seed=42):
    positions = _random_positions(n_frames, n_atoms, box=box, seed=seed)

    lines = []
    for frame_positions in positions:
        lines.append(f"{n_atoms} {box} {box} {box}\n\n")
        for i, (x, y, z) in enumerate(frame_positions):
            name = "O" if i % 2 == 0 else "H"
            lines.append(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)


@pytest.mark.benchmark(group="RDF")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkRDF:

    def benchmark_run(self, benchmark):
        traj = _random_trajectory(200, 100)

        benchmark(
            lambda: RDF(traj, "O", "H", delta_r=0.05, r_max=8.0).run()
        )

    def benchmark_run_no_intra_molecular(self, benchmark):
        traj = _random_trajectory(200, 100)

        benchmark(
            lambda: RDF(
                traj,
                "O",
                "H",
                delta_r=0.05,
                r_max=8.0,
                no_intra_molecular=True,
            ).run()
        )

    def benchmark_run_reader_fast_path(self, benchmark, tmp_path):
        # end-to-end fast path: header-only cell scan, raw-frame
        # streaming from file and the Cython (or fallback)
        # distance-histogram kernel
        filename = _write_random_file(tmp_path / "traj.xyz", 200, 100)

        benchmark(
            lambda: RDF(
                TrajectoryReader(filename),
                "O",
                "H",
                delta_r=0.05,
                r_max=8.0,
            ).run()
        )

    def benchmark_run_reader_fast_path_large(self, benchmark, tmp_path):
        filename = _write_random_file(tmp_path / "traj.xyz", 1000, 100)

        benchmark(
            lambda: RDF(
                TrajectoryReader(filename),
                "O",
                "H",
                delta_r=0.05,
                r_max=8.0,
            ).run()
        )



WATER = np.array([[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0], [-0.2399, 0.9266, 0.0]])


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


def _water_frames(n_molecules, n_frames, seed=42):
    """Seeded random waters in a box of liquid density, atoms wrapped."""

    rng = np.random.default_rng(seed)
    box = (n_molecules * 30.0)**(1.0 / 3.0)
    frames = []

    for _ in range(n_frames):
        centers = rng.uniform(0.0, box, (n_molecules, 1, 3))
        matrices = Rotation.from_quat(rng.normal(size=(n_molecules, 4))).as_matrix()
        molecules = np.einsum("nij,kj->nki", matrices, WATER) + centers
        frames.append(molecules.reshape(-1, 3) % box)

    return frames, box


@pytest.mark.benchmark(group="RDF")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkRDFCenters:

    def benchmark_run_centers_reader_fast_path(self, benchmark, tmp_path):
        # centers of mass of the molecules, streamed from file
        frames, box = _water_frames(500, 10)
        topology = _water_topology(500)
        lines = []
        for positions in frames:
            lines.append(f"{len(positions)} {box} {box} {box}\n\n")
            for name, (x, y, z) in zip("OHH" * 500, positions):
                lines.append(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")
        filename = tmp_path / "water.xyz"
        filename.write_text("".join(lines), encoding="utf-8")

        benchmark(
            lambda: RDF(
                TrajectoryReader(str(filename), topology=topology),
                "*",
                "*",
                delta_r=0.05,
                r_max=6.0,
                reference_com=True,
                target_com=True,
            ).run()
        )

    def benchmark_run_centers_in_memory(self, benchmark):
        frames, box = _water_frames(100, 5)
        topology = _water_topology(100)
        cell = Cell(box, box, box)
        traj = Trajectory([
            AtomicSystem(topology=topology, pos=positions, cell=cell)
            for positions in frames
        ])

        benchmark(
            lambda: RDF(
                traj,
                "*",
                "*",
                delta_r=0.05,
                r_max=3.0,
                reference_com=True,
                target_com=True,
            ).run()
        )
