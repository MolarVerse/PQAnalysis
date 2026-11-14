import numpy as np
import pytest

from PQAnalysis.analysis.msd import MSD
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory


def _random_walk_positions(n_frames, n_atoms, box=20.0, seed=42):
    rng = np.random.default_rng(seed)

    return np.cumsum(
        rng.normal(0.0, 0.1, (n_frames, n_atoms, 3)),
        axis=0,
    ) % box


def _random_walk_trajectory(n_frames, n_atoms, box=20.0, seed=42):
    cell = Cell(box, box, box)
    atoms = [Atom("O") if i % 2 == 0 else Atom("H") for i in range(n_atoms)]

    positions = _random_walk_positions(n_frames, n_atoms, box=box, seed=seed)

    return Trajectory([
        AtomicSystem(atoms=atoms, pos=pos, cell=cell) for pos in positions
    ])


def _write_random_walk_file(path, n_frames, n_atoms, box=20.0, seed=42):
    positions = _random_walk_positions(n_frames, n_atoms, box=box, seed=seed)

    lines = []
    for frame_positions in positions:
        lines.append(f"{n_atoms} {box} {box} {box}\n\n")
        for i, (x, y, z) in enumerate(frame_positions):
            name = "O" if i % 2 == 0 else "H"
            lines.append(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)


@pytest.mark.benchmark(group="MSD")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkMSD:

    def benchmark_run(self, benchmark):
        traj = _random_walk_trajectory(2000, 100)

        benchmark(lambda: MSD(traj, "O", window=200, gap=10).run())

    def benchmark_run_many_origins(self, benchmark):
        traj = _random_walk_trajectory(2000, 100)

        benchmark(lambda: MSD(traj, "O", window=1000, gap=10).run())

    def benchmark_run_reader_fast_path(self, benchmark, tmp_path):
        # end-to-end fast path: raw-frame streaming from file plus
        # the Cython (or fallback) accumulation kernel
        filename = _write_random_walk_file(
            tmp_path / "traj.xyz", 2000, 100
        )

        benchmark(
            lambda:
            MSD(TrajectoryReader(filename), "O", window=200, gap=10).run()
        )

    def benchmark_run_reader_fast_path_many_origins(self, benchmark, tmp_path):
        filename = _write_random_walk_file(
            tmp_path / "traj.xyz", 2000, 100
        )

        benchmark(
            lambda:
            MSD(TrajectoryReader(filename), "O", window=1000, gap=10).run()
        )



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


def _water_walk(n_molecules, n_frames, box=20.0, seed=42):
    """Waters whose centers do a seeded random walk, atoms wrapped."""

    rng = np.random.default_rng(seed)
    water = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])
    centers = np.cumsum(
        rng.normal(0.0, 0.1, (n_frames, n_molecules, 1, 3)), axis=0
    ) + rng.uniform(0.0, box, (1, n_molecules, 1, 3))

    return (centers + water).reshape(n_frames, -1, 3) % box, box


@pytest.mark.benchmark(group="MSD")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkMSDCenters:

    def benchmark_run_centers_reader_fast_path(self, benchmark, tmp_path):
        positions, box = _water_walk(100, 1000)
        topology = _water_topology(100)
        lines = []
        for frame in positions:
            lines.append(f"{len(frame)} {box} {box} {box}\n\n")
            for name, (x, y, z) in zip("OHH" * 100, frame):
                lines.append(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")
        filename = tmp_path / "water.xyz"
        filename.write_text("".join(lines), encoding="utf-8")

        benchmark(
            lambda: MSD(
                TrajectoryReader(str(filename), topology=topology),
                "*",
                window=100,
                gap=10,
                target_com=True,
            ).run()
        )

    def benchmark_run_centers_in_memory(self, benchmark):
        positions, box = _water_walk(50, 300)
        topology = _water_topology(50)
        cell = Cell(box, box, box)
        traj = Trajectory([
            AtomicSystem(topology=topology, pos=frame, cell=cell)
            for frame in positions
        ])

        benchmark(
            lambda: MSD(traj, "*", window=50, gap=10, target_com=True).run()
        )
