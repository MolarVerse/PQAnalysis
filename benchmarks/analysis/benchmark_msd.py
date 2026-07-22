import numpy as np
import pytest

from PQAnalysis.analysis.msd import MSD
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import TrajectoryReader
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
