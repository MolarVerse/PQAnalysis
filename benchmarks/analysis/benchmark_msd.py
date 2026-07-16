import numpy as np
import pytest

from PQAnalysis.analysis.msd import MSD
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.traj import Trajectory


def _random_walk_trajectory(n_frames, n_atoms, box=20.0, seed=42):
    rng = np.random.default_rng(seed)
    cell = Cell(box, box, box)
    atoms = [Atom("O") if i % 2 == 0 else Atom("H") for i in range(n_atoms)]

    positions = np.cumsum(
        rng.normal(0.0, 0.1, (n_frames, n_atoms, 3)),
        axis=0,
    ) % box

    return Trajectory([
        AtomicSystem(atoms=atoms, pos=pos, cell=cell) for pos in positions
    ])


@pytest.mark.benchmark(group="MSD")
class BenchmarkMSD:

    def benchmark_run(self, benchmark):
        traj = _random_walk_trajectory(2000, 100)

        benchmark(lambda: MSD(traj, "O", window=200, gap=10).run())

    def benchmark_run_many_origins(self, benchmark):
        traj = _random_walk_trajectory(2000, 100)

        benchmark(lambda: MSD(traj, "O", window=1000, gap=10).run())
