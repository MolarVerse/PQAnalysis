import numpy as np
import pytest

from PQAnalysis.analysis.momentum import Momentum
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory


def _random_velocities(n_frames, n_atoms, seed=42):
    rng = np.random.default_rng(seed)

    return rng.uniform(-1.0, 1.0, (n_frames, n_atoms, 3))


def _write_velocity_file(path, n_frames, n_atoms, seed=42):
    velocities = _random_velocities(n_frames, n_atoms, seed=seed)

    lines = []
    for frame_velocities in velocities:
        lines.append(f"{n_atoms}\n\n")
        for x, y, z in frame_velocities:
            lines.append(f"H {x:.6f} {y:.6f} {z:.6f}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)


@pytest.mark.benchmark(group="Momentum")
class BenchmarkMomentum:

    def benchmark_run(self, benchmark):
        atoms = [Atom("H") for _ in range(100)]
        frames = [
            AtomicSystem(atoms=atoms, vel=vel)
            for vel in _random_velocities(50, 100)
        ]
        traj = Trajectory(frames)

        benchmark(lambda: Momentum(traj).run())

    def benchmark_run_reader_fast_path(self, benchmark, tmp_path):
        # end-to-end fast path: raw-frame streaming from a velocity
        # file without per-frame AtomicSystem construction
        filename = _write_velocity_file(tmp_path / "traj.vel", 500, 100)

        benchmark(lambda: Momentum(TrajectoryReader(filename)).run())
