import numpy as np
import pytest

from PQAnalysis.analysis.green_kubo import GreenKubo
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory


def _random_velocities(n_frames, n_atoms, seed=42):
    rng = np.random.default_rng(seed)

    return rng.uniform(-1.0e12, 1.0e12, (n_frames, n_atoms, 3))


def _make_trajectory(n_frames=400, n_atoms=50, seed=42):
    atoms = [Atom("Ar") for _ in range(n_atoms)]

    return Trajectory([
        AtomicSystem(atoms=atoms, vel=vel)
        for vel in _random_velocities(n_frames, n_atoms, seed=seed)
    ])


def _write_velocity_file(path, n_frames, n_atoms, seed=42):
    velocities = _random_velocities(n_frames, n_atoms, seed=seed)

    lines = []
    for frame in velocities:
        lines.append(f"{n_atoms} 100.0 100.0 100.0\n\n")
        for x, y, z in frame:
            lines.append(f"Ar {x:.7e} {y:.7e} {z:.7e}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)


@pytest.mark.benchmark(group="GreenKubo")
class BenchmarkGreenKubo:

    def benchmark_run_fft(self, benchmark):
        traj = _make_trajectory()

        benchmark(
            lambda: GreenKubo(
                traj,
                time_step=0.002,
                window_size=100,
                method="fft",
            ).run()
        )

    def benchmark_run_direct(self, benchmark):
        traj = _make_trajectory()

        benchmark(
            lambda: GreenKubo(
                traj,
                time_step=0.002,
                window_size=100,
                method="direct",
            ).run()
        )

    def benchmark_run_reader_fast_path(self, benchmark, tmp_path):
        # end-to-end fast path: raw-frame streaming from a velocity file
        filename = _write_velocity_file(tmp_path / "traj.vel", 2000, 100)

        benchmark(
            lambda: GreenKubo(
                TrajectoryReader(filename),
                time_step=0.002,
                window_size=200,
            ).run()
        )
