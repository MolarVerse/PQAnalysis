import numpy as np
import pytest

from PQAnalysis.analysis.vacf import VACF, vacf_spectrum
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.traj import Trajectory


def _make_trajectory(n_frames=400, n_atoms=50, seed=42):
    rng = np.random.default_rng(seed)
    atoms = [Atom("H") for _ in range(n_atoms)]
    frames = [
        AtomicSystem(atoms=atoms, vel=rng.uniform(-1.0, 1.0, (n_atoms, 3)))
        for _ in range(n_frames)
    ]
    return Trajectory(frames)


@pytest.mark.benchmark(group="VACF")
class BenchmarkVACF:

    def benchmark_run_direct(self, benchmark):
        traj = _make_trajectory()

        benchmark(
            lambda: VACF(
                traj,
                window_size=100,
                time_step=0.002,
                gap=5,
            ).run()
        )

    def benchmark_run_fft(self, benchmark):
        traj = _make_trajectory()

        benchmark(
            lambda: VACF(
                traj,
                window_size=100,
                time_step=0.002,
                method="fft",
            ).run()
        )

    def benchmark_spectrum(self, benchmark):
        time = np.arange(1001) * 0.002
        correlation = np.cos(2.0 * np.pi * 25.0 * time) * np.exp(-2.0 * time)

        benchmark(
            lambda: vacf_spectrum(
                time,
                correlation,
                ftsize=2000,
                window_function="blackman",
                window_start=0.5,
                window_stop=1.5,
            )
        )
