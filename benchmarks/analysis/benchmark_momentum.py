import numpy as np
import pytest

from PQAnalysis.analysis.momentum import Momentum
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.traj import Trajectory



@pytest.mark.benchmark(group="Momentum")
class BenchmarkMomentum:

    def benchmark_run(self, benchmark):
        rng = np.random.default_rng(42)
        atoms = [Atom("H") for _ in range(100)]
        frames = [
            AtomicSystem(atoms=atoms, vel=rng.uniform(-1.0, 1.0, (100, 3)))
            for _ in range(50)
        ]
        traj = Trajectory(frames)

        benchmark(lambda: Momentum(traj).run())
