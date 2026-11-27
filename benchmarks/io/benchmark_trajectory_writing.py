import numpy as np
import pytest

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import TrajectoryWriter
from PQAnalysis.io.conversion_api import traj2box, traj2extxyz, traj2qmcfc
from PQAnalysis.traj import Trajectory

SIZES = [(100, 200), (1000, 50), (5000, 10)]
SIZE_IDS = [f"{n}atoms-{f}frames" for n, f in SIZES]


def _trajectory(n_atoms, n_frames, seed=42):
    rng = np.random.default_rng(seed)
    atoms = [Atom(name) for name in rng.choice(["O", "H", "C"], n_atoms)]
    cell = Cell(30.0, 30.0, 30.0)

    return Trajectory([
        AtomicSystem(
            atoms=atoms, pos=rng.uniform(0.0, 30.0, (n_atoms, 3)), cell=cell
        ) for _ in range(n_frames)
    ])


def _write_xyz(path, n_atoms, n_frames):
    TrajectoryWriter(str(path)).write(_trajectory(n_atoms, n_frames))

    return str(path)


@pytest.mark.benchmark(group="TrajectoryWriter")
class BenchmarkTrajectoryWriter:

    @pytest.mark.parametrize("size", SIZES, ids=SIZE_IDS)
    def benchmark_write_xyz(self, benchmark, tmp_path, size):
        trajectory = _trajectory(*size)
        path = str(tmp_path / "out.xyz")

        benchmark(lambda: TrajectoryWriter(path, mode="o").write(trajectory))

    @pytest.mark.parametrize("engine_format", ["PQ", "QMCFC"])
    def benchmark_write_engine_format(self, benchmark, tmp_path, engine_format):
        trajectory = _trajectory(1000, 50)
        path = str(tmp_path / "out.xyz")

        benchmark(
            lambda: TrajectoryWriter(
                path, engine_format=engine_format, mode="o"
            ).write(trajectory)
        )


@pytest.mark.benchmark(group="TrajectoryConversion")
class BenchmarkTrajectoryConversion:

    @pytest.fixture(name="source")
    def source_fixture(self, tmp_path):
        return _write_xyz(tmp_path / "in.xyz", 1000, 50)

    @pytest.mark.parametrize("profile", ["ase", "nep", "pq"])
    def benchmark_traj2extxyz(self, benchmark, tmp_path, source, profile):
        output = str(tmp_path / "out.extxyz")

        benchmark(
            traj2extxyz, [source], output, "o", extxyz_profile=profile
        )

    def benchmark_traj2qmcfc(self, benchmark, tmp_path, source):
        output = str(tmp_path / "out.xyz")

        benchmark(traj2qmcfc, [source], output, "o")

    def benchmark_traj2box(self, benchmark, tmp_path, source):
        output = str(tmp_path / "out.box")

        benchmark(traj2box, [source], False, output, "o")
