import numpy as np
import pytest

from PQAnalysis.analysis import RDF
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import TrajectoryReader
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
