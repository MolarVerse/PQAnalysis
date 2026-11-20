"""Scaling series: one analysis at growing numbers of atoms and frames.

A complexity regression (for example a quadratic loop over frames) shows up
as a changed growth between the points of a series, not only as a constant
factor of a single size.
"""

import numpy as np
import pytest

from PQAnalysis.analysis import MSD, RDF, VACF
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.traj import Trajectory

#: Number of atoms per 20 Angstrom box, so the density is the same at every size.
ATOMS_PER_BOX = 100

N_ATOMS = [100, 200, 400]
N_FRAMES = [500, 1000, 2000]


def _box(n_atoms):
    return 20.0 * (n_atoms / ATOMS_PER_BOX)**(1.0 / 3.0)


def _atoms(n_atoms):
    return [Atom("O") if i % 2 == 0 else Atom("H") for i in range(n_atoms)]


def _random_trajectory(n_frames, n_atoms, seed=42):
    rng = np.random.default_rng(seed)
    box = _box(n_atoms)
    cell = Cell(box, box, box)
    atoms = _atoms(n_atoms)

    return Trajectory([
        AtomicSystem(
            atoms=atoms, pos=rng.uniform(0.0, box, (n_atoms, 3)), cell=cell
        ) for _ in range(n_frames)
    ])


def _walk_trajectory(n_frames, n_atoms, seed=42):
    rng = np.random.default_rng(seed)
    box = _box(n_atoms)
    cell = Cell(box, box, box)
    atoms = _atoms(n_atoms)
    positions = np.cumsum(
        rng.normal(0.0, 0.1, (n_frames, n_atoms, 3)), axis=0
    ) % box

    return Trajectory([
        AtomicSystem(atoms=atoms, pos=pos, cell=cell) for pos in positions
    ])


def _velocity_trajectory(n_frames, n_atoms, seed=42):
    rng = np.random.default_rng(seed)
    atoms = [Atom("H") for _ in range(n_atoms)]

    return Trajectory([
        AtomicSystem(atoms=atoms, vel=rng.uniform(-1.0, 1.0, (n_atoms, 3)))
        for _ in range(n_frames)
    ])


@pytest.mark.benchmark(group="ScalingRDF")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkScalingRDF:

    @pytest.mark.parametrize("n_atoms", N_ATOMS)
    def benchmark_atoms(self, benchmark, n_atoms):
        traj = _random_trajectory(50, n_atoms)

        benchmark(lambda: RDF(traj, "O", "H", delta_r=0.05, r_max=8.0).run())

    @pytest.mark.parametrize("n_frames", [50, 100, 200])
    def benchmark_frames(self, benchmark, n_frames):
        traj = _random_trajectory(n_frames, 100)

        benchmark(lambda: RDF(traj, "O", "H", delta_r=0.05, r_max=8.0).run())


@pytest.mark.benchmark(group="ScalingMSD")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkScalingMSD:

    @pytest.mark.parametrize("n_atoms", N_ATOMS)
    def benchmark_atoms(self, benchmark, n_atoms):
        traj = _walk_trajectory(1000, n_atoms)

        benchmark(lambda: MSD(traj, "O", window=200, gap=10).run())

    @pytest.mark.parametrize("n_frames", N_FRAMES)
    def benchmark_frames(self, benchmark, n_frames):
        traj = _walk_trajectory(n_frames, 100)

        benchmark(lambda: MSD(traj, "O", window=200, gap=10).run())


@pytest.mark.benchmark(group="ScalingVACF")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkScalingVACF:

    @pytest.mark.parametrize("n_atoms", N_ATOMS)
    def benchmark_direct_atoms(self, benchmark, n_atoms):
        traj = _velocity_trajectory(400, n_atoms)

        benchmark(
            lambda: VACF(traj, window_size=100, time_step=0.002, gap=5).run()
        )

    @pytest.mark.parametrize("n_frames", [200, 400, 800])
    def benchmark_direct_frames(self, benchmark, n_frames):
        traj = _velocity_trajectory(n_frames, 50)

        benchmark(
            lambda: VACF(traj, window_size=100, time_step=0.002, gap=5).run()
        )

    @pytest.mark.parametrize("n_frames", [200, 400, 800])
    def benchmark_fft_frames(self, benchmark, n_frames):
        traj = _velocity_trajectory(n_frames, 50)

        benchmark(
            lambda: VACF(
                traj, window_size=100, time_step=0.002, method="fft"
            ).run()
        )
