import numpy as np
import pytest

from PQAnalysis.analysis.vacf import VACF, vacf_spectrum
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Element, Residue
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory, TrajectoryFormat


def _random_velocities(n_frames, n_atoms, seed=42):
    rng = np.random.default_rng(seed)

    return rng.uniform(-1.0, 1.0, (n_frames, n_atoms, 3))


def _make_trajectory(n_frames=400, n_atoms=50, seed=42):
    atoms = [Atom("H") for _ in range(n_atoms)]

    return Trajectory([
        AtomicSystem(atoms=atoms, vel=vel)
        for vel in _random_velocities(n_frames, n_atoms, seed=seed)
    ])


def _write_velocity_file(path, n_frames, n_atoms, seed=42):
    velocities = _random_velocities(n_frames, n_atoms, seed=seed)

    lines = []
    for frame in velocities:
        lines.append(f"{n_atoms} 20.0 22.0 24.0\n\n")
        for i, (x, y, z) in enumerate(frame):
            name = "O" if i % 2 == 0 else "H"
            lines.append(f"{name} {x:.7f} {y:.7f} {z:.7f}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)


def _write_charge_file(path, n_frames, n_atoms, seed=4711):
    rng = np.random.default_rng(seed)
    charges = 0.5 + 0.1 * rng.standard_normal((n_frames, n_atoms))

    lines = []
    for frame in charges:
        lines.append(f"{n_atoms}\n\n")
        for i, charge in enumerate(frame):
            name = "O" if i % 2 == 0 else "H"
            lines.append(f"{name} {charge:.7f}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)


@pytest.mark.benchmark(group="VACF")
@pytest.mark.usefixtures("kernel_impl")
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

    def benchmark_run_reader_fast_path(self, benchmark, tmp_path):
        # end-to-end fast path: raw-frame streaming from file plus
        # the Cython (or fallback) accumulation kernel
        filename = _write_velocity_file(tmp_path / "traj.vel", 2000, 100)

        benchmark(
            lambda: VACF(
                TrajectoryReader(filename),
                window_size=200,
                time_step=0.002,
                gap=10,
            ).run()
        )

    def benchmark_run_reader_fast_path_flux_static(self, benchmark, tmp_path):
        filename = _write_velocity_file(tmp_path / "traj.vel", 2000, 100)
        charges = np.tile([-0.8, 0.4], 50)

        benchmark(
            lambda: VACF(
                TrajectoryReader(filename),
                window_size=200,
                time_step=0.002,
                gap=10,
                charges=charges,
            ).run()
        )

    def benchmark_run_reader_fast_path_flux_charge_traj(
        self,
        benchmark,
        tmp_path,
    ):
        # charge-flux mode with the lockstep raw charge stream
        filename = _write_velocity_file(tmp_path / "traj.vel", 2000, 100)
        charge_filename = _write_charge_file(
            tmp_path / "traj.chrg", 2000, 100
        )

        benchmark(
            lambda: VACF(
                TrajectoryReader(filename),
                window_size=200,
                time_step=0.002,
                gap=10,
                charge_traj=TrajectoryReader(
                    charge_filename,
                    traj_format=TrajectoryFormat.CHARGE,
                ),
            ).run()
        )


@pytest.mark.benchmark(group="VACF")
class BenchmarkVACFSpectrum:

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


@pytest.mark.benchmark(group="VACF")
@pytest.mark.usefixtures("kernel_impl")
class BenchmarkVACFCenters:

    def benchmark_run_centers_reader_fast_path(self, benchmark, tmp_path):
        n_molecules = 200
        topology = _water_topology(n_molecules)
        velocities = _random_velocities(1000, 3 * n_molecules)
        lines = []
        for frame in velocities:
            lines.append(f"{3 * n_molecules}\n\n")
            for name, (x, y, z) in zip("OHH" * n_molecules, frame):
                lines.append(f"{name} {x:.7f} {y:.7f} {z:.7f}\n")
        filename = tmp_path / "water.vel"
        filename.write_text("".join(lines), encoding="utf-8")

        benchmark(
            lambda: VACF(
                TrajectoryReader(str(filename), topology=topology),
                window_size=100,
                time_step=0.002,
                gap=10,
                target_com=True,
            ).run()
        )

    def benchmark_run_centers_in_memory(self, benchmark):
        n_molecules = 100
        topology = _water_topology(n_molecules)
        traj = Trajectory([
            AtomicSystem(topology=topology, vel=vel)
            for vel in _random_velocities(300, 3 * n_molecules)
        ])

        benchmark(
            lambda: VACF(
                traj, window_size=50, time_step=0.002, gap=5,
                target_com=True,
            ).run()
        )
