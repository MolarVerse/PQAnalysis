"""
Integration test of the VACF of the centers of mass of molecules through the
input file and the vacf API, on the 2000 molecule TIP3P water box of
examples/h2o_mm-mshake.

The velocity of every molecule rotates in the xy plane at a common frequency
(random amplitude and phase) and every atom adds an internal motion with zero
total momentum, so the normalized VACF of the centers is exactly
cos(omega * lag), while the VACF of the atoms is not.
"""

import shutil

import numpy as np
import pytest

from PQAnalysis.analysis.vacf import vacf
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Cell
from PQAnalysis.io import TrajectoryWriter, read_restart_file
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.tools import traj2com
from PQAnalysis.traj import Trajectory

pytestmark = pytest.mark.integration

N_MOLECULES = 2000
N_FRAMES = 60
WINDOW = 10
GAP = 2
OMEGA = 0.3


def write_input(name, **keys):
    lines = {
        "traj_files": "water.vel",
        "target_selection": "all",
        "window": WINDOW,
        "gap": GAP,
        "time_step": 0.5,
        "out_file": f"{name}.dat",
        "log_file": f"{name}.log",
    }
    lines.update(keys)

    with open(f"{name}.in", "w", encoding="utf-8") as file:
        for key, value in lines.items():
            if value is not None:
                print(f"{key} = {value}", file=file)

    return f"{name}.in"


def read_table(name):
    """The columns time and correlation of the output file."""

    return np.loadtxt(f"{name}.dat", comments=("#", "!"))


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestVACFCentersInputFile:

    @pytest.fixture(autouse=True)
    def rotating_water(self, test_integration_folder):
        system = read_restart_file("input.rst", "moldescriptor.dat")
        topology = system.topology
        masses = system.atomic_masses
        molecule = topology.residue_numbers
        rng = np.random.default_rng(6)
        amplitude = rng.uniform(0.5, 2.0, N_MOLECULES)
        phase = rng.uniform(0.0, 2.0 * np.pi, N_MOLECULES)

        frames = []
        for frame in range(N_FRAMES):
            angle = OMEGA * frame + phase
            center = np.stack(
                [amplitude * np.cos(angle), amplitude * np.sin(angle),
                 np.zeros(N_MOLECULES)], axis=1,
            )
            internal = rng.normal(size=(len(masses), 3))
            momentum = np.add.reduceat(
                internal * masses[:, None], np.arange(0, len(masses), 3)
            )
            internal -= np.repeat(
                momentum / np.add.reduceat(masses, np.arange(0, len(masses), 3))[:, None],
                3, axis=0,
            )
            frames.append(
                AtomicSystem(
                    topology=topology, vel=center[molecule] + internal,
                    cell=Cell(),
                )
            )
        TrajectoryWriter("water.vel").write(Trajectory(frames), traj_type="vel")

    def centers_keys(self, **extra):
        keys = {
            "restart_file": "input.rst",
            "moldescriptor_file": "moldescriptor.dat",
            "target_com": "True",
        }
        keys.update(extra)

        return keys

    def test_the_keys_give_the_analytic_vacf(self):
        vacf(write_input("com", **self.centers_keys()))

        table = read_table("com")

        assert np.allclose(
            table[:, 1], np.cos(OMEGA * np.arange(WINDOW + 1)), atol=1e-5
        )
        assert np.allclose(table[:, 0], np.arange(WINDOW + 1) * 0.5)

    def test_the_vacf_of_the_atoms_is_not_the_one_of_the_centers(self):
        vacf(write_input("atoms"))

        assert not np.allclose(
            read_table("atoms")[:, 1],
            np.cos(OMEGA * np.arange(WINDOW + 1)), atol=1e-2,
        )

    def test_equals_the_manual_center_of_mass_workflow(self):
        # the manual workflow the keys replace: a velocity file of the centers
        traj2com(
            ["water.vel"], output="com.vel", by_residue=True,
            restart_file="input.rst", moldescriptor_file="moldescriptor.dat",
        )
        vacf(write_input("manual", traj_files="com.vel"))
        vacf(write_input("keys", **self.centers_keys()))

        assert np.allclose(
            read_table("keys"), read_table("manual"), atol=1e-6
        )

    def test_the_log_counts_molecules(self):
        vacf(write_input("com", **self.centers_keys()))

        log = open("com.log", encoding="utf-8").read()

        assert (
            "total number of molecules (centers of mass) in target "
            "selection: 2000" in log
        )

    def test_the_fft_method(self):
        vacf(write_input("fft", method="fft", **self.centers_keys()))

        # the fft estimator averages over all origins: the same rotation
        assert np.allclose(
            read_table("fft")[:, 1], np.cos(OMEGA * np.arange(WINDOW + 1)),
            atol=1e-2,
        )

    def test_the_topology_files_are_inferred(self):
        shutil.copy("input.rst", "water.rst")

        vacf(write_input("inferred", target_com="True"))
        vacf(
            write_input(
                "explicit", **self.centers_keys(restart_file="water.rst")
            )
        )

        assert np.array_equal(read_table("inferred"), read_table("explicit"))

    def test_without_the_key_the_atom_vacf_is_unchanged(self):
        vacf(write_input("atoms"))
        vacf(write_input("false", target_com="False"))

        assert np.array_equal(read_table("atoms"), read_table("false"))

    def test_the_key_needs_the_topology_files(self):
        with pytest.raises(InputFileError, match="target_com key requires"):
            vacf(write_input("none", target_com="True"))

    def test_the_key_cannot_be_used_with_the_charge_flux(self):
        with pytest.raises(InputFileError, match="no charge flux"):
            vacf(
                write_input(
                    "flux", charge_file="charges.dat", **self.centers_keys()
                )
            )
