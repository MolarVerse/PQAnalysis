"""
Integration test of the MSD of the centers of mass of molecules through the
input file and the msd API, on the 2000 molecule TIP3P water box of
examples/h2o_mm-mshake.

The molecules of the box move with a constant velocity of up to 0.3 box per
frame and cross the periodic boundary. In addition every oxygen drifts away
from its molecule with a small constant velocity, so atoms and centers move
differently. The velocity of a center is the mass weighted mean of the
velocities of its atoms, so the MSD of the centers is known exactly:
MSD(lag) = mean over the molecules of (v_center * lag)^2, and so is that of
the oxygens, which have the velocity of their molecule plus the drift.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.msd import msd
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.io import TrajectoryWriter, read_restart_file
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.tools import traj2com
from PQAnalysis.traj import Trajectory

pytestmark = pytest.mark.integration

N_FRAMES = 14
WINDOW = 10
GAP = 2


def write_input(name, **keys):
    lines = {
        "traj_files": "water.xyz",
        "target_selection": "all",
        "window": WINDOW,
        "gap": GAP,
        "time_step": 0.002,
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
    """The columns lag, x, y and z of the output file."""

    return np.loadtxt(f"{name}.dat", comments=("#", "!"))


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestMSDCentersInputFile:

    @pytest.fixture(autouse=True)
    def moving_water(self, test_integration_folder):
        system = read_restart_file("input.rst", "moldescriptor.dat")
        cell = system.cell
        rng = np.random.default_rng(8)
        fractional = rng.uniform(-0.3, 0.3, (2000, 3))
        molecule_velocities = fractional * cell.box_lengths
        drift = rng.uniform(-0.02, 0.02, (2000, 3)) * cell.box_lengths

        masses = system.atomic_masses.reshape(-1, 3)
        oxygen_weight = masses[:, 0] / masses.sum(axis=1)

        per_atom = np.repeat(molecule_velocities, 3, axis=0)
        per_atom[::3] += drift
        # the velocity of a center and of an oxygen
        self.velocities = molecule_velocities + oxygen_weight[:, None] * drift
        self.oxygen_velocities = molecule_velocities + drift

        frames = [
            AtomicSystem(
                topology=system.topology,
                pos=cell.image(system.pos + frame * per_atom),
                cell=cell,
            ) for frame in range(N_FRAMES)
        ]
        TrajectoryWriter("water.xyz").write(Trajectory(frames))

    def analytic(self, velocities=None):
        lags = np.arange(WINDOW + 1)
        velocities = self.velocities if velocities is None else velocities

        return (velocities**2).mean(axis=0)[None, :] * lags[:, None]**2

    def centers_keys(self, **extra):
        keys = {
            "restart_file": "input.rst",
            "moldescriptor_file": "moldescriptor.dat",
            "target_com": "True",
        }
        keys.update(extra)

        return keys

    def test_the_keys_give_the_analytic_msd(self):
        msd(write_input("com", **self.centers_keys()))

        table = read_table("com")

        assert np.allclose(table[:, 1:4], self.analytic(), rtol=1e-4, atol=1e-6)
        # many box lengths squared
        assert table[-1, 1:4].sum() > 100.0

    def test_equals_the_manual_center_of_mass_workflow(self):
        # the manual workflow the keys replace
        traj2com(
            ["water.xyz"], output="com.xyz", by_residue=True,
            restart_file="input.rst", moldescriptor_file="moldescriptor.dat",
        )
        msd(write_input("manual", traj_files="com.xyz"))
        msd(write_input("keys", **self.centers_keys()))

        manual, keys = read_table("manual"), read_table("keys")

        assert keys.shape == manual.shape
        assert np.allclose(keys, manual, rtol=1e-6, atol=1e-8)

    def test_the_log_counts_molecules(self):
        msd(write_input("com", **self.centers_keys()))

        log = open("com.log", encoding="utf-8").read()

        assert (
            "total number of molecules (centers of mass) in target "
            "selection: 2000" in log
        )

    def test_the_topology_files_are_inferred(self):
        import shutil
        shutil.copy("input.rst", "water.rst")

        msd(write_input("inferred", target_com="True"))
        msd(write_input("explicit", **self.centers_keys(restart_file="water.rst")))

        assert np.array_equal(read_table("inferred"), read_table("explicit"))

    def test_the_msd_of_the_oxygens_is_not_that_of_the_centers(self):
        # without target_com the same input gives the atom MSD
        msd(write_input("atoms", target_selection="O"))
        msd(write_input("com", **self.centers_keys()))

        atoms, centers = read_table("atoms"), read_table("com")
        oxygen = self.analytic(self.oxygen_velocities)

        assert np.allclose(atoms[:, 1:4], oxygen, rtol=1e-4, atol=1e-6)
        assert np.allclose(centers[:, 1:4], self.analytic(), rtol=1e-4, atol=1e-6)
        assert np.abs(atoms[-1, 1:4] / centers[-1, 1:4] - 1.0).max() > 1e-3

    def test_the_selection_defines_the_center(self):
        # the center of the oxygens of a molecule is the oxygen
        msd(write_input("atoms", target_selection="O"))
        msd(write_input("com", **self.centers_keys(target_selection="O")))

        assert np.allclose(
            read_table("atoms"), read_table("com"), rtol=1e-6, atol=1e-8
        )

    def test_without_the_key_the_atom_msd_is_unchanged(self):
        msd(write_input("atoms"))
        msd(write_input("false", target_com="False"))

        assert np.array_equal(read_table("atoms"), read_table("false"))

    def test_the_key_needs_the_topology_files(self):
        with pytest.raises(InputFileError, match="target_com key requires"):
            msd(write_input("none", target_com="True"))
