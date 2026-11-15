"""
End-to-end tests of the analysis commands on the tutorial fixture
examples/water (one isolated H2O in a 10 A box, 25 frames), run through the
real ``pqanalysis`` command line.

The expected numbers are computed in the tests from the raw trajectory files
with plain NumPy, independent of PQAnalysis (the RDF histogram, the MSD and
VACF origins as documented, the momentum, and the vibrational frequencies from
the diagonalized mass weighted Hessian). Only the apodized VACF spectrum is a
committed reference, ref_data/spectrum.dat; regenerate it after a deliberate
change with

    cd examples/water && pqanalysis vacf vacf.in

and copy the written spectrum.dat to ref_data/. Every command runs with the
compiled and the fallback kernels and through the CSV and XVG exports.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.core import Atom

pytestmark = pytest.mark.integration

REFERENCE = Path(__file__).parent / "ref_data"
BOX = 10.0
N_FRAMES = 25
SPEED_OF_LIGHT = 2.99792458e10  # cm/s


def run_cli(monkeypatch, *arguments):
    monkeypatch.setattr(
        sys, "argv", ["pqanalysis", "--log-file", "off", *arguments]
    )
    pqanalysis_main()


def read_frames(filename):
    """Positions or velocities of all frames, shape (n_frames, n_atoms, 3)."""

    lines = Path(filename).read_text(encoding="utf-8").splitlines()
    n_atoms = int(lines[0].split()[0])
    size = n_atoms + 2

    return np.array([
        [[float(value) for value in line.split()[1:4]]
         for line in lines[i * size + 2:i * size + 2 + n_atoms]]
        for i in range(len(lines) // size)
    ])


def read_table(filename):
    return np.loadtxt(filename, comments="#")


def read_csv(filename):
    return np.loadtxt(filename, delimiter=",", skiprows=1)


def read_xvg(filename):
    """
    The table of an XVG file: every column is one data set (x against that
    column), the sets are separated by '&'.
    """

    columns, current = [], []
    for line in Path(filename).read_text(encoding="utf-8").splitlines():
        if line.startswith("&"):
            columns.append(current)
            current = []
        elif line.strip() and line[0] not in "#@":
            current.append(float(line.split()[1]))
    if current:
        columns.append(current)

    return np.array(columns).T


def sliding_origins(n_frames, window, gap):
    """The time origins: every gap-th frame (1-based) that has a full window."""

    return list(range(gap - 1, n_frames - window, gap))


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestRDF:

    def expected(self):
        positions = read_frames("trajectory.xyz")
        # the oxygen (atom 0) and its two hydrogens, minimum image
        difference = positions[:, [0]] - positions[:, 1:]
        difference -= BOX * np.round(difference / BOX)
        distance = np.linalg.norm(difference, axis=-1).ravel()
        edges = np.arange(0.0, 4.0 + 1e-9, 0.5)
        counts = np.histogram(distance, edges)[0]

        shell = 4.0 / 3.0 * np.pi * (edges[1:]**3 - edges[:-1]**3)
        target_density = 2.0 / BOX**3
        g = counts / (shell * target_density * 1 * N_FRAMES)

        return (edges[:-1] + edges[1:]) / 2, g, np.cumsum(counts) / N_FRAMES

    def test_matches_an_independent_calculation(
        self, test_integration_folder, kernel_impl, monkeypatch
    ):
        run_cli(monkeypatch, "rdf", "rdf.in")

        table = read_table("rdf.dat")
        r, g, coordination = self.expected()

        assert np.allclose(table[:, 0], r)
        assert np.allclose(table[:, 1], g, rtol=1e-10, atol=1e-12)
        assert np.allclose(table[:, 2], coordination)
        # each oxygen has two hydrogens
        assert table[-1, 2] == pytest.approx(2.0)
        assert g.max() > 100.0

    def test_exports(self, test_integration_folder, monkeypatch):
        run_cli(
            monkeypatch, "rdf", "rdf.in", "--export", "rdf.csv", "--export",
            "rdf.xvg",
        )

        native = read_table("rdf.dat")
        assert np.allclose(read_csv("rdf.csv"), native, rtol=1e-12)
        assert np.allclose(read_xvg("rdf.xvg"), native, rtol=1e-12)


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestMSD:

    def expected(self):
        positions = read_frames("trajectory.xyz")[:, 0]
        step = positions[1:] - positions[:-1]
        step -= BOX * np.round(step / BOX)
        unwrapped = np.vstack([positions[:1], positions[0] + np.cumsum(step, 0)])
        origins = sliding_origins(N_FRAMES, 8, 2)

        return np.array([
            np.mean([(unwrapped[t + lag] - unwrapped[t])**2 for t in origins], 0)
            for lag in range(9)
        ])

    def test_matches_an_independent_calculation(
        self, test_integration_folder, kernel_impl, monkeypatch
    ):
        run_cli(monkeypatch, "msd", "msd.in")

        table = read_table("msd.dat")

        assert np.array_equal(table[:, 0], np.arange(9))
        # the file has 8 decimals
        assert np.allclose(table[:, 1:], self.expected(), atol=1e-7)
        assert np.all(table[0, 1:] == 0.0)
        assert table[-1, 1:].sum() > 0.01

    def test_exports(self, test_integration_folder, monkeypatch):
        run_cli(
            monkeypatch, "msd", "msd.in", "--export", "msd.csv", "--export",
            "msd.xvg",
        )

        native = read_table("msd.dat")
        assert np.allclose(read_csv("msd.csv"), native, atol=1e-7)
        assert np.allclose(read_xvg("msd.xvg"), native, atol=1e-7)

    def test_the_diffusion_coefficients_are_those_of_the_fit(
        self, test_integration_folder, kernel_impl, monkeypatch, capsys
    ):
        run_cli(monkeypatch, "msd", "msd.in")

        table = read_table("msd.dat")
        time = table[-4:, 0] * 0.001  # ps, the last 4 points of the window
        printed = {}
        for line in capsys.readouterr().out.splitlines():
            parts = line.split()
            if parts and parts[0] in ("D_x", "D_y", "D_z", "D_total"):
                printed[parts[0]] = float(parts[2])

        # Einstein relation: D = slope / (2 d), 1 A^2/ps = 1e-8 m^2/s
        for axis, name in enumerate(("D_x", "D_y", "D_z")):
            slope = np.polyfit(time, table[-4:, 1 + axis], 1)[0]
            assert printed[name] == pytest.approx(slope / 2 * 1e-8, rel=1e-4)
        slope_total = np.polyfit(time, table[-4:, 1:].sum(axis=1), 1)[0]
        assert printed["D_total"] == pytest.approx(slope_total / 6 * 1e-8, rel=1e-4)


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestVACF:

    def expected(self):
        velocities = read_frames("trajectory.vel")
        origins = sliding_origins(N_FRAMES, 8, 2)

        return np.array([
            np.mean([
                (velocities[t + lag] * velocities[t]).sum()
                / (velocities[t]**2).sum() for t in origins
            ])
            for lag in range(9)
        ])

    def test_matches_an_independent_calculation(
        self, test_integration_folder, kernel_impl, monkeypatch
    ):
        run_cli(monkeypatch, "vacf", "vacf.in")

        table = read_table("vacf.dat")

        assert np.allclose(table[:, 0], np.arange(9) * 0.001)
        assert np.allclose(table[:, 1], self.expected(), atol=1e-7)
        assert table[0, 1] == 1.0

    def test_the_spectrum_equals_the_committed_reference(
        self, test_integration_folder, kernel_impl, monkeypatch
    ):
        run_cli(monkeypatch, "vacf", "vacf.in")

        spectrum = read_table("spectrum.dat")
        reference = read_table(REFERENCE / "spectrum.dat")

        assert spectrum.shape == reference.shape
        # the wavenumber axis follows from ftsize and the time step:
        # 2 * (ftsize - 1) samples of 0.001 ps
        assert np.allclose(
            spectrum[:, 0],
            np.arange(1, 17) / (2 * 15 * 0.001e-12 * SPEED_OF_LIGHT),
            rtol=1e-6,
        )
        assert np.allclose(spectrum, reference, rtol=1e-6, atol=1e-9)

    def test_exports(self, test_integration_folder, monkeypatch):
        run_cli(
            monkeypatch, "vacf", "vacf.in", "--export", "vacf.csv",
            "--export", "vacf.xvg",
        )

        native = read_table("vacf.dat")
        assert np.allclose(read_csv("vacf.csv"), native, atol=1e-7)
        assert np.allclose(read_xvg("vacf.xvg"), native, atol=1e-7)


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestCheckMomentum:

    ARGUMENTS = [
        "check_momentum", "trajectory.vel", "--selection", "all",
        "--output", "momentum.dat",
    ]

    def test_matches_an_independent_calculation(
        self, test_integration_folder, kernel_impl, monkeypatch
    ):
        run_cli(monkeypatch, *self.ARGUMENTS)

        velocities = read_frames("trajectory.vel")
        masses = np.array([Atom(name).mass for name in "OHH"])
        momentum = (velocities * masses[None, :, None]).sum(axis=1)
        expected = np.linalg.norm(momentum, axis=1) * 1e-15
        table = read_table("momentum.dat")

        assert np.array_equal(table[:, 0], np.arange(1, N_FRAMES + 1))
        # the total momentum of the fixture is zero up to rounding, so
        # compare against the size of the single atom terms
        scale = np.abs(velocities * masses[None, :, None]).max() * 1e-15
        assert np.allclose(table[:, 1], expected, rtol=1e-6, atol=1e-6 * scale)

    def test_the_scale_is_applied(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, *self.ARGUMENTS)
        run_cli(
            monkeypatch, "check_momentum", "trajectory.vel", "--selection",
            "all", "--scale", "1e-14", "--output", "scaled.dat",
        )

        assert np.allclose(
            read_table("scaled.dat")[:, 1], 10.0 * read_table("momentum.dat")[:, 1],
            rtol=1e-6, atol=1e-30,
        )

    def test_exports(self, test_integration_folder, monkeypatch):
        run_cli(
            monkeypatch, *self.ARGUMENTS, "--export", "momentum.csv",
            "--export", "momentum.xvg",
        )

        native = read_table("momentum.dat")
        assert np.allclose(read_csv("momentum.csv"), native, rtol=1e-10)
        assert np.allclose(read_xvg("momentum.xvg"), native, rtol=1e-10)


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
class TestVibrations:

    def expected(self):
        # the Hessian file holds force constants with the opposite sign
        # (hessian_sign = auto), in kcal mol^-1 A^-2
        hessian = -np.loadtxt("hessian.dat")
        masses = np.repeat([Atom(name).mass for name in "OHH"], 3)
        mass_weighted = (hessian + hessian.T) / 2 / np.sqrt(
            np.outer(masses, masses)
        )
        eigenvalues = np.linalg.eigvalsh(mass_weighted)
        # kcal -> J (4184), amu -> kg (1e-3 per mol), A^2 -> m^2 (1e-20)
        omega = np.sqrt(np.abs(eigenvalues) * 4184.0 / (1e-3 * 1e-20))

        return np.sort(np.sign(eigenvalues) * omega / (2 * np.pi * SPEED_OF_LIGHT))

    def test_matches_an_independent_calculation(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(monkeypatch, "vibrations", "vibrations.in")

        wavenumbers = read_table("wavenumbers.dat")[:, 0]
        expected = self.expected()

        assert len(wavenumbers) == 9
        # the three vibrations of water
        assert np.allclose(wavenumbers[-3:], expected[-3:], rtol=1e-4)
        assert 1400.0 < wavenumbers[-3] < 1700.0
        assert 3500.0 < wavenumbers[-2] < 3900.0
        assert 3500.0 < wavenumbers[-1] < 3900.0
        # the six translations and rotations are near zero
        assert np.abs(wavenumbers[:6]).max() < 100.0
        assert np.allclose(wavenumbers[:6], expected[:6], atol=1.0)

    def test_the_columns_are_complete(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "vibrations", "vibrations.in")

        table = read_table("wavenumbers.dat")

        # wavenumber, IR intensity, force constant, reduced mass
        assert table.shape == (9, 4)
        assert np.all(table[-3:, 2] > 0.0)  # force constants of the vibrations
        assert np.all(table[-3:, 3] > 0.0)  # reduced masses

    def test_exports(self, test_integration_folder, monkeypatch):
        run_cli(
            monkeypatch, "vibrations", "vibrations.in", "--export",
            "wavenumbers.csv", "--export", "wavenumbers.xvg",
        )

        native = read_table("wavenumbers.dat")
        assert np.allclose(read_csv("wavenumbers.csv"), native, rtol=1e-8)
        assert np.allclose(read_xvg("wavenumbers.xvg"), native, rtol=1e-8)
