"""
End-to-end tests of build_spectrum, build_nep_traj and continue_input through
the real command line.

The expected numbers are computed in the tests with plain NumPy from the
documented formulas and from the input files:

* build_spectrum: Gaussian ``sum_k I_k exp(-alpha (g - nu_k)^2)`` and
  Lorentzian ``sum_k I_k gamma^2 / ((g - nu_k)^2 + gamma^2)`` with
  ``gamma^2 = ln(2) / alpha``, ``alpha = 4 ln(2) / fwhm^2`` for the full width
  at half maximum.
* build_nep_traj: energies and forces in eV and eV/A from kcal/mol, the
  lattice from the cell and the split of the frames into the files.
* continue_input: the next input file differs from the previous one only in
  the file names, whose counters are incremented.
"""

import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.io.exceptions import FileWritingModeError

pytestmark = pytest.mark.integration

#: kcal/mol in eV: 4184 J per kcal, divided by Avogadro and the elementary charge
EV_PER_KCAL_PER_MOL = 4184.0 / (6.02214076e23 * 1.602176634e-19)


def run_cli(monkeypatch, *arguments):
    monkeypatch.setattr(
        sys, "argv", ["pqanalysis", "--log-file", "off", *arguments]
    )
    pqanalysis_main()


def read_table(filename):
    return np.loadtxt(filename, comments="#")


def gaussian(grid, sticks, alpha):
    delta = grid[:, None] - sticks[None, :, 0]

    return (sticks[None, :, 1] * np.exp(-alpha * delta**2)).sum(axis=1)


def lorentzian(grid, sticks, alpha):
    gamma_squared = math.log(2.0) / alpha
    delta = grid[:, None] - sticks[None, :, 0]

    return (
        sticks[None, :, 1] * gamma_squared / (delta**2 + gamma_squared)
    ).sum(axis=1)


STICKS = np.array([
    [1000.0, 1.0],
    [1500.0, 3.0],
    [3000.0, 0.5],
])


@pytest.fixture(name="sticks")
def sticks_fixture(tmpdir):  # pylint: disable=unused-argument
    Path("sticks.dat").write_text(
        "\n".join(f"{nu} {intensity}" for nu, intensity in STICKS) + "\n",
        encoding="utf-8",
    )

    return STICKS


class TestBuildSpectrum:

    def test_the_default_gaussian(self, sticks, monkeypatch):
        run_cli(monkeypatch, "build_spectrum", "sticks.dat", "-o", "out.dat")

        table = read_table("out.dat")
        grid = np.arange(10.0, 4000.0, 0.25)

        assert np.allclose(table[:, 0], grid)
        assert np.allclose(table[:, 1], gaussian(grid, sticks, 0.0025), rtol=1e-9, atol=1e-12)
        # the maximum of a stick at 1500 cm^-1 with the intensity 3
        assert table[:, 1].max() == pytest.approx(3.0, rel=1e-3)
        assert grid[table[:, 1].argmax()] == pytest.approx(1500.0)

    @pytest.mark.parametrize(
        "kernel, function", [("--lorentzian", lorentzian), (None, gaussian)]
    )
    def test_the_alpha_and_the_grid(self, sticks, monkeypatch, kernel, function):
        arguments = [
            "build_spectrum", "sticks.dat", "--alpha", "0.01", "--min", "900",
            "--max", "1600", "--step", "5", "-o", "out.dat",
        ]
        if kernel:
            arguments.append(kernel)
        run_cli(monkeypatch, *arguments)

        table = read_table("out.dat")
        grid = np.arange(900.0, 1600.0, 5.0)

        assert np.allclose(table[:, 0], grid)
        assert np.allclose(table[:, 1], function(grid, sticks, 0.01), rtol=1e-9, atol=1e-12)

    def test_the_width_as_a_full_width_at_half_maximum(self, sticks, monkeypatch):
        run_cli(
            monkeypatch, "build_spectrum", "sticks.dat", "--fwhm", "50",
            "--min", "1000", "--max", "1100", "--step", "1", "-o", "out.dat",
        )

        table = read_table("out.dat")
        alpha = 4.0 * math.log(2.0) / 50.0**2

        assert np.allclose(table[:, 1], gaussian(table[:, 0], sticks, alpha), rtol=1e-9, atol=1e-12)
        # a Gaussian of the width 50 has half its maximum 25 cm^-1 away
        assert gaussian(np.array([1025.0]), STICKS[:1], alpha)[0] == pytest.approx(0.5)

    def test_the_stick_file_may_have_a_header_and_more_columns(
        self, tmpdir, monkeypatch
    ):
        # the output of the vibrations command: header lines and 4 columns
        Path("vibrations.dat").write_text(
            "# PQAnalysis: Vibrational analysis\n"
            "# FIELDS wavenumber ir_intensity force_constant reduced_mass\n"
            "1492.5\t150.5\t1.4\t1.08\n3669.0\t86.4\t8.3\t1.05\n",
            encoding="utf-8",
        )

        run_cli(
            monkeypatch, "build_spectrum", "vibrations.dat", "-o", "out.dat"
        )

        table = read_table("out.dat")
        expected = gaussian(
            table[:, 0], np.array([[1492.5, 150.5], [3669.0, 86.4]]), 0.0025
        )
        assert np.allclose(table[:, 1], expected, rtol=1e-9, atol=1e-12)

    def test_the_exports(self, sticks, monkeypatch):
        run_cli(
            monkeypatch, "build_spectrum", "sticks.dat", "--max", "200", "-o",
            "out.dat", "--export", "out.csv", "--export", "out.xvg",
        )

        native = read_table("out.dat")
        csv = np.loadtxt("out.csv", delimiter=",", skiprows=1)

        assert np.allclose(csv, native, rtol=1e-10)
        assert Path("out.xvg").is_file()

    def test_the_writing_modes(self, sticks, monkeypatch):
        run_cli(monkeypatch, "build_spectrum", "sticks.dat", "--max", "50", "-o", "out.dat")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "build_spectrum", "sticks.dat", "--max", "50", "-o", "out.dat")

        run_cli(
            monkeypatch, "build_spectrum", "sticks.dat", "--max", "50", "-o",
            "out.dat", "--mode", "o",
        )


@pytest.mark.parametrize("example_dir", ["water"], indirect=False)
def test_the_spectrum_of_the_vibrations_of_water(
    test_integration_folder, monkeypatch
):
    # the pipeline: vibrations of water, then a broadened spectrum of its sticks
    run_cli(monkeypatch, "vibrations", "vibrations.in")
    run_cli(
        monkeypatch, "build_spectrum", "wavenumbers.dat", "--fwhm", "30",
        "--min", "1000", "--max", "4000", "--step", "2", "-o", "ir.dat",
    )

    sticks = read_table("wavenumbers.dat")[:, :2]
    table = read_table("ir.dat")
    alpha = 4.0 * math.log(2.0) / 30.0**2

    assert np.allclose(table[:, 1], gaussian(table[:, 0], sticks, alpha), rtol=1e-8, atol=1e-10)
    # the bending band of water near 1500 cm^-1
    region = table[(table[:, 0] > 1400) & (table[:, 0] < 1600)]
    assert region[:, 0][region[:, 1].argmax()] == pytest.approx(1492.5, abs=2.0)


N_FRAMES = 6
TOLERANCE = {"rtol": 1e-5, "atol": 1e-8}


@pytest.fixture(name="nep_data")
def nep_data_fixture(tmpdir):  # pylint: disable=unused-argument
    """Six frames of C and H with energies and forces, as PQ writes them."""

    rng = np.random.default_rng(1)
    positions = rng.uniform(0.0, 5.0, (N_FRAMES, 2, 3))
    forces = rng.normal(size=(N_FRAMES, 2, 3))
    energies = -185900.0 - 10.0 * np.arange(1, N_FRAMES + 1)

    def frames(values):
        return "".join(
            "2 10.0 10.0 10.0 90.0 90.0 90.0\n#\n"
            + "".join(
                f"{name} {x:.6f} {y:.6f} {z:.6f}\n"
                for name, (x, y, z) in zip("CH", frame)
            ) for frame in values
        )

    Path("md.xyz").write_text(frames(positions), encoding="utf-8")
    Path("md.frc").write_text(frames(forces), encoding="utf-8")
    Path("md.info").write_text(
        Path(__file__).parents[2].joinpath(
            "data", "readEnergyFile", "md-01.info"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    Path("md.instant_en").write_text(
        "".join(
            f"{i} 298.0 23000.0 -185898.0 {energy:.1f} 0.0 5.6 0.0 3729.1 "
            "0.83 1.1e-16 1.0\n" for i, energy in enumerate(energies, 1)
        ),
        encoding="utf-8",
    )
    # the rounded values of the files
    return (
        np.round(positions, 6), np.round(forces, 6), np.round(energies, 1)
    )


def read_nep(filename):
    """The frames (energy, lattice, species, data columns) of a NEP file."""

    lines = Path(filename).read_text(encoding="utf-8").splitlines()
    frames, start = [], 0
    while start < len(lines):
        n_atoms = int(lines[start])
        header = lines[start + 1]
        energy = float(re.search(r"energy=(\S+)", header).group(1))
        lattice = np.array(
            re.search(r'lattice="([^"]*)"', header).group(1).split(), dtype=float
        ).reshape(3, 3)
        rows = [line.split() for line in lines[start + 2:start + 2 + n_atoms]]
        frames.append((
            energy,
            lattice,
            [row[0] for row in rows],
            np.array([[float(v) for v in row[1:]] for row in rows]),
            header,
        ))
        start += 2 + n_atoms

    return frames


class TestBuildNepTraj:

    def test_every_frame_with_its_energy(self, nep_data, monkeypatch):
        positions, _, energies = nep_data

        run_cli(monkeypatch, "build_nep_traj", "md", "-o", "nep.xyz")

        frames = read_nep("nep.xyz")

        assert len(frames) == N_FRAMES
        for frame, position, energy in zip(frames, positions, energies):
            assert frame[0] == pytest.approx(energy * EV_PER_KCAL_PER_MOL, rel=1e-5)
            assert np.allclose(frame[1], np.diag([10.0, 10.0, 10.0]), atol=1e-9)
            assert frame[2] == ["C", "H"]
            assert np.allclose(frame[3], position, atol=1e-7)
            assert "forces" not in frame[4]

    def test_the_forces(self, nep_data, monkeypatch):
        _, forces, _ = nep_data

        run_cli(monkeypatch, "build_nep_traj", "md", "--use-forces", "-o", "nep.xyz")

        frames = read_nep("nep.xyz")

        assert len(frames) == N_FRAMES
        for frame, force in zip(frames, forces):
            assert "forces:R:3" in frame[4]
            # columns: the positions, then the forces in eV/A
            assert np.allclose(
                frame[3][:, 3:6], force * EV_PER_KCAL_PER_MOL, **TOLERANCE
            )

    def test_the_test_ratio_splits_the_frames(self, nep_data, monkeypatch):
        _, _, energies = nep_data

        run_cli(
            monkeypatch, "build_nep_traj", "md", "--test-ratio", "0.5", "-o",
            "nep.xyz",
        )

        train = read_nep("nep.xyz_train")
        test = read_nep("nep.xyz_test")
        expected = sorted(energies * EV_PER_KCAL_PER_MOL)

        assert len(train) == len(test) == 3
        # every frame is in exactly one of the files
        assert sorted(f[0] for f in train + test) == pytest.approx(expected, rel=1e-5)

    def test_the_total_ratios(self, nep_data, monkeypatch):
        _, _, energies = nep_data

        run_cli(
            monkeypatch, "build_nep_traj", "md", "--total-ratios", "4:1:1",
            "-o", "nep.xyz",
        )

        parts = {
            name: read_nep(f"nep.xyz_{name}")
            for name in ("train", "test", "validation")
        }

        assert {name: len(frames) for name, frames in parts.items()} == {
            "train": 4, "test": 1, "validation": 1,
        }
        every = sorted(f[0] for frames in parts.values() for f in frames)
        assert every == pytest.approx(
            sorted(energies * EV_PER_KCAL_PER_MOL), rel=1e-5
        )

    def test_the_writing_modes(self, nep_data, monkeypatch):
        run_cli(monkeypatch, "build_nep_traj", "md", "-o", "nep.xyz")

        with pytest.raises(FileWritingModeError):
            run_cli(monkeypatch, "build_nep_traj", "md", "-o", "nep.xyz")


FILE_KEYS = (
    "start_file", "output_file", "info_file", "energy_file", "traj_file",
    "vel_file", "charge_file", "force_file", "restart_file", "file_prefix",
)


def numbered(text):
    """The numbers in a value such as pq-md-08.rst: (prefix, number, rest)."""

    match = re.fullmatch(r"(.*?)(\d+)(\D*)", text)

    return match.group(1), int(match.group(2)), match.group(3)


@pytest.mark.parametrize("example_dir", ["continue_input"], indirect=False)
class TestContinueInput:

    def test_pq_the_next_input_file(self, test_integration_folder, monkeypatch):
        run_cli(monkeypatch, "continue_input", "pq/run-08.in")

        previous = Path("pq/run-08.in").read_text(encoding="utf-8").splitlines()
        following = Path("pq/run-09.in").read_text(encoding="utf-8").splitlines()

        assert len(previous) == len(following)
        changed = [
            (old, new) for old, new in zip(previous, following) if old != new
        ]
        assert changed
        for old, new in changed:
            key = old.split("=")[0].strip()
            assert key in FILE_KEYS
            old_value = old.split("=")[1].strip().rstrip(";")
            new_value = new.split("=")[1].strip().rstrip(";")
            old_parts, new_parts = numbered(old_value), numbered(new_value)
            assert new_parts[0] == old_parts[0]
            assert new_parts[2] == old_parts[2]
            # the counter of the file names goes up by one
            assert new_parts[1] == old_parts[1] + 1 or key == "start_file"

        # the start file of the next run is the restart file of the previous run
        values = {
            line.split("=")[0].strip(): line.split("=")[1].strip().rstrip(";")
            for line in following if "=" in line
        }
        old_values = {
            line.split("=")[0].strip(): line.split("=")[1].strip().rstrip(";")
            for line in previous if "=" in line
        }
        assert values["start_file"] == old_values["restart_file"]
        assert values["restart_file"] == "pq-md-09.rst"

    def test_the_result_equals_the_reference_of_the_example(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(monkeypatch, "continue_input", "pq/run-08.in")

        assert Path("pq/run-09.in").read_text(encoding="utf-8") == Path(
            "pq/run-09.in.ref"
        ).read_text(encoding="utf-8")

    def test_qmcfc_with_its_format(self, test_integration_folder, monkeypatch):
        run_cli(
            monkeypatch, "continue_input", "qmcfc/run-01.in", "--input-format",
            "QMCFC",
        )

        assert Path("qmcfc/run-02.in").read_text(encoding="utf-8") == Path(
            "qmcfc/run-02.in.ref"
        ).read_text(encoding="utf-8")

    def test_several_continuations_in_a_row(
        self, test_integration_folder, monkeypatch
    ):
        run_cli(monkeypatch, "continue_input", "pq/run-08.in", "-n", "3")

        for number in (9, 10, 11):
            assert Path(f"pq/run-{number:02d}.in").is_file()
        text = Path("pq/run-11.in").read_text(encoding="utf-8")

        assert "start_file = pq-md-10.rst;" in text
        assert "restart_file = pq-md-11.rst;" in text
        assert "traj_file = pq-md-11.xyz;" in text
