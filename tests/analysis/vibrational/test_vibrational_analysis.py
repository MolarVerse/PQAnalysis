"""
Tests for vibrational analysis numerical routines.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.vibrational.vibrational_analysis import (
    calculate,
    mass_weighted_hessian,
    read_hessian_file,
    select_mode_indices,
    write_extxyz_modes,
    write_xyz_modes,
)
from PQAnalysis.io import RestartFileReader

from .. import pytestmark  # pylint: disable=unused-import



class TestVibrationalAnalysis:

    """
    Tests for the vibrational analysis routines.
    """

    def test_mass_weighted_hessian_preserves_sign(self):
        hessian = np.diag([-4.0, 2.0, 8.0])

        assert np.allclose(
            np.linalg.eigvalsh(
                mass_weighted_hessian(hessian, np.array([1.0]))
            ),
            np.array([-4.0, 2.0, 8.0]),
        )
        assert np.allclose(
            np.linalg.eigvalsh(
                mass_weighted_hessian(hessian, np.array([1.0]), sign=-1.0)
            ),
            np.array([-8.0, -2.0, 4.0]),
        )

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_calculate_h2o(self, test_with_data_dir):
        system = RestartFileReader("h2o.rst").read()
        hessian = read_hessian_file("hessian.dat")

        result = calculate(system.atomic_masses, system.pos, hessian)

        assert result.wavenumbers.shape == (9, )
        assert result.force_constants.shape == (9, )
        assert result.reduced_masses.shape == (9, )
        assert result.normal_modes.shape == (9, 9)
        assert result.intensities is None
        assert np.all(np.isfinite(result.wavenumbers))

    def test_select_mode_indices(self):
        wavenumbers = np.array([-2.0, 0.0, 0.1, 10.0])

        assert select_mode_indices(wavenumbers) == [0, 1, 2, 3]
        assert select_mode_indices(wavenumbers, "nonzero") == [0, 2, 3]
        assert select_mode_indices(wavenumbers, "positive",
                                   threshold=1.0) == [3]
        assert select_mode_indices(wavenumbers, [1, 4, 4]) == [0, 3]

    def test_write_xyz_modes(self, tmpdir):  # pylint: disable=unused-argument
        normal_modes = np.eye(6)
        atom_coords = np.zeros((2, 3))
        atom_names = ["o", "h"]
        wavenumbers = np.array([0.0, 50.0, 100.0, 150.0, 200.0, 250.0])

        write_xyz_modes(
            normal_modes,
            atom_coords,
            atom_names,
            filename="mode",
            wavenumbers=wavenumbers,
            modes=[2],
            n_frames=4,
            amplitude=0.5,
        )

        mode_file = "mode-2.xyz"
        with open(mode_file, "r", encoding="utf-8") as file:
            lines = file.readlines()

        assert len(lines) == 16
        assert lines[1].startswith(
            "mode=2 frequency_cm-1=5.00000000e+01 frame=1/4"
        )
        assert lines[5].startswith(
            "mode=2 frequency_cm-1=5.00000000e+01 frame=2/4"
        )
        assert lines[6].split()[0] == "O"
        assert lines[7].split()[0] == "H"
        assert lines[6].split()[2] == "0.5"

        write_xyz_modes(
            normal_modes,
            atom_coords,
            atom_names,
            "legacy",
            0.5,
            0.5,
            wavenumbers=wavenumbers,
            modes=[2],
        )

        with open("legacy-2.xyz", "r", encoding="utf-8") as file:
            lines = file.readlines()

        assert len(lines) == 12

    def test_write_extxyz_modes(self, tmpdir):  # pylint: disable=unused-argument
        normal_modes = np.eye(6)
        atom_coords = np.zeros((2, 3))
        atom_names = ["o", "h"]
        wavenumbers = np.array([0.0, 50.0, 100.0, 150.0, 200.0, 250.0])
        intensities = np.arange(6, dtype=float)

        write_extxyz_modes(
            normal_modes,
            atom_coords,
            atom_names,
            filename="modes.xyz",
            wavenumbers=wavenumbers,
            intensities=intensities,
            modes=[6],
        )

        with open("modes.xyz", "r", encoding="utf-8") as file:
            lines = file.readlines()

        assert len(lines) == 4
        assert lines[1].startswith(
            "Properties=species:S:1:pos:R:3:mode:R:3 "
            "mode=6 frequency_cm-1=2.50000000e+02"
        )
        assert "IR_intensity=5.00000000e+00" in lines[1]
        assert lines[2].split()[0] == "O"
        assert lines[3].split()[0] == "H"
        assert len(lines[2].split()) == 7
