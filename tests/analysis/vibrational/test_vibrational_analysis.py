"""
Tests for vibrational analysis numerical routines.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.vibrational.vibrational_analysis import (
    calculate,
    hessian_sign_factor,
    mode_displacement,
    mass_weighted_hessian,
    read_hessian_file,
    select_mode_indices,
    wavenumber,
    write_calculate_output,
    write_extxyz_modes,
    write_normal_modes,
    write_xyz_modes,
    _mode_frame_count,
    _mode_wavenumbers,
    _scale_mode_to_amplitude,
    _xyz_atom_symbol,
)
from PQAnalysis.analysis.vibrational.exceptions import VibrationalAnalysisError
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

    def test_read_hessian_file_errors(self, tmpdir):  # pylint: disable=unused-argument
        with pytest.raises(VibrationalAnalysisError) as exception:
            read_hessian_file("missing.dat")

        assert str(exception.value) == "Hessian file 'missing.dat' not found."

        with open("bad.dat", "w", encoding="utf-8") as file:
            file.write("1.0 nope\n")

        with pytest.raises(VibrationalAnalysisError) as exception:
            read_hessian_file("bad.dat")

        assert str(
            exception.value
        ) == "Hessian file 'bad.dat' contains non-numeric data."

        with open("rectangular.dat", "w", encoding="utf-8") as file:
            file.write("1.0 2.0 3.0\n4.0 5.0 6.0\n")

        with pytest.raises(VibrationalAnalysisError) as exception:
            read_hessian_file("rectangular.dat")

        assert "Hessian matrix must be square" in str(exception.value)

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

    def test_calculate_validates_shapes(self):
        hessian = np.eye(6)

        with pytest.raises(VibrationalAnalysisError) as exception:
            calculate(np.ones((2, 1)), np.zeros((2, 3)), hessian)

        assert str(
            exception.value
        ) == "Atom masses must be a one-dimensional array."

        with pytest.raises(VibrationalAnalysisError) as exception:
            calculate(np.ones(2), np.zeros((2, 2)), hessian)

        assert str(
            exception.value
        ) == "Atom coordinates must have shape (n_atoms, 3)."

        with pytest.raises(VibrationalAnalysisError) as exception:
            calculate(np.ones(2), np.zeros((2, 3)), np.eye(3))

        assert "Hessian shape must be" in str(exception.value)

        with pytest.raises(VibrationalAnalysisError) as exception:
            calculate(
                np.ones(2),
                np.zeros((2, 3)),
                hessian,
                atom_charges=np.ones(3),
            )

        assert str(
            exception.value
        ) == "The number of atom charges must match the number of atoms."

    def test_hessian_sign_factor_options(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        masses = np.ones(2)
        hessian = np.eye(6)

        assert hessian_sign_factor(coords, masses, hessian, 1.0) == 1.0
        assert hessian_sign_factor(coords, masses, hessian, -1.0) == -1.0
        assert hessian_sign_factor(coords, masses, hessian, "positive") == 1.0
        assert hessian_sign_factor(coords, masses, hessian, "negative") == -1.0
        assert hessian_sign_factor(coords, masses, hessian, "auto") == 1.0
        assert hessian_sign_factor(coords, masses, -hessian, "auto") == -1.0

        one_atom_coords = np.zeros((1, 3))
        assert hessian_sign_factor(
            one_atom_coords,
            np.ones(1),
            np.eye(3),
            "auto",
        ) == 1.0
        nonlinear_coords = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        )
        assert hessian_sign_factor(
            nonlinear_coords,
            np.ones(3),
            np.zeros((9, 9)),
            "auto",
        ) == 1.0

        with pytest.raises(VibrationalAnalysisError) as exception:
            hessian_sign_factor(coords, masses, hessian, "sideways")

        assert str(
            exception.value
        ) == "hessian_sign must be auto, positive, negative, 1, or -1."

    def test_wavenumber_units(self):
        eigenvalues = np.array([1.0, -1.0])

        for unit in ("kcal", "hartree", "ev"):
            values, omega = wavenumber(eigenvalues, unit=unit)
            assert values[0] > 0.0
            assert omega[1] < 0.0

        with pytest.raises(VibrationalAnalysisError) as exception:
            wavenumber(eigenvalues, unit="kj")

        assert str(
            exception.value
        ) == "Invalid unit. Options are kcal, hartree and ev."

    def test_select_mode_indices(self):
        wavenumbers = np.array([-2.0, 0.0, 0.1, 10.0])

        assert select_mode_indices(wavenumbers) == [0, 1, 2, 3]
        assert select_mode_indices(wavenumbers, "nonzero") == [0, 2, 3]
        assert select_mode_indices(wavenumbers, "positive",
                                   threshold=1.0) == [3]
        assert select_mode_indices(wavenumbers, [1, 4, 4]) == [0, 3]

        with pytest.raises(VibrationalAnalysisError) as exception:
            select_mode_indices(wavenumbers, "invalid")

        assert str(
            exception.value
        ) == "modes must be all, nonzero, positive or one-based mode numbers."

        with pytest.raises(VibrationalAnalysisError) as exception:
            select_mode_indices(wavenumbers, [5])

        assert str(
            exception.value
        ) == "Mode 5 is outside the available range 1..4."

    def test_write_calculate_output_without_intensities(self, capsys):
        result = calculate(np.ones(1), np.zeros((1, 3)), np.eye(3))

        write_calculate_output(result)

        assert capsys.readouterr(
        ).out.startswith("# Wavenumbers (cm-1)  Force constants")

    def test_write_normal_modes_stdout(self, capsys):
        write_normal_modes(np.eye(2))

        assert capsys.readouterr().out.splitlines()[0] == "1.0 0.0"

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

        write_xyz_modes(normal_modes, atom_coords, atom_names, filename="zero")

        assert "frequency_cm-1=0.00000000e+00" in open(
            "zero-1.xyz",
            "r",
            encoding="utf-8",
        ).read()

    def test_mode_writer_validation(self, tmpdir):  # pylint: disable=unused-argument
        normal_modes = np.eye(3)
        atom_coords = np.zeros((1, 3))
        atom_names = ["x"]

        with pytest.raises(VibrationalAnalysisError) as exception:
            write_xyz_modes(normal_modes, atom_coords, atom_names, n_frames=0)

        assert str(
            exception.value
        ) == "Number of mode frames must be positive."

        with pytest.raises(VibrationalAnalysisError) as exception:
            write_xyz_modes(
                normal_modes, atom_coords, atom_names, amplitude=-1.0
            )

        assert str(exception.value) == "Mode amplitude must not be negative."

        with pytest.raises(VibrationalAnalysisError) as exception:
            write_xyz_modes(
                normal_modes,
                atom_coords,
                atom_names,
                temperature=0.0,
            )

        assert str(exception.value) == "Mode temperature must be positive."

        with pytest.raises(VibrationalAnalysisError) as exception:
            write_xyz_modes(
                normal_modes, atom_coords, atom_names, threshold=-1.0
            )

        assert str(exception.value) == "Mode threshold must not be negative."

        with pytest.raises(VibrationalAnalysisError) as exception:
            _mode_frame_count(None, 0.25, 0.0)

        assert str(exception.value) == "Mode step must be positive."

        with pytest.raises(VibrationalAnalysisError) as exception:
            _mode_wavenumbers(np.ones(2), normal_modes)

        assert str(
            exception.value
        ) == "The number of wavenumbers must match the number of normal modes."

    def test_mode_scaling(self):
        mode = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]])

        assert np.allclose(_scale_mode_to_amplitude(mode, 0.0), 0.0)
        assert np.allclose(
            _scale_mode_to_amplitude(np.zeros((1, 3)), 0.3), 0.0
        )
        assert np.isclose(
            np.linalg.norm(_scale_mode_to_amplitude(mode, 0.5)[0]), 0.5
        )
        with pytest.raises(VibrationalAnalysisError) as exception:
            _scale_mode_to_amplitude(mode, -1.0)

        assert str(exception.value) == "Mode amplitude must not be negative."

        scaled = mode_displacement(mode, 100.0, temperature=300.0)
        assert scaled[0, 0] > 3.0

        with pytest.raises(VibrationalAnalysisError) as exception:
            mode_displacement(mode, 100.0, temperature=0.0)

        assert str(exception.value) == "Mode temperature must be positive."

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

        write_extxyz_modes(
            normal_modes,
            atom_coords,
            atom_names,
            filename="modes_no_intensity.xyz",
            modes=[1],
        )

        with open("modes_no_intensity.xyz", "r", encoding="utf-8") as file:
            assert "IR_intensity" not in file.readlines()[1]

        with pytest.raises(VibrationalAnalysisError) as exception:
            write_extxyz_modes(
                normal_modes,
                atom_coords,
                atom_names,
                threshold=-1.0,
            )

        assert str(exception.value) == "Mode threshold must not be negative."

    def test_xyz_atom_symbol(self):
        assert _xyz_atom_symbol(None) == "X"
        assert _xyz_atom_symbol("") == "X"
        assert _xyz_atom_symbol("cl") == "Cl"
