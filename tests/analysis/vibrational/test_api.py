"""
Tests for the vibrational analysis API.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from PQAnalysis.analysis.vibrational import vibrations
from PQAnalysis.analysis.vibrational import api as vibrational_api
from PQAnalysis.analysis.vibrational.api import (
    _read_atom_charges,
    _read_structure_file,
)
from PQAnalysis.analysis.vibrational.exceptions import VibrationalAnalysisError

from .. import pytestmark  # pylint: disable=unused-import



class TestVibrationalAnalysisAPI:

    """
    Tests for the vibrational analysis API.
    """

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_vibrations(self, test_with_data_dir):
        vibrations("input.in")

        output = Path("wavenumbers.dat")
        normal_modes = Path("normal_modes.dat")

        assert output.is_file()
        assert normal_modes.is_file()
        assert output.read_text(
            encoding="utf-8"
        ).startswith("# Wavenumbers (cm-1)  Intensities (km mol-1)")

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_vibrations_with_mode_output(self, test_with_data_dir):
        vibrations("mode_output.in")

        output = Path("wavenumbers.dat")
        normal_modes = Path("normal_modes.dat")
        modes_file = Path("modes.xyz")
        mode_6 = Path("mode-6.xyz")

        assert output.is_file()
        assert normal_modes.is_file()
        assert modes_file.is_file()
        assert mode_6.is_file()
        assert Path("mode-1.xyz").exists() is False
        assert modes_file.read_text(encoding="utf-8").splitlines(
        )[1].startswith("Properties=species:S:1:pos:R:3:mode:R:3 mode=6")
        assert "frame=1/8" in mode_6.read_text(encoding="utf-8")

    def test_read_single_frame_xyz(self, tmpdir):  # pylint: disable=unused-argument
        Path("structure.xyz").write_text(
            "1\ncomment\nH 0.0 0.0 0.0\n",
            encoding="utf-8",
        )

        system = _read_structure_file("structure.xyz")

        assert system.n_atoms == 1
        assert _read_atom_charges(system, None) is None

    def test_rejects_multi_frame_xyz(self, tmpdir):  # pylint: disable=unused-argument
        Path("structure.xyz").write_text(
            "1\nfirst\nH 0.0 0.0 0.0\n"
            "1\nsecond\nH 0.0 0.0 0.0\n",
            encoding="utf-8",
        )

        with pytest.raises(VibrationalAnalysisError) as exception:
            _read_structure_file("structure.xyz")

        assert str(
            exception.value
        ) == "XYZ structure input must contain exactly one frame."

    def test_read_atom_charges_for_xyz(self, tmpdir):  # pylint: disable=unused-argument
        Path("structure.xyz").write_text(
            "1\ncomment\nH 0.0 0.0 0.0\n",
            encoding="utf-8",
        )
        Path("moldescriptor.dat").write_text(
            "# Molecule 1\nH 1 0.2\nH 1 0.2\n",
            encoding="utf-8",
        )

        system = _read_structure_file("structure.xyz")

        assert _read_atom_charges(system,
                                  "moldescriptor.dat").tolist() == [0.2]

    def test_rejects_xyz_with_multiple_molecule_types(self, tmpdir):
        Path("structure.xyz").write_text(
            "1\ncomment\nH 0.0 0.0 0.0\n",
            encoding="utf-8",
        )
        Path("moldescriptor.dat").write_text(
            "# Molecule 1\nH 1 0.2\nH 1 0.2\n"
            "# Molecule 2\nO 1 -0.4\nO 1 -0.4\n",
            encoding="utf-8",
        )

        system = _read_structure_file("structure.xyz")

        with pytest.raises(VibrationalAnalysisError) as exception:
            _read_atom_charges(system, "moldescriptor.dat")

        assert str(
            exception.value
        ) == "XYZ input requires a moldescriptor file with exactly one molecule type."

    def test_rejects_xyz_moldescriptor_size_mismatch(self, tmpdir):
        Path("structure.xyz").write_text(
            "2\ncomment\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0\n",
            encoding="utf-8",
        )
        Path("moldescriptor.dat").write_text(
            "# Molecule 1\nH 1 0.2\nH 1 0.2\n",
            encoding="utf-8",
        )

        system = _read_structure_file("structure.xyz")

        with pytest.raises(VibrationalAnalysisError) as exception:
            _read_atom_charges(system, "moldescriptor.dat")

        assert str(
            exception.value
        ) == "The moldescriptor molecule size must match the XYZ structure."

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_rejects_missing_restart_residue_charge(
        self, test_with_data_dir, monkeypatch
    ):
        Path("moldescriptor.dat").write_text(
            "# Molecule 2\nH 1 0.2\nH 1 0.2\n",
            encoding="utf-8",
        )

        class FakeMoldescriptorReader:

            def __init__(self, filename):
                self.filename = filename

            def read(self):
                return [
                    SimpleNamespace(
                        id=2, n_atoms=1, partial_charges=np.array([0.2])
                    )
                ]

        monkeypatch.setattr(
            vibrational_api, "MoldescriptorReader", FakeMoldescriptorReader
        )

        system = SimpleNamespace(
            n_atoms=1,
            topology=SimpleNamespace(
                residue_ids=np.array([1]),
                residues=[
                    SimpleNamespace(
                        n_atoms=1, partial_charges=np.array([0.2])
                    )
                ],
                residue_atom_indices=[np.array([0])],
            ),
        )

        with pytest.raises(VibrationalAnalysisError) as exception:
            _read_atom_charges(system, "moldescriptor.dat")

        assert str(
            exception.value
        ) == "Residue id 1 has no moldescriptor entry."

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_rejects_restart_residue_size_mismatch(self, test_with_data_dir):
        Path("moldescriptor.dat").write_text(
            "# Molecule 1\nH 1 0.2\nH 1 0.2\n",
            encoding="utf-8",
        )
        residue = SimpleNamespace(n_atoms=3, partial_charges=np.array([0.2]))
        system = SimpleNamespace(
            n_atoms=1,
            topology=SimpleNamespace(
                residue_ids=np.array([1]),
                residues=[residue],
                residue_atom_indices=[np.array([0])],
            ),
        )

        with pytest.raises(VibrationalAnalysisError) as exception:
            _read_atom_charges(system, "moldescriptor.dat")

        assert str(
            exception.value
        ) == "The moldescriptor residue size does not match the structure."
