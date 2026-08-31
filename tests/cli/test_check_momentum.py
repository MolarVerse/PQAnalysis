"""
Tests for the check_momentum CLI.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.cli import check_momentum as check_momentum_cli
from PQAnalysis.cli.check_momentum import CheckMomentumCLI
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.traj import MDEngineFormat

from ..analysis.momentum.test_momentum import GAS3_LEGACY_NORMS



class TestCheckMomentumCLI:

    """
    Tests for CheckMomentumCLI.
    """

    def test_program_name(self):
        """
        The program name is check_momentum.
        """
        assert CheckMomentumCLI.program_name() == "check_momentum"

    def test_main_dispatches_defaults(self, monkeypatch):
        """
        main() dispatches the positional trajectory files and the
        default settings to the check_momentum API function.
        """
        called = []

        monkeypatch.setattr(
            check_momentum_cli,
            "check_momentum",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["check_momentum", "run-01.vel", "--log-file", "off"],
        )

        check_momentum_cli.main()

        assert called == [
            {
                "trajectory_files": ["run-01.vel"],
                "output": None,
                "selection": None,
                "use_full_atom_info": False,
                "scale": 1e-15,
                "md_format": MDEngineFormat("PQ"),
                "mode": FileWritingMode("w"),
            }
        ]

    def test_main_dispatches_options(self, monkeypatch):
        """
        main() dispatches multiple trajectory files, output, selection,
        scale and engine options.
        """
        called = []

        monkeypatch.setattr(
            check_momentum_cli,
            "check_momentum",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "check_momentum",
                "run-01.vel",
                "run-02.vel",
                "-o",
                "momentum.dat",
                "--selection",
                "O",
                "--use-full-atom-info",
                "--scale",
                "1.0",
                "--engine",
                "qmcfc",
                "--log-file",
                "off",
            ],
        )

        check_momentum_cli.main()

        assert called == [
            {
                "trajectory_files": ["run-01.vel", "run-02.vel"],
                "output": "momentum.dat",
                "selection": "O",
                "use_full_atom_info": True,
                "scale": 1.0,
                "md_format": MDEngineFormat("QMCFC"),
                "mode": FileWritingMode("w"),
            }
        ]

    def test_main_dispatches_repeated_exports(self, monkeypatch):
        called = []
        monkeypatch.setattr(
            check_momentum_cli,
            "check_momentum",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "check_momentum",
                "run-01.vel",
                "--export",
                "momentum.tsv",
                "--export",
                "momentum.xvg",
                "--log-file",
                "off",
            ],
        )

        check_momentum_cli.main()

        assert called[0]["export_files"] == ["momentum.tsv", "momentum.xvg"]

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_main_writes_bitwise_legacy_output(
        self,
        test_with_data_dir,
        monkeypatch,
    ):
        """The real CLI writes float64 values that round-trip exactly."""
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "check_momentum",
                "gas3.vel",
                "--output",
                "momentum.dat",
                "--engine",
                "qmcfc",
                "--log-file",
                "off",
            ],
        )

        check_momentum_cli.main()

        stored = np.loadtxt("momentum.dat", comments="#", usecols=1)
        assert np.array_equal(stored, GAS3_LEGACY_NORMS)
