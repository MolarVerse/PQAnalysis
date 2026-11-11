"""
Tests for the traj2com CLI.
"""

import sys
from unittest import mock

import numpy as np
import pytest

from PQAnalysis.cli import traj2com as traj2com_cli
from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.cli.traj2com import Traj2COMCLI

from . import ArgparseNamespace


def _write_xyz(filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as file:
        print("4 10.0 10.0 10.0 90.0 90.0 90.0", file=file)
        print("", file=file)
        print("C 0.0 0.0 0.0", file=file)
        print("C 2.0 0.0 0.0", file=file)
        print("H 0.0 2.0 0.0", file=file)
        print("H 0.0 4.0 0.0", file=file)


class TestTraj2COMCLI:

    def test_program_name(self):
        assert Traj2COMCLI.program_name() == "traj2com"

    def test_run_passes_the_arguments(self):
        args = ArgparseNamespace(
            trajectory_file=["a.xyz", "b.xyz"],
            output="com.xyz",
            mode="w",
            engine="PQ",
            selection="*|X",
            use_full_atom_info=True,
            group=None,
            by_residue=True,
            restart_file="input.rst",
            moldescriptor_file="moldescriptor.dat",
            continuous=True,
        )

        with mock.patch.object(traj2com_cli, "traj2com") as tool:
            Traj2COMCLI.run(args)

        tool.assert_called_once_with(
            trajectory_files=["a.xyz", "b.xyz"],
            output="com.xyz",
            mode="w",
            md_format="PQ",
            selection="*|X",
            use_full_atom_info=True,
            group=None,
            by_residue=True,
            restart_file="input.rst",
            moldescriptor_file="moldescriptor.dat",
            continuous=True,
        )

    def test_main_parses_the_options(self, monkeypatch):
        called = []
        monkeypatch.setattr(
            traj2com_cli, "traj2com", lambda **kwargs: called.append(kwargs)
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "traj2com", "traj.xyz", "--group", "3", "--selection", "O",
                "--use-full-atom-info", "--continuous", "-o", "com.xyz",
                "--log-file", "off",
            ],
        )

        traj2com_cli.main()

        assert called[0]["trajectory_files"] == ["traj.xyz"]
        assert called[0]["group"] == 3
        assert called[0]["by_residue"] is False
        assert called[0]["continuous"] is True
        assert called[0]["selection"] == "O"
        assert called[0]["use_full_atom_info"] is True
        assert called[0]["output"] == "com.xyz"

    def test_group_and_by_residue_exclude_each_other(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys,
            "argv",
            ["traj2com", "traj.xyz", "--group", "3", "--by-residue", "--log-file", "off"],
        )

        with pytest.raises(SystemExit):
            traj2com_cli.main()

        assert "not allowed with argument" in capsys.readouterr().err

    @pytest.mark.usefixtures("tmpdir")
    def test_pqanalysis_main_dispatches_the_command(self, monkeypatch):
        _write_xyz("traj.xyz")
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "pqanalysis", "--log-file", "off", "traj2com", "traj.xyz",
                "--group", "2", "-o", "com.xyz",
            ],
        )

        pqanalysis_main()

        lines = open("com.xyz", encoding="utf-8").read().splitlines()
        assert lines[0].split()[0] == "2"
        assert np.allclose(
            [float(value) for value in lines[2].split()[1:]], [1.0, 0.0, 0.0],
            atol=1e-6
        )
