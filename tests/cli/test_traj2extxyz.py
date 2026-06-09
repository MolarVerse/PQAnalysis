from unittest import mock

import numpy as np
import pytest

from PQAnalysis.cli.traj2extxyz import Traj2ExtXYZCLI, main, traj2extxyz
from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import read_trajectory
from PQAnalysis.traj import TrajectoryFormat

from . import ArgparseNamespace



def _write_pq_xyz(filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as file:
        print("2 10.0 11.0 12.0 90.0 90.0 90.0", file=file)
        print("", file=file)
        print("H 0.0 0.0 0.0", file=file)
        print("O 0.0 0.0 1.0", file=file)


def _assert_extxyz_output(filename: str) -> None:
    from ase.io import read as ase_read

    with open(filename, "r", encoding="utf-8") as file:
        lines = file.read().splitlines()

    assert lines[0] == "2"
    assert "Lattice=" in lines[1]
    assert 'Properties="species:S:1:pos:R:3"' in lines[1]

    trajectory = read_trajectory(filename, traj_format=TrajectoryFormat.EXTXYZ)
    frame = trajectory[0]

    assert frame.atoms == [Atom("h"), Atom("o")]
    assert np.allclose(frame.pos, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert frame.cell == Cell(10.0, 11.0, 12.0)

    ase_frame = ase_read(filename, index=0)
    assert ase_frame.get_chemical_symbols() == ["H", "O"]
    assert np.allclose(ase_frame.cell.lengths(), [10.0, 11.0, 12.0])


@pytest.mark.usefixtures("tmpdir")
def test_traj2extxyz():
    _write_pq_xyz("input.xyz")

    traj2extxyz(["input.xyz"], output="output.extxyz")

    _assert_extxyz_output("output.extxyz")


@pytest.mark.usefixtures("tmpdir")
def test_main():
    _write_pq_xyz("input.xyz")

    main_traj2extxyz()


@pytest.mark.usefixtures("tmpdir")
def test_pqanalysis_main():
    _write_pq_xyz("input.xyz")

    main_pqanalysis_traj2extxyz()


@mock.patch(
    "argparse.ArgumentParser.parse_args",
    return_value=ArgparseNamespace(
        trajectory_file=["input.xyz"],
        output="output.extxyz",
        engine="PQ",
        mode="w",
        log_file=None,
        logging_level="INFO",
    )
)
def main_traj2extxyz(mock_args):
    assert Traj2ExtXYZCLI.program_name() == "traj2extxyz"

    main()

    _assert_extxyz_output("output.extxyz")


@mock.patch(
    "argparse.ArgumentParser.parse_args",
    return_value=ArgparseNamespace(
        cli_command="traj2extxyz",
        trajectory_file=["input.xyz"],
        output="output.extxyz",
        engine="PQ",
        mode="w",
        log_file=None,
        logging_level="INFO",
    )
)
def main_pqanalysis_traj2extxyz(mock_args):
    pqanalysis_main()

    _assert_extxyz_output("output.extxyz")
