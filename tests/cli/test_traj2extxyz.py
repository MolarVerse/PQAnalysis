from unittest import mock

import numpy as np
import pytest

from PQAnalysis.cli.traj2extxyz import Traj2ExtXYZCLI, main, traj2extxyz
from PQAnalysis.cli.main import main as pqanalysis_main
from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import read_trajectory
from PQAnalysis.traj import TrajectoryFormat

from PQAnalysis.exceptions import PQValueError
from PQAnalysis.utils.units import eV, kcal_per_mol

from . import ArgparseNamespace
from ..conftest import assert_logging, assert_logging_with_exception



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


def test_run_passes_extxyz_profile():
    args = ArgparseNamespace(
        trajectory_file=["input.xyz"],
        output="output.extxyz",
        engine="PQ",
        mode="w",
        extxyz_profile="pq",
    )

    with mock.patch("PQAnalysis.cli.traj2extxyz.traj2extxyz") as convert:
        Traj2ExtXYZCLI.run(args)

    convert.assert_called_once_with(
        ["input.xyz"],
        "output.extxyz",
        "w",
        "PQ",
        "pq",
    )


@mock.patch(
    "argparse.ArgumentParser.parse_args",
    return_value=ArgparseNamespace(
        trajectory_file=["input.xyz"],
        output="output.extxyz",
        engine="PQ",
        mode="w",
        extxyz_profile="ase",
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
        extxyz_profile="ase",
        log_file=None,
        logging_level="INFO",
    )
)
def main_pqanalysis_traj2extxyz(mock_args):
    pqanalysis_main()

    _assert_extxyz_output("output.extxyz")


def _write_pq_extxyz(filename: str, units: str | None = "pq") -> None:
    units_key = f" units={units}" if units is not None else ""
    with open(filename, "w", encoding="utf-8") as file:
        print("2", file=file)
        print(
            'Lattice="10 0 0 0 11 0 0 0 12" '
            "Properties=species:S:1:pos:R:3:vel:R:3:forces:R:3:charge:R:1 "
            f'pbc="T T T" step=10 energy=-1.5{units_key}',
            file=file,
        )
        print("H 0.0 0.0 0.0 1.0 2.0 3.0 0.5 0.0 0.0 0.41", file=file)
        print("O 0.0 0.0 1.0 4.0 5.0 6.0 -0.5 0.0 0.0 -0.82", file=file)


@pytest.mark.usefixtures("tmpdir")
def test_traj2extxyz_converts_pq_extxyz_to_ase_units():
    from ase.io import read as ase_read

    _write_pq_extxyz("run.extxyz")

    traj2extxyz(["run.extxyz"], output="ase.extxyz", extxyz_profile="ase")

    with open("ase.extxyz", "r", encoding="utf-8") as file:
        assert "units=ase" in file.read().splitlines()[1]

    frame = ase_read("ase.extxyz", index=0)
    factor = eV / kcal_per_mol

    assert np.isclose(frame.get_potential_energy(), -1.5 * factor)
    assert np.allclose(frame.get_forces()[0], [0.5 * factor, 0.0, 0.0])
    assert np.allclose(frame.arrays["vel"][1], [4.0, 5.0, 6.0])
    assert np.allclose(frame.cell.lengths(), [10.0, 11.0, 12.0])


@pytest.mark.usefixtures("tmpdir")
def test_traj2extxyz_warns_without_units(caplog):
    _write_pq_extxyz("run.extxyz", units=None)

    assert_logging(
        caplog,
        "io.conversion_api",
        "WARNING",
        "run.extxyz has no units key; it is read as PQ units.",
        traj2extxyz,
        ["run.extxyz"],
        output="nep.extxyz",
        extxyz_profile="nep",
    )

    with open("nep.extxyz", "r", encoding="utf-8") as file:
        assert "units=nep" in file.read().splitlines()[1]


@pytest.mark.usefixtures("tmpdir")
def test_traj2extxyz_rejects_converted_files(caplog):
    _write_pq_extxyz("ase.extxyz", units="ase")

    assert_logging_with_exception(
        caplog,
        "io.conversion_api",
        "ERROR",
        (
            "ase.extxyz is in ase units, but only files in PQ units "
            "(units=pq) can be converted."
        ),
        PQValueError,
        traj2extxyz,
        ["ase.extxyz"],
        output="twice.extxyz",
    )
