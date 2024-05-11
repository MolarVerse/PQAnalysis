import pytest

from filecmp import cmp as filecmp
from unittest import mock

from PQAnalysis.cli.traj2box import traj2box, main

from . import ArgparseNamespace



@pytest.mark.parametrize("example_dir", ["traj2box"], indirect=False)
def test_traj2box(test_with_data_dir):
    traj2box(trajectory_files=["test.xyz"], vmd=False, output="test_box.dat")

    assert filecmp("box.dat", "test_box.dat")

    traj2box(
        trajectory_files=["test.xyz"],
        vmd=True,
        output="test_box.vmd.xyz"
    )

    assert filecmp("box.vmd.xyz", "test_box.vmd.xyz")



@pytest.mark.parametrize("example_dir", ["traj2box"], indirect=False)
def test_main(test_with_data_dir):
    main_box_file()
    main_vmd()



@mock.patch(
    'argparse.ArgumentParser.parse_args',
    return_value=ArgparseNamespace(
    trajectory_file=["test.xyz"],
    vmd=False,
    output="test_box.dat",
    log_file=None,
    logging_level="INFO",
    mode='w'
    )
)
def main_box_file(mock_args):
    main()
    assert filecmp("box.dat", "test_box.dat")



@mock.patch(
    'argparse.ArgumentParser.parse_args',
    return_value=ArgparseNamespace(
    trajectory_file=["test.xyz"],
    vmd=True,
    output="test_box.vmd.xyz",
    log_file=None,
    logging_level="INFO",
    mode='w'
    )
)
def main_vmd(mock_args):
    main()
    assert filecmp("box.vmd.xyz", "test_box.vmd.xyz")
