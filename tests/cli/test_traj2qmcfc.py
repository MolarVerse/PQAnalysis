import pytest
import sys
import argparse

from filecmp import cmp as filecmp
from unittest import mock

from PQAnalysis.cli.traj2qmcfc import traj2qmcfc, main


@pytest.mark.parametrize("example_dir", ["traj2qmcfc"], indirect=False)
def test_traj2qmcfc(test_with_data_dir):
    traj2qmcfc(trajectory_files=[
        "acof_triclinic.xyz", "acof_triclinic_2.xyz"], output="test_traj.qmcfc.xyz")

    assert filecmp("traj.qmcfc.xyz", "test_traj.qmcfc.xyz")


@pytest.mark.parametrize("example_dir", ["traj2qmcfc"], indirect=False)
def test_main(test_with_data_dir):
    main_traj2qmcfc()


@mock.patch('argparse.ArgumentParser.parse_args',
            return_value=argparse.Namespace(trajectory_file=["acof_triclinic.xyz",
                                                             "acof_triclinic_2.xyz"], vmd=False, output="test_traj.qmcfc.xyz"))
def main_traj2qmcfc(mock_args):
    main()
    assert filecmp("traj.qmcfc.xyz", "test_traj.qmcfc.xyz")
