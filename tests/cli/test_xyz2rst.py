import pytest

from unittest import mock
from filecmp import cmp as filecmp

from . import ArgparseNamespace
from PQAnalysis.cli.xyz2rst import main
from PQAnalysis.io.conversion_api import xyz2rst


@pytest.mark.parametrize("example_dir", ["xyz2rst"], indirect=False)
def test_xyz2rst(test_with_data_dir, capsys):
    print()
    xyz2rst("md-01.xyz", force_file="md-01.force", velocity_file="md-01.vel")

    captured = capsys.readouterr()
    assert captured.out == """
Box  15.0623 15.0964 20.0232  89.9232 90.2261 120.324
C    0    0    1.0 1.0 1.0 1.0 2.0 3.0 1.0 2.0 3.0
H    1    0    2.0 2.0 2.0 2.0 2.0 2.0 4.0 5.0 6.0
N    2    0    3.0 3.0 3.0 3.0 3.0 3.0 7.0 8.0 9.0
N    3    0    4.0 4.0 4.0 4.0 4.0 4.0 10.0 11.0 12.0
"""
    print()
    xyz2rst("md-01_nobox.xyz")

    captured = capsys.readouterr()
    assert captured.out == """
C    0    0    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
H    1    0    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0
N    2    0    3.0 3.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0
N    3    0    4.0 4.0 4.0 0.0 0.0 0.0 0.0 0.0 0.0
"""


@pytest.mark.parametrize("example_dir", ["xyz2rst"], indirect=False)
def test_main(test_with_data_dir, capsys):
    main_xyz2rst()

    captured = capsys.readouterr()
    assert captured.out == """
Box  15.0623 15.0964 20.0232  89.9232 90.2261 120.324
C    0    0    1.0 1.0 1.0 1.0 2.0 3.0 1.0 2.0 3.0
H    1    0    2.0 2.0 2.0 2.0 2.0 2.0 4.0 5.0 6.0
N    2    0    3.0 3.0 3.0 3.0 3.0 3.0 7.0 8.0 9.0
N    3    0    4.0 4.0 4.0 4.0 4.0 4.0 10.0 11.0 12.0
"""

    main_xyz2rst_qmcfc()

    captured = capsys.readouterr()
    assert captured.out == """
Box  15.0623 15.0964 20.0232  89.9232 90.2261 120.324
C    0    0    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
H    1    0    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0
N    2    0    3.0 3.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 3.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0
N    3    0    4.0 4.0 4.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 4.0 4.0 0.0 0.0 0.0 0.0 0.0 0.0
"""


@mock.patch(
    'argparse.ArgumentParser.parse_args',
    return_value=ArgparseNamespace(
        xyz_file="md-01.xyz",
        velocity_file="md-01.vel",
        force_file="md-01.force",
        output=None,
        engine="PQ",
        mode="w",
        log_file=None,
        logging_level="INFO",
    )
)
def main_xyz2rst(mock_args):
    print()
    main()


@mock.patch(
    'argparse.ArgumentParser.parse_args',
    return_value=ArgparseNamespace(
        xyz_file="md-01_qmcfc.xyz",
        velocity_file=None,
        force_file=None,
        output=None,
        engine="qmcfc",
        mode="w",
        log_file=None,
        logging_level="INFO",
    )
)
def main_xyz2rst_qmcfc(mock_args):
    print()
    main()
