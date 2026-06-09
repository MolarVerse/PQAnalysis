import pytest

from unittest import mock
from filecmp import cmp as filecmp

from . import ArgparseNamespace
from PQAnalysis.cli.xyz2rst import main, XYZ2RstCLI
from PQAnalysis.io.conversion_api import xyz2rst
from PQAnalysis.io.restart_file.exceptions import RestartFileWriterError


def test_name():
    assert XYZ2RstCLI.program_name() == "xyz2rst"


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
C    0    0    1.0 1.0 1.0
H    1    0    2.0 2.0 2.0
N    2    0    3.0 3.0 3.0
N    3    0    4.0 4.0 4.0
"""

    print()
    xyz2rst("md-01_null.xyz", random_seed=1, randomize=0.1)

    captured = capsys.readouterr()
    assert captured.out == """
C    0    0    0.16243453323841095 -0.0611756406724453 -0.05281717702746391
H    1    0    -0.10729686170816422 0.08654076606035233 -0.2301538735628128
N    2    0    0.17448118329048157 -0.07612068951129913 0.031903911381959915
N    3    0    -0.024937037378549576 0.14621078968048096 -0.20601406693458557
"""

    print()
    xyz2rst(
        "md-01.xyz",
        force_file="md-01.force",
        velocity_file="md-01.vel",
        moldescriptor_file="moldescriptor.dat",
    )

    captured = capsys.readouterr()
    assert captured.out == """
Box  15.0623 15.0964 20.0232  89.9232 90.2261 120.324
C    1    1    1.0 1.0 1.0 1.0 2.0 3.0 1.0 2.0 3.0
H    2    1    2.0 2.0 2.0 2.0 2.0 2.0 4.0 5.0 6.0
N    1    2    3.0 3.0 3.0 3.0 3.0 3.0 7.0 8.0 9.0
N    2    2    4.0 4.0 4.0 4.0 4.0 4.0 10.0 11.0 12.0
"""

    with open("repeated.xyz", "w", encoding="utf-8") as file:
        print("4", file=file)
        print("# comment", file=file)
        print("C 0.0 0.0 0.0", file=file)
        print("H 0.0 0.0 1.0", file=file)
        print("C 1.0 0.0 0.0", file=file)
        print("H 1.0 0.0 1.0", file=file)

    print()
    xyz2rst("repeated.xyz", moldescriptor_file="moldescriptor.dat")

    captured = capsys.readouterr()
    assert captured.out == """
C    1    1    0.0 0.0 0.0
H    2    1    0.0 0.0 1.0
C    1    1    1.0 0.0 0.0
H    2    1    1.0 0.0 1.0
"""

    with pytest.raises(RestartFileWriterError) as exception:
        xyz2rst("md-01_qmcfc.xyz", moldescriptor_file="moldescriptor.dat")
    assert str(exception.value) == (
        "Could not infer moltype from moldescriptor at atom index 0."
    )

    with open("truncated.xyz", "w", encoding="utf-8") as file:
        print("3", file=file)
        print("# comment", file=file)
        print("C 0.0 0.0 0.0", file=file)
        print("H 0.0 0.0 1.0", file=file)
        print("N 0.0 0.0 2.0", file=file)

    with pytest.raises(RestartFileWriterError) as exception:
        xyz2rst("truncated.xyz", moldescriptor_file="moldescriptor.dat")
    assert str(exception.value) == (
        "Could not infer moltype from moldescriptor at atom index 2."
    )


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

    main_xyz2rst_moldescriptor()

    captured = capsys.readouterr()
    assert captured.out == """
Box  15.0623 15.0964 20.0232  89.9232 90.2261 120.324
C    1    1    1.0 1.0 1.0
H    2    1    2.0 2.0 2.0
N    1    2    3.0 3.0 3.0
N    2    2    4.0 4.0 4.0
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
        moldescriptor_file=None,
        randomize=0.0,
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
        xyz_file="md-01.xyz",
        velocity_file=None,
        force_file=None,
        moldescriptor_file="moldescriptor.dat",
        randomize=0.0,
        output=None,
        engine="PQ",
        mode="w",
        log_file=None,
        logging_level="INFO",
    )
)
def main_xyz2rst_moldescriptor(mock_args):
    print()
    main()


@mock.patch(
    'argparse.ArgumentParser.parse_args',
    return_value=ArgparseNamespace(
        xyz_file="md-01_qmcfc.xyz",
        velocity_file=None,
        force_file=None,
        moldescriptor_file=None,
        randomize=0.0,
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
