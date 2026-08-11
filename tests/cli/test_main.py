"""Tests for lazy unified CLI dispatch."""

import pytest

from PQAnalysis.cli.main import _detect_command



@pytest.mark.parametrize(
    ("arguments", "expected"),
    [
        (["msd", "input.in"], "msd"),
        (["--progress", "check_momentum", "traj.vel"], "check_momentum"),
        (["--logging-level", "DEBUG", "rdf", "input.in"], "rdf"),
        (["--logging-level=INFO", "vacf", "input.in"], "vacf"),
        (["--log-file", "off", "convert", "rdf.dat"], "convert"),
        (["--log-file=run.log", "vibrations", "input.in"], "vibrations"),
        (["--help"], None),
        (["--help", "rdf"], None),
        (["--version", "msd"], None),
        ([], None),
    ],
)
def test_detect_command(arguments, expected):
    assert _detect_command(arguments) == expected
