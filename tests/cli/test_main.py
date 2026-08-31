"""Tests for lazy unified CLI dispatch."""

import sys

from importlib import import_module
from types import SimpleNamespace

import pytest

import PQAnalysis.cli._argument_parser as argument_parser

from PQAnalysis.cli.main import _COMMANDS, _detect_command



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



@pytest.mark.parametrize(
    ("command", "module_name", "class_name"),
    [
        (command, module_name, class_name)
        for command, (
            module_name,
            class_name,
            _description,
        ) in _COMMANDS.items()
    ],
    ids=_COMMANDS,
)
def test_deferred_command_contract(command, module_name, class_name):
    module = import_module(module_name, "PQAnalysis.cli")
    command_cli = getattr(module, class_name)
    parser = argument_parser._ArgumentParser(
        prog=f"pqanalysis-{command}"
    )

    assert command_cli.program_name() == command
    command_cli.add_arguments(parser)



def test_argcomplete_is_loaded_only_for_shell_completion(monkeypatch):
    calls = []
    fake_argcomplete = SimpleNamespace(
        autocomplete=lambda parser: calls.append(parser)
    )

    monkeypatch.setenv("_ARGCOMPLETE", "1")
    monkeypatch.setitem(sys.modules, "argcomplete", fake_argcomplete)
    monkeypatch.setattr(argument_parser, "print_header", lambda: None)

    parser = argument_parser._ArgumentParser(prog="pqanalysis-test")
    root_logger = argument_parser.logging.getLogger()
    original_level = root_logger.level

    try:
        args = parser.parse_args(["--log-file", "off"])
    finally:
        root_logger.setLevel(original_level)

    assert calls == [parser]
    assert args.progress is True



@pytest.mark.parametrize(
    ("flags", "expected"),
    [
        ([], True),
        (["--progress"], True),
        (["--no-progress"], False),
    ],
)
def test_progress_flag_matches_its_help_text(monkeypatch, flags, expected):
    # --progress used to be a store_false flag, so passing it hid the
    # progress bar although its help text reads "Show progress bar."
    monkeypatch.setattr(argument_parser, "print_header", lambda: None)
    monkeypatch.setattr(
        argument_parser.config,
        "with_progress_bar",
        argument_parser.config.with_progress_bar,
    )

    parser = argument_parser._ArgumentParser(prog="pqanalysis-test")
    root_logger = argument_parser.logging.getLogger()
    original_level = root_logger.level

    try:
        args = parser.parse_args([*flags, "--log-file", "off"])
    finally:
        root_logger.setLevel(original_level)

    assert args.progress is expected
    assert argument_parser.config.with_progress_bar is expected
