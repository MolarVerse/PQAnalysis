"""
Tests for the add_molecules CLI argument parser.
"""

import pytest

import PQAnalysis.cli._argument_parser as argument_parser

from PQAnalysis.cli.add_molecules import AddMoleculesCLI



class TestAddMoleculesParser:

    """
    Tests for the add_molecules argument parser.
    """

    @pytest.fixture
    def parser(self, monkeypatch):
        """
        An add_molecules parser with quiet parse_args side effects.
        """
        monkeypatch.setattr(argument_parser, "print_header", lambda: None)

        parser = argument_parser._ArgumentParser(prog="add_molecules")
        AddMoleculesCLI.add_arguments(parser)

        return parser

    def _parse(self, parser, arguments):
        root_logger = argument_parser.logging.getLogger()
        original_level = root_logger.level

        try:
            return parser.parse_args(
                ["md-01.rst", "mol.xyz", "--log-file", "off"] + arguments
            )
        finally:
            root_logger.setLevel(original_level)

    def test_n_molecules_option_strings(self, parser):
        """
        The n_molecules argument registers both aliases separately.
        """
        option_strings = [
            action.option_strings
            for action in parser._actions
            if action.dest == "n_molecules"
        ]

        assert option_strings == [["-n", "--n-molecules"]]

    def test_n_molecules_long_flag(self, parser):
        """
        The advertised --n-molecules long flag is accepted.
        """
        args = self._parse(parser, ["--n-molecules", "2"])

        assert args.n_molecules == 2

    def test_n_molecules_short_flag(self, parser):
        """
        The -n short flag is accepted.
        """
        args = self._parse(parser, ["-n", "3"])

        assert args.n_molecules == 3

    def test_n_molecules_default(self, parser):
        """
        The n_molecules argument defaults to 1.
        """
        args = self._parse(parser, [])

        assert args.n_molecules == 1
