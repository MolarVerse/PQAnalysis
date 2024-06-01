"""
This module provides a Wrapper class for the argparse.ArgumentParser.

The _ArgumentParser class is a subclass of the argparse.ArgumentParser class.
It provides a set of methods that are used in many cli scripts of this package.
Therefore, in order to be able to change the behavior of all scripts at once,
the most common methods are implemented here.
"""

import logging
import argparse
import argcomplete

from beartype.typing import Sequence
from rich_argparse import ArgumentDefaultsRichHelpFormatter

import PQAnalysis.config as config  # pylint: disable=consider-using-from-import # here needed to set the config attributes

from PQAnalysis.utils.common import __header__
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis._version import __version__



class _ArgumentParser(argparse.ArgumentParser):

    """
    The _ArgumentParser class is a subclass of the argparse.ArgumentParser class.
    It provides a set of methods that are used in many cli scripts of this package.
    Therefore, in order to be able to change the behavior of all scripts at once,
    the most common methods are implemented here.

    Furthermore, it provides two special methods to parse the progress and version
    arguments. These arguments are used in all scripts and therefore implemented here.
    """

    def __init__(self, *args, **kwargs):
        """
        The initialization of the ArgumentParser class is the same as the initialization
        of the argparse.ArgumentParser class. The only difference is that the progress
        and version arguments are parsed here automatically.

        Parameters
        ----------
        args : list
            The list of positional arguments for the argparse.ArgumentParser class.
        kwargs : dict
            The dictionary of keyword arguments for the argparse.ArgumentParser class.
        """

        super().__init__(
            formatter_class=ArgumentDefaultsRichHelpFormatter, *args, **kwargs
        )

        # To remove the ".py" ending from the prog argument
        if 'prog' not in kwargs:
            kwargs['prog'] = self.prog.split(".")[0]
            super().__init__(
                formatter_class=ArgumentDefaultsRichHelpFormatter,
                *args,
                **kwargs
            )

    def parse_args(
        self,
        args: Sequence[str] | None = None,
        namespace: None = None
    ) -> argparse.Namespace:
        """
        The parse_args method is the same as the parse_args method of the
        argparse.ArgumentParser class. The only difference is that the progress
        argument is parsed here automatically and the config.with_progress_bar
        attribute is set accordingly.

        Returns
        -------
        argparse.Namespace
            The parsed arguments.
        """
        self._parse_progress()
        self._parse_version()
        self._parse_log_file()
        self._parse_logging_level()

        argcomplete.autocomplete(self)
        args = super().parse_args(args, namespace)

        config.with_progress_bar = args.progress

        if args.log_file is None:
            args.log_file = "__LOG_DEFINED_BY_TIMESTAMP__"

        if args.log_file.lower() != "off":
            config.use_log_file = True
            if args.log_file != "__LOG_DEFINED_BY_TIMESTAMP__":
                config.log_file_name = args.log_file

        logger = logging.getLogger()
        logger.setLevel(logging.getLevelName(args.logging_level))

        return args

    def parse_input_file(self):
        """
        The parse_input_file method adds the input_file argument to the parser.
        The input_file argument is a positional argument and is required.
        """
        super().add_argument('input_file', type=str, help='The input file.')

    def parse_trajectory_file(self):
        """
        The parse_trajectory_file method adds the trajectory_file argument to the parser.
        The trajectory_file argument is a positional argument and is required.
        """
        super().add_argument(
            'trajectory_file',
            type=str,
            nargs='+',
            help='The trajectory file(s).'
        )

    def parse_output_file(self):
        """
        The parse_output_file method adds the output_file argument to the parser.
        The output_file argument is an optional argument and defaults to None.
        """
        super().add_argument(
            '-o',
            '--output',
            type=str,
            default=None,
            help=
            'The output file. If not specified, the output is printed to stdout.'
        )

    def parse_engine(self):
        """
        The parse_engine method adds the engine argument to the parser.
        The engine argument is an optional argument and defaults to PQ.
        """
        super().add_argument(
            '--engine',
            choices=MDEngineFormat.values(),
            type=MDEngineFormat,
            default=MDEngineFormat("PQ"),
            help='The case-insensitive MDEngineFormat of the input trajectory.'
        )

    def parse_mode(self):
        """
        The parse_mode method adds the mode argument to the parser.
        The mode argument is an optional argument and defaults to "w".
        Possible values are "w", "a" and "o".
        """
        super().add_argument(
            '--mode',
            choices=FileWritingMode.values(),
            type=FileWritingMode,
            default=FileWritingMode("w"),
            help=(
                'The writing mode. The following modes '
                'are available: "w": write, "a": append, '
                '"o": overwrite.'
            )
        )

    def _parse_progress(self):
        """
        The _parse_progress method adds the progress argument to the parser.
        The progress argument is an optional argument and defaults to True.
        """
        super().add_argument(
            '--progress', action='store_false', help='Show progress bar.'
        )

    def _parse_version(self):
        """
        The _parse_version method adds the version argument to the parser.
        The version argument is an optional argument and defaults to the current version.
        """
        super().add_argument(
            '--version', action='version', version=__version__
        )

    def _parse_log_file(self):
        """
        The _parse_log_file method adds the log_file argument to the parser.
        The log_file argument is an optional argument and defaults to None.
        If only the option log_file is given without a value the log will
        be printed to an automatically generated file.
        """
        super().add_argument(
            '--log-file',
            type=str,
            default=None,
            const="__LOG_DEFINED_BY_TIMESTAMP__",
            nargs='?',
            help=(
                "This options can be used to set wether "
                "a log file should be used or not. "
                "If the option is not given or without a value, "
                "the log will be printed to an "
                "automatically generated file. If the option "
                "is given with the case-insensitive value \"off\""
                ", the log will be only printed "
                "to stderr. If the option is given with a "
                "value, the log will be printed to the given file."
            )
        )

    def _parse_logging_level(self):
        """
        The _parse_logging_level method adds the logging_level argument to the parser.
        The logging_level argument is an optional argument and defaults to "INFO".
        """
        super().add_argument(
            '--logging-level',
            type=str,
            choices=logging.getLevelNamesMapping().keys(),
            default="INFO",
            help='The logging level.'
        )
