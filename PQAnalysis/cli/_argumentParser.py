"""
This module provides a Wrapper class for the argparse.ArgumentParser.

The _ArgumentParser class is a subclass of the argparse.ArgumentParser class.
It provides a set of methods that are used in many cli scripts of this package.
Therefore, in order to be able to change the behavior of all scripts at once,
the most common methods are implemented here.
"""

import argparse
import argcomplete

import PQAnalysis.config as config
from PQAnalysis.traj import MDEngineFormat
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
        """

        # To remove the ".py" ending from the prog argument
        super().__init__(*args, **kwargs)

        if 'prog' not in kwargs:
            kwargs['prog'] = self.prog.split(".")[0]
            super().__init__(*args, **kwargs)

        self._parse_progress()
        self._parse_version()

    def parse_args(self) -> argparse.Namespace:
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
        argcomplete.autocomplete(self)
        args = super().parse_args()
        config.with_progress_bar = args.progress
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
        super().add_argument('trajectory_file', type=str,
                             nargs='+', help='The trajectory file(s).')

    def parse_output_file(self):
        """
        The parse_output_file method adds the output_file argument to the parser.
        The output_file argument is an optional argument and defaults to None.
        """
        super().add_argument('-o', '--output', type=str, default=None,
                             help='The output file. If not specified, the output is printed to stdout.')

    def parse_engine(self):
        """
        The parse_engine method adds the engine argument to the parser.
        The engine argument is an optional argument and defaults to pimd-qmcf.
        """
        super().add_argument('--engine', choices=MDEngineFormat.values(), type=MDEngineFormat, default=MDEngineFormat("pimd-qmcf"),
                             help='The case-insensitive MDEngineFormat of the input trajectory.')

    def _parse_progress(self):
        """
        The _parse_progress method adds the progress argument to the parser.
        The progress argument is an optional argument and defaults to True.
        """
        super().add_argument('--progress', action='store_false', help='Show progress bar.')

    def _parse_version(self):
        """
        The _parse_version method adds the version argument to the parser.
        The version argument is an optional argument and defaults to the current version.
        """
        super().add_argument('--version', action='version', version=__version__)
