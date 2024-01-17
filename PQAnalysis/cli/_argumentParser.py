import argparse

import PQAnalysis.config as config
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis._version import __version__


class ArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parse_progress()
        self._parse_version()

    def parse_args(self):
        args = super().parse_args()
        config.with_progress_bar = args.progress
        return args

    def parse_input_file(self):
        super().add_argument('input_file', type=str, help='The input file.')

    def parse_trajectory_file(self):
        super().add_argument('trajectory_file', type=str,
                             nargs='+', help='The trajectory file(s).')

    def parse_output_file(self):
        super().add_argument('-o', '--output', type=str, default=None,
                             help='The output file. If not specified, the output is printed to stdout.')

    def parse_engine_format(self):
        super().add_argument('--engine-format', choices=MDEngineFormat.values(), type=MDEngineFormat, default="pimd-qmcf",
                             help='The case-insensitive MDEngineFormat of the input trajectory.')

    def _parse_progress(self):
        argparse.Namespace(progress=False)
        super().add_argument('--progress', action='store_false', help='Show progress bar.')

    def _parse_version(self):
        argparse.Namespace(version=__version__)
        super().add_argument('--version', action='version', version=__version__)

    @property
    def engine_format(self):
        return MDEngineFormat(self.md_format)
