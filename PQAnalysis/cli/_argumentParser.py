import argparse

from PQAnalysis.traj import MDEngineFormat


class ArgumentParser(argparse.ArgumentParser):
    def __init__(self, description: str) -> None:
        super().__init__(description=description)

    def add_argument(self, *args, **kwargs):
        super().add_argument(*args, **kwargs)

    def parse_arguments(self):
        super().parse_args()

    def parse_output_file(self):
        super().add_argument('-o', '--output', type=str, default=None,
                             help='The output file. If not specified, the output is printed to stdout.')

    def parse_input_file(self):
        super().add_argument('input_file', type=str, help='The input file.')

    def parse_md_format(self):
        super().add_argument('--md-format', type=str, default="pimd-qmcf",
                             help=f'The case-insensitive md-format of the input trajectory. Default is "pimd-qmcf" and possible oprions are: {MDEngineFormat.value_repr()}.')
