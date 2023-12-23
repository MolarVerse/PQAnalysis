"""
Creates n new input files.

This command line tool generates n new input files by increasing the number in the filename by one.
Furthermore, all other numbers in the start- and output-files within the input file are increased by one as well.

With this feature it is possible to continue a simulation with a new input file without having to change the input file manually. 
"""

import argparse

from ..types import PositiveInt
from ..io import InputFileFormat
from ..io import PIMD_QMCF_InputFileReader as Reader


def main():
    """
    Wrapper for the command line interface of continue_input.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file', type=str, help='The input file.')
    parser.add_argument('-n', '--number', type=int, default=1,
                        help='The number of times the input file should be continued.')
    parser.add_argument('-f', '--format', type=str, default='pimd-qmcf',
                        help='The format of the input file. Default is pimd-qmcf.')

    args = parser.parse_args()

    format = InputFileFormat(args.format)

    if format != InputFileFormat.PIMD_QMCF:
        raise NotImplementedError(
            f"Format {args.format} not implemented yet for continuing input file.")

    reader = Reader(args.input_file)
    reader.read()
    reader.continue_input_file(args.number)
