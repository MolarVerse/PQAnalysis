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
            f"Format {format} not implemented yet for continuing input file.")

    continue_input_file(args.input_file, args.number, format)


def continue_input_file(filename: str, n: PositiveInt, format: InputFileFormat):
    """
    Creates n new input files by increasing the number in the filename by one.

    The input file must contain a number before the file extension from which the
    new filenames are created. The number is increased by one for each new file.
    All other numbers in the start- and output-files within the input file are increased by one as well.

    Parameters
    ----------
    filename : str
        the input file
    n : PositiveInt
        number of new input files to be created
    format : InputFileFormat
        the format of the input file
    """
    reader = Reader(filename)
    reader.read()
    reader.continue_input_file(n)
