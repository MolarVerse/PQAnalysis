"""
.. _cli.continue_input:

Command Line Tool for Extending MD Simulation Input Files
---------------------------------------------------------

"""
from ._argumentParser import ArgumentParser
from PQAnalysis.io import InputFileFormat, continue_input_file

__outputdoc__ = """

This command line tool generates n new input files by increasing the number in the filename by one.
Furthermore, all other numbers in the start- and output-files within the input file are increased by one as well.

With this feature it is possible to continue a simulation with a new input file without having to change the input file manually. 
"""

__doc__ += __outputdoc__


def main():
    """
    The main function of the continue input command line tool.
    """
    parser = ArgumentParser(description=__outputdoc__)
    parser.parse_input_file()
    parser.add_argument('-n', '--number', type=int, default=1,
                        help='The number of times the input file should be continued.')
    parser.add_argument('--input-format', type=str, default='pimd-qmcf',
                        help='The format of the input file. Default is pimd-qmcf.')

    args = parser.parse_args()

    input_format = InputFileFormat(args.input_format)

    continue_input_file(args.input_file, args.number, input_format)
