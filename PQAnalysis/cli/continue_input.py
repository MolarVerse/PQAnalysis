"""
.. _cli.continue_input:

Command Line Tool for Extending MD Simulation Input Files
---------------------------------------------------------

"""
import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io import InputFileFormat, continue_input_file

__outputdoc__ = """

This command line tool generates n new input files by increasing the number in the filename by one.
Furthermore, all other numbers in the start- and output-files within the input file are increased by one as well.

With this feature it is possible to continue a simulation with a new input file without having to change the input file manually. 
"""

__doc__ += __outputdoc__

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.continue_input.html.
"""


def main():
    """
    The main function of the continue input command line tool which is basically just a wrapper for the continue_input_file function. For more information on the continue_input_file function please visit :py:func:`PQAnalysis.io.api.continue_input_file`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)
    parser.parse_input_file()
    parser.add_argument('-n', '--number', type=int, default=1,
                        help='The number of times the input file should be continued.')
    parser.add_argument('--input-format', type=str, default='PQ',
                        help='The format of the input file. Default is PQ.')

    args = parser.parse_args()

    input_format = InputFileFormat(args.input_format)

    continue_input_file(args.input_file, args.number, input_format)
