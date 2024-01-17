"""
.. _cli.rst2xyz:

Command Line Tool for Converting Restart Files to XYZ Files
-----------------------------------------------------------

If the box information from the restart file should not be included in the xyz file, 
please use the --nobox option.
"""

from ._argumentParser import ArgumentParser
from PQAnalysis.io import rst2xyz


def main():
    """
    Wrapper for the command line interface of rst2xyz.
    """
    parser = ArgumentParser(description=__doc__)

    parser.parse_output_file()

    parser.add_argument('restart_file', type=str,
                        help='The restart file to be converted.')
    parser.add_argument('--nobox', action='store_true',
                        help='Do not print the box.')
    args = parser.parse_args()

    rst2xyz(args.restart_file, args.output, not args.nobox)
