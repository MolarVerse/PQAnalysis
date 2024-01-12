"""
.. _cli.rst2xyz:

Command Line Tool for Converting Restart Files to XYZ Files
-----------------------------------------------------------

If the box information from the restart file should not be included in the xyz file, 
please use the --nobox option.
"""

import argparse

from ..io import rst2xyz


def main():
    """
    Wrapper for the command line interface of rst2xyz.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('restart_file', type=str,
                        help='The restart file to be converted.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='The output file. If not specified, the output is printed to stdout.')
    parser.add_argument('--nobox', action='store_true',
                        help='Do not print the box.')
    args = parser.parse_args()

    rst2xyz(args.restart_file, args.output, not args.nobox)
