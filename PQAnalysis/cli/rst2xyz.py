"""
.. _cli.rst2xyz:

Command Line Tool for Converting Restart Files to XYZ Files
-----------------------------------------------------------


"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io import rst2xyz


__outputdoc__ = """

This command line tool can be used to convert restart files to xyz files. 

If the box information from the restart file should not be included in the xyz file, 
please use the --nobox option.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.rst2xyz.html.
"""

__doc__ += __outputdoc__


def main():
    """
    Main function of the rst2xyz command line tool, which is basically just a wrapper for the rst2xyz function. For more information on the rst2xyz function please visit :py:func:`PQAnalysis.io.api.rst2xyz`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    parser.parse_output_file()

    parser.add_argument('restart_file', type=str,
                        help='The restart file to be converted.')
    parser.add_argument('--nobox', action='store_true',
                        help='Do not print the box.')
    args = parser.parse_args()

    rst2xyz(args.restart_file, args.output, not args.nobox)
