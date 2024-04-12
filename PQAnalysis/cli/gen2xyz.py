"""
.. _cli.gen2xyz:

Command Line Tool for Converting GEN Files to XYZ Files
-----------------------------------------------------------


"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io import gen2xyz


__outputdoc__ = """

This command line tool can be used to convert gen files to xyz files. 

If the box information from the gen file should not be included in the xyz file, 
please use the --nobox option.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.gen2xyz.html.
"""

__doc__ += __outputdoc__


def main():
    """
    Main function of the gen2xyz command line tool, which is basically just a wrapper for the gen2xyz function. For more information on the gen2xyz function please visit :py:func:`PQAnalysis.io.api.gen2xyz`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    parser.add_argument('gen_file', type=str,
                        help='The gen file to be converted.')

    parser.parse_output_file()

    parser.add_argument('--nobox', action='store_true',
                        help='Do not print the box.')

    parser.parse_engine()
    parser.parse_mode()

    args = parser.parse_args()

    gen2xyz(
        gen_file=args.gen_file,
        output=args.output,
        print_box=not args.nobox,
        engine=args.engine,
        mode=args.mode,
    )
