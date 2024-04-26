"""
.. _cli.rst2xyz:

Command Line Tool for Converting Restart Files to XYZ Files
-----------------------------------------------------------


"""


from PQAnalysis.io import rst2xyz
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser


__outputdoc__ = """

This command line tool can be used to convert restart files to xyz files. 

If the box information from the restart file should not be included in the xyz file, 
please use the --nobox option.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.rst2xyz.html."
__epilog__ += "\n"

__doc__ += __outputdoc__


def main():
    """
    Main function of the rst2xyz command line tool, which is basically just a wrapper 
    for the rst2xyz function. For more information on the rst2xyz function 
    please visit :py:func:`PQAnalysis.io.api.rst2xyz`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    parser.parse_output_file()

    parser.add_argument(
        'restart_file',
        type=str,
        help='The restart file to be converted.'
    )

    parser.add_argument(
        '--nobox',
        action='store_true',
        help='Do not print the box.'
    )

    parser.parse_engine()
    parser.parse_mode()

    args = parser.parse_args()

    rst2xyz(
        restart_file=args.restart_file,
        output=args.output,
        print_box=not args.nobox,
        md_format=args.engine,
        mode=args.mode,
    )
