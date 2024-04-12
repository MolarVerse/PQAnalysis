"""
.. _cli.xyz2gen:

Command Line Tool for Converting XYZ Files to GEN Files
-----------------------------------------------------------


"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io import xyz2gen


__outputdoc__ = """

This command line tool can be used to convert xyz files to gen files. 
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.xyz2gen.html.
"""

__doc__ += __outputdoc__


def main():
    """
    Main function of the xyz2gen command line tool, which is basically just a wrapper for the xyz2gen function. For more information on the xyz2gen function please visit :py:func:`PQAnalysis.io.api.xyz2gen`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    parser.add_argument('xyz_file', type=str,
                        help='The gen file to be converted.')

    parser.parse_output_file()

    parser.add_argument('periodic', choices=[True, False, None], default=None,
                        help='If True, the box is printed. If False, the box is not printed. If None, the box is printed if it is present in the xyz file.')

    parser.parse_engine()
    parser.parse_mode()

    args = parser.parse_args()

    xyz2gen(
        xyz_file=args.xyz_file,
        output=args.output,
        periodic=args.periodic,
        md_format=args.engine,
        mode=args.mode,
    )
