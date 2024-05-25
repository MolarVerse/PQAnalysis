"""
.. _cli.xyz2gen:

Command Line Tool for Converting XYZ Files to GEN Files
-----------------------------------------------------------


"""

from PQAnalysis.config import code_base_url
from PQAnalysis.io import xyz2gen
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to convert xyz files to gen files.
"""

__epilog__ = '\n'
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.xyz2gen.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class XYZ2GENCLI(CLIBase):

    """
    Command Line Tool for Converting XYZ Files to GEN Files
    """

    @classmethod
    def program_name(cls) -> str:
        """
        Returns the name of the program.

        Returns
        -------
        str
            The name of the program.
        """
        return 'xyz2gen'

    @classmethod
    def add_arguments(cls, parser: _ArgumentParser) -> None:
        """
        Adds the arguments to the parser.

        Parameters
        ----------
        parser : _ArgumentParser
            The parser to which the arguments should be added.
        """
        parser.parse_output_file()

        parser.add_argument(
            'xyz_file',
            type=str,
            help='The gen file to be converted.',
        )

        parser.add_argument(
            '--periodic',
            choices=[True, False, None],
            default=None,
            help=(
                'If True, the box is printed. If False, the box is not printed. '
                'If None, the box is printed if it is present in the xyz file.'
            )
        )

        parser.parse_engine()
        parser.parse_mode()

    @classmethod
    def run(cls, args):
        """
        Runs the command line tool.

        Parameters
        ----------
        args : _Namespace
            The arguments from the command line.
        """
        xyz2gen(
            xyz_file=args.xyz_file,
            output=args.output,
            periodic=args.periodic,
            md_format=args.engine,
            mode=args.mode,
        )



def main():
    """
    Main function of the xyz2gen command line tool, which is basically just
    a wrapper for the xyz2gen function. For more information on the xyz2gen
    function please visit :py:func:`PQAnalysis.io.api.xyz2gen`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    XYZ2GENCLI.add_arguments(parser)

    args = parser.parse_args()

    XYZ2GENCLI.run(args)
