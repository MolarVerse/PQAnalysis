"""
.. _cli.rst2xyz:

Command Line Tool for Converting XYZ Files to Restart Files
-----------------------------------------------------------


"""

from PQAnalysis.io import xyz2rst
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to convert xyz files to restart files. 

"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.xyz2rst.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__


class XYZ2RstCLI(CLIBase):

    """
    Command Line Tool for Converting XYZ Files to Restart Files
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
        return 'xyz2rst'

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
            help='The xyz file to be converted.'
        )

        parser.add_argument(
            '--velocity_file',
            type=str,
            default=None,
            help='The velocity file to be converted.'
        )

        parser.add_argument(
            '--force_file',
            type=str,
            default=None,
            help='The force file to be converted.'
        )

        parser.add_argument(
            '--randomize',
            type=float,
            default=0.0,
            help='Randomize the atom order. Default is 0.0.'
        )

        parser.parse_engine()
        parser.parse_mode()

    @classmethod
    def run(cls, args):
        """
        Runs the command line tool.

        Parameters
        ----------
        args : argparse.Namespace
            The arguments parsed by the parser.
        """
        xyz2rst(
            xyz_file=args.xyz_file,
            velocity_file=args.velocity_file,
            force_file=args.force_file,
            randomize=args.randomize,
            output=args.output,
            md_format=args.engine,
            mode=args.mode,
        )


def main():
    """
    Main function of the xyz2rst command line tool, which is basically just a wrapper 
    for the xyz2rst function. For more information on the xyz2rst function 
    please visit :py:func:`PQAnalysis.io.api.xyz2rst`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    XYZ2RstCLI.add_arguments(parser)

    args = parser.parse_args()

    XYZ2RstCLI.run(args)
