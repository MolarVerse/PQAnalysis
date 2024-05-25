"""
.. _cli.gen2xyz:

Command Line Tool for Converting GEN Files to XYZ Files
-----------------------------------------------------------


"""

from PQAnalysis.io import gen2xyz
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to convert gen files to xyz files. 

If the box information from the gen file should not be included in the xyz file, 
please use the --nobox option.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.gen2xyz.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class GEN2XYZCLI(CLIBase):

    """
    Command Line Tool for Converting GEN Files to XYZ Files
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
        return 'gen2xyz'

    @classmethod
    def add_arguments(cls, parser: _ArgumentParser) -> None:
        """
        Adds the arguments to the parser.

        Parameters
        ----------
        parser : _ArgumentParser
            The parser to which the arguments should be added.
        """
        parser.add_argument(
            'gen_file', type=str, help='The gen file to be converted.'
        )

        parser.parse_output_file()

        parser.add_argument(
            '--nobox',
            action='store_true',
            default=False,
            help='Do not print the box.'
        )

        parser.parse_engine()
        parser.parse_mode()

    @classmethod
    def run(cls, args):
        """
        Runs the command line tool.

        Parameters
        ----------
        args : _ArgumentParser
            The arguments that were parsed by the parser.
        """
        gen2xyz(
            gen_file=args.gen_file,
            output=args.output,
            print_box=not args.nobox,
            md_format=args.engine,
            mode=args.mode,
        )



def main():
    """
    Main function of the gen2xyz command line tool, which is basically just 
    a wrapper for the gen2xyz function. For more information on the gen2xyz
    function please visit :py:func:`PQAnalysis.io.api.gen2xyz`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    GEN2XYZCLI.add_arguments(parser)

    args = parser.parse_args()

    GEN2XYZCLI.run(args)
