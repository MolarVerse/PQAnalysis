"""
.. _cli.traj2extxyz:

Command Line Tool for Converting PQ Trajectory Files to Extended XYZ
--------------------------------------------------------------------

"""

from PQAnalysis.config import code_base_url
from PQAnalysis.io import traj2extxyz
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """
Converts PQ xyz trajectory files to extended xyz format.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.traj2extxyz.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class Traj2ExtXYZCLI(CLIBase):

    """
    Command Line Tool for Converting PQ Trajectory Files to Extended XYZ
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
        return 'traj2extxyz'

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
            'trajectory_file',
            type=str,
            nargs='+',
            help='The trajectory file(s) to be converted.'
        )

        parser.parse_engine()
        parser.parse_mode()

    @classmethod
    def run(cls, args) -> None:
        """
        Runs the command line tool.

        Parameters
        ----------
        args : _Namespace
            The arguments from the command line.
        """
        traj2extxyz(
            args.trajectory_file,
            args.output,
            args.mode,
            args.engine,
        )



def main():
    """
    Main function of the traj2extxyz command line tool, which is basically just
    a wrapper for the traj2extxyz function. For more information on the
    traj2extxyz function please visit :py:func:`PQAnalysis.io.api.traj2extxyz`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    Traj2ExtXYZCLI.add_arguments(parser)

    args = parser.parse_args()

    Traj2ExtXYZCLI.run(args)
