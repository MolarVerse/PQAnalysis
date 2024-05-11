"""
.. _cli.traj2qmcfc:

Command Line Tool for Converting PQ to QMCFC Trajectory Files
---------------------------------------------------------------------------

"""

from PQAnalysis.config import code_base_url
from PQAnalysis.io import traj2qmcfc
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """
Converts a PQ trajectory to a QMCFC trajectory format output.

Both formats are adapted xyz file formats with the box dimensions and box angles
being placed in the same line after the number of atoms. The QMCFC contains an 
additional dummy 'X' particle as first entry of the coordinates section for
visibility and debugging reasons in conjunction with vmd.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.traj2qmcfc.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class Traj2QMCFCCLI(CLIBase):

    """
    Command Line Tool for Converting PQ to QMCFC Trajectory Files
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
        return 'traj2qmcfc'

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
            help='The trajectory file to be converted.'
        )

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
        traj2qmcfc(args.trajectory_file, args.output)



def main():
    """
    Main function of the traj2qmcfc command line tool, which is basically just a 
    wrapper for the traj2qmcfc function. For more information on the traj2qmcfc 
    function please visit :py:func:`PQAnalysis.io.api.traj2qmcfc`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    Traj2QMCFCCLI.add_arguments(parser)

    args = parser.parse_args()

    Traj2QMCFCCLI.run(args)
