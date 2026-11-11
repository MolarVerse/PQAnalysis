"""
.. _cli.traj2com:

Command Line Tool for Computing Center of Mass Trajectories
-----------------------------------------------------------

"""

from PQAnalysis.config import code_base_url
from PQAnalysis.tools import traj2com
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """
Computes the center of mass trajectory of trajectory files. The center of
mass is computed for all selected atoms, for groups of consecutive atoms
(--group) or for every molecule (--by-residue). Molecules are defined by a
restart file together with a moldescriptor file. The format is inferred from
the file extension: an xyz file gives the positions of the centers of mass,
a vel file their velocities, a force file their forces and a charge file
their charges.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.traj2com.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class Traj2COMCLI(CLIBase):

    """
    Command Line Tool for Computing Center of Mass Trajectories
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
        return 'traj2com'

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

        parser.parse_trajectory_file()

        parser.add_argument(
            '--selection',
            type=str,
            default=None,
            help=(
                'The atoms that define the centers of mass, e.g. "O" or '
                '"*|X" (all atoms except X). If not specified, all atoms '
                'are used.'
            )
        )

        parser.add_argument(
            '--use-full-atom-info',
            action='store_true',
            default=False,
            help='Use the full atom information for the selection.'
        )

        grouping = parser.add_mutually_exclusive_group()

        grouping.add_argument(
            '--group',
            type=int,
            default=None,
            metavar='N',
            help=(
                'Compute one center of mass for every N consecutive atoms '
                'of the selection. Without --group and --by-residue one '
                'center of mass of all selected atoms is computed.'
            )
        )

        grouping.add_argument(
            '--by-residue',
            action='store_true',
            default=False,
            help=(
                'Compute one center of mass for every molecule. Needs '
                '--restart-file and --moldescriptor-file.'
            )
        )

        parser.add_argument(
            '--restart-file',
            type=str,
            default=None,
            help='The restart file that defines the molecules.'
        )

        parser.add_argument(
            '--moldescriptor-file',
            type=str,
            default=None,
            help='The moldescriptor file with the molecule types.'
        )

        parser.add_argument(
            '--continuous',
            action='store_true',
            default=False,
            help=(
                'Follow the centers over the periodic boundary instead of '
                'folding them back into the cell, e.g. for the mean square '
                'displacement. Needs positions and a center that moves less '
                'than half of the box between two frames.'
            )
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
        traj2com(
            trajectory_files=args.trajectory_file,
            output=args.output,
            mode=args.mode,
            md_format=args.engine,
            selection=args.selection,
            use_full_atom_info=args.use_full_atom_info,
            group=args.group,
            by_residue=args.by_residue,
            restart_file=args.restart_file,
            moldescriptor_file=args.moldescriptor_file,
            continuous=args.continuous,
        )



def main():
    """
    Main function of the traj2com command line tool, which is basically just
    a wrapper for the traj2com function. For more information on the
    traj2com function please visit
    :py:func:`PQAnalysis.tools.traj_to_com_traj.traj2com`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    Traj2COMCLI.add_arguments(parser)

    args = parser.parse_args()

    Traj2COMCLI.run(args)
