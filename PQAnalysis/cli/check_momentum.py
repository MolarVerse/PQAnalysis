"""
.. _cli.check_momentum:

Command Line Tool for Checking the Total Linear Momentum
========================================================


"""

from PQAnalysis.analysis.momentum import check_momentum
from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to calculate the norm of the total
linear momentum P = sum_i m_i * v_i of (a selection of) atoms for
every frame of a velocity trajectory. It can be used to check a
simulation for center of mass drift.

For every frame one row containing the one-based frame index and the
scaled momentum norm is written. With velocities in Angstrom/s (PQ
velocity trajectories) the default scaling factor of 1e-15 converts
the momentum norm from amu*Angstrom/s to amu*Angstrom/fs.

Note that the velocities are parsed from file in single precision, so
reported norms below roughly 1e-7 * sum_i m_i * |v_i| * scale are
parsing noise, not physical center of mass drift. The legacy
equipartition.jl tool parses the velocities in double precision and
therefore resolves correspondingly smaller drift for
momentum-conserving trajectories.
"""

__epilog__ = "\n"
__epilog__ += "For more information on the command line options of this tool please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.check_momentum.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class CheckMomentumCLI(CLIBase):

    """
    Command Line Tool for Checking the Total Linear Momentum
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
        return 'check_momentum'

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
                'The selection of atoms to include in the total '
                'momentum. If not specified, all atoms are included.'
            )
        )

        parser.add_argument(
            '--use-full-atom-info',
            action='store_true',
            default=False,
            help='Use the full atom information for the selection.'
        )

        parser.add_argument(
            '--scale',
            type=float,
            default=1e-15,
            help=(
                'The scaling factor applied to the momentum norm '
                'before output. The default converts amu*Angstrom/s '
                'to amu*Angstrom/fs.'
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
        args : argparse.Namespace
            The arguments parsed by the parser.
        """
        check_momentum(
            trajectory_files=args.trajectory_file,
            output=args.output,
            selection=args.selection,
            use_full_atom_info=args.use_full_atom_info,
            scale=args.scale,
            md_format=args.engine,
            mode=args.mode,
        )



def main():
    """
    Main function of the check_momentum command line tool, which is
    basically just a wrapper for the check_momentum function. For more
    information on the check_momentum function please visit
    :py:func:`PQAnalysis.analysis.momentum.api.check_momentum`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    CheckMomentumCLI.add_arguments(parser)

    args = parser.parse_args()

    CheckMomentumCLI.run(args)
