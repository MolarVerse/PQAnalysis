"""
.. _cli.parm2pq:

Command Line Tool for Importing AMBER, CHARMM and GROMACS Force Fields
----------------------------------------------------------------------


"""

from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool converts a force field read by ParmEd (AMBER
prmtop + inpcrd, GROMACS top + gro, CHARMM psf + parameter files +
coordinates) into PQ setup files: moldescriptor.dat, topology.top,
parameter.prm, intra_nonbonded.dat, start.rst and the input snippet
pq_input.in with the matching keys.

Each connected molecule becomes a molecule of PQ (protein chains are
split); identical molecules share one moltype. Bonds and angles get the
PQ force constant (2 k), torsion phases the PQ sign convention, and every
pair of Lennard-Jones types one NONCOULOMBICS LJ line. Exclusions and
scaled 1-4 pairs go into the intra non-bonded file; bonded terms between
molecules carry the linker mark.

CMAP, Urey-Bradley terms, harmonic impropers, non-uniform 1-4 scaling,
special 1-4 Lennard-Jones parameters and virtual sites are not supported
and stop the conversion with an error naming the term.

The tool needs ParmEd: pip install PQAnalysis[ff].
"""

__epilog__ = "\n"
__epilog__ += "For more information on the options please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.parm2pq.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class Parm2PQCLI(CLIBase):

    """
    Command Line Tool for Importing Force Fields into PQ Setup Files
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
        return 'parm2pq'

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
            'topology_file',
            type=str,
            help='The topology: AMBER prmtop, GROMACS top or CHARMM psf.',
        )

        parser.add_argument(
            '-c',
            '--coordinates',
            type=str,
            default=None,
            help=(
                'The coordinate file (inpcrd/rst7, gro, crd, pdb). Its box '
                'is used; without box a cubic vacuum box of the extent plus '
                'two Coulomb cutoffs is written.'
            ),
        )

        parser.add_argument(
            '-p',
            '--parameters',
            type=str,
            nargs='+',
            default=[],
            help='CHARMM parameter and stream files for a psf topology.',
        )

        parser.add_argument(
            '-d',
            '--output-dir',
            type=str,
            default='.',
            help='The directory for the PQ setup files (default: .).',
        )

        parser.add_argument(
            '--shake-h',
            action='store_true',
            help=(
                'Constrain bonds to hydrogen with SHAKE instead of writing '
                'them as bond terms.'
            ),
        )

        parser.add_argument(
            '--rcoulomb',
            type=float,
            default=12.0,
            help=(
                'The Coulomb cutoff in Angstrom for the input snippet '
                '(default: 12.0).'
            ),
        )

        parser.add_argument(
            '--long-range',
            type=str,
            choices=['none', 'wolf', 'reaction-field'],
            default=None,
            help=(
                'The long-range Coulomb method for the input snippet '
                '(default: wolf with a box from the coordinates, none for '
                'the vacuum box).'
            ),
        )

        parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Replace existing output files.',
        )

    @classmethod
    def run(cls, args):
        """
        Runs the command line tool.

        Parameters
        ----------
        args : _Namespace
            The arguments from the command line.
        """
        from PQAnalysis.io.force_field import parm2pq  # pylint: disable=import-outside-toplevel

        parm2pq(
            topology_file=args.topology_file,
            coordinate_file=args.coordinates,
            parameter_files=args.parameters,
            output_dir=args.output_dir,
            shake_h=args.shake_h,
            rcoulomb=args.rcoulomb,
            overwrite=args.overwrite,
            long_range=args.long_range,
        )



def main():
    """
    Main function of the parm2pq command line tool, which is a wrapper for
    :py:func:`PQAnalysis.io.force_field.api.parm2pq`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    Parm2PQCLI.add_arguments(parser)

    args = parser.parse_args()

    Parm2PQCLI.run(args)
