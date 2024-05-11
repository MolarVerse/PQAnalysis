"""
.. _cli.add_molecules:

Command Line Tool for Converting Restart Files to XYZ Files
-----------------------------------------------------------


"""

from PQAnalysis.config import code_base_url
from PQAnalysis.tools.add_molecule import add_molecule
from PQAnalysis.io.formats import OutputFileFormat

from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to add molecules to a restart file.

The molecules are added by fitting the molecule to the restart file. The fitting is done randomly by rotating the molecule and translating it to a random position. After the fitting, the molecule is added to the restart file. The class can add multiple molecules to the restart file.

This tool provides also the possibility of not only extending the restart file with the molecule but also the topology file. Therefore, the topology file of the restart file and the topology file of the molecule file have to be provided as well as the desired output topology file.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input "
__epilog__ += "file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.add_molecules.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class AddMoleculesCLI(CLIBase):

    """
    Command Line Tool for Adding Molecules to Restart Files
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
        return 'add_molecules'

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
            'restart_file',
            type=str,
            help='The restart file where the molecules should be added.'
        )

        parser.add_argument(
            'molecule_file',
            type=str,
            help=(
            "The molecule file that contains the coordinates "
            "of the molecule that should be added. Can be in "
            "any format that is supported by the PQAnalysis library."
            )
        )

        parser.parse_output_file()
        parser.parse_mode()

        parser.add_argument(
            '--mol-file-type',
            dest='mol_file_type',
            type=OutputFileFormat,
            default=OutputFileFormat.AUTO,
            choices=OutputFileFormat.__members__.values(),
            help=(
            'The file format of the molecule file. '
            'If not specified, the file format will '
            'be inferred from the file extension.'
            )
        )

        parser.add_argument(
            "--rst-mol-desc-file",
            dest='rst_mol_desc_file',
            type=str,
            help=(
            "The moldescriptor file that is associated with the "
            "restart file. If not specified, the moldescriptor "
            "file will not be used."
            ),
            default=None
        )

        parser.add_argument(
            "--molecule-mol-desc-file",
            dest='molecule_mol_desc_file',
            type=str,
            help=(
            "The moldescriptor file that is associated with "
            "the molecule file. If not specified, the moldescriptor "
            "file will not be used. Can only be used if the "
            "molecule file is a restart file type."
            ),
            default=None
        )

        parser.add_argument(
            "-n, --n-molecules",
            dest='n_molecules',
            type=int,
            default=1,
            help=
            "The number of molecules that should be added to the restart file."
        )

        parser.add_argument(
            "--max-iter",
            dest='max_iter',
            type=int,
            default=100,
            help=(
            "The maximum number of iterations that should "
            "be used to fit the molecule to the restart file."
            )
        )

        parser.add_argument(
            "-c",
            "--cut",
            dest='cut',
            type=float,
            default=1.0,
            help=(
            "The distance cutoff that should be used "
            "to fit the molecule to the restart file in Angstrom."
            )
        )

        parser.add_argument(
            "--max-disp",
            "--max-displacement",
            dest='max_disp',
            type=float,
            default=0.1,
            help=(
            "The maximum displacement that should be applied "
            "to the given molecule geometry relative to "
            "its center of mass in percentage."
            )
        )

        parser.add_argument(
            "--rot",
            "--rotation-angle-step",
            dest='rot',
            type=int,
            default=10,
            help=(
            "If the randomly placed molecule does not "
            "fit into the restart file, the molecule "
            "is rotated by the given angle step in degrees."
            )
        )

        parser.add_argument(
            "--topology-file",
            "--top-file",
            dest='top_file',
            type=str,
            help=(
            "The topology file that is associated with "
            "the restart file. If not specified, "
            "the topology file will not be used."
            ),
            default=None
        )

        parser.add_argument(
            "--added-topology-file",
            "--added-top-file",
            dest='added_top_file',
            type=str,
            help=(
            "The topology file that is associated with "
            "the molecule file. If not specified, the "
            "topology file will not be used."
            ),
            default=None
        )

        parser.add_argument(
            "--output-topology-file",
            "--output-top-file",
            dest='output_top_file',
            type=str,
            help=
            "The output topology file. If not specified, the output is printed to stdout.",
            default=None
        )

        parser.parse_engine()

    @classmethod
    def run(cls, args):
        """
        Runs the add_molecules function.

        Parameters
        ----------
        args : argparse.Namespace
            The parsed arguments.
        """
        add_molecule(
            args.restart_file,
            args.molecule_file,
            output_file=args.output,
            molecule_file_type=args.mol_file_type,
            restart_moldescriptor_file=args.rst_mol_desc_file,
            molecule_moldescriptor_file=args.molecule_mol_desc_file,
            number_of_additions=args.n_molecules,
            max_iterations=args.max_iter,
            distance_cutoff=args.cut,
            max_displacement=args.max_disp,
            rotation_angle_step=args.rot,
            md_engine_format=args.engine,
            mode=args.mode,
            topology_file=args.top_file,
            topology_file_to_add=args.added_top_file,
            topology_file_output=args.output_top_file
        )



def main():
    """
    Main function of the add_molecules command line tool,
    which is basically just a wrapper for the add_molecules
    function. For more information on the add_molecules
    function please visit :py:func:`PQAnalysis.tools.add_molecule`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    AddMoleculesCLI.add_arguments(parser)

    args = parser.parse_args()

    AddMoleculesCLI.run(args)
