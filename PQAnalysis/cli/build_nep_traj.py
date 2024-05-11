"""
.. _cli.build_nep_traj:

Command Line Tool for Building Neuroevolution Potential (NEP) training/test trajectories
----------------------------------------------------------------------------------------


"""

from PQAnalysis.io.nep.nep_writer import NEPWriter
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to convert output of PQ of QMCFC simulations to training and test files for the Neuroevolution Potential (NEP) method. The output is written to a xyz file.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.build_nep_traj.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class BuildNEPTrajCLI(CLIBase):

    """
    Command Line Tool for Building Neuroevolution Potential (NEP) training/test trajectories
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
        return 'build_nep_traj'

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
            'file_prefixes',
            type=str,
            nargs='+',
            help='The file prefixes for all input files.'
        )

        parser.parse_output_file()

        parser.add_argument(
            '--test-ratio',
            type=float,
            default=0.0,
            help=(
            "The ratio of testing frames to the total number of "
            "frames, by default 0.0. If the test_ratio is 0.0 no "
            "train and test files are created. If the test_ratio "
            "is larger not equal to 0.0, the test_ratio is used "
            "to determine the number of training and testing frames. "
            "The final ratio will be as close to the test_ratio as "
            "possible, but if it is not possible to have the exact "
            "ratio, always the higher next higher ratio is chosen. "
            "As output filenames the original filename is used with "
            "the suffix _train or _test appended and the same "
            "FileWritingMode as the original file is used."
            )
        )

        parser.add_argument(
            '--total-ratios',
            type=str,
            default=None,
            help=(
            "The total_ratios keyword argument is used to "
            "describe frame ratios including validation frames "
            "in the format train_ratio:test_ratio:validation_ratio. "
            "The validation_ratio is optional and if not given, "
            "no validation frames are written. The total sum of "
            "the integer values provided do not have to add up "
            "to the total number of frames in the input trajectory "
            "files. The ratios are used to determine the ratios of "
            "the training, testing, and validation frames. The final"
            "ratio will be as close to the given ratios as possible, "
            "but if it is not possible to have the exact ratio, "
            "always the next higher ratio is chosen. As output "
            "filenames the original filename is used with the suffix "
            "_train, _test, or _validation appended and the same "
            "FileWritingMode as the original file is used. The "
            "validation frames are written to a file with the "
            "suffix _validation and a file with the suffix _validation.ref. "
            "The _validation file contains only the coordinates and "
            "box information to function as crude testing input and "
            "the _validation.ref file contains all information "
            "additionally provided in the original files. "
            "Pay Attention: This keyword argument is mutually exclusive "
            "with the test_ratio keyword argument. If both are given, "
            "a ValueError is raised."
            )
        )

        parser.add_argument(
            '--use-forces',
            action='store_true',
            default=False,
            help='Whether to include forces in the output file.'
        )

        parser.add_argument(
            '--use-virial',
            action='store_true',
            default=False,
            help='Whether to include the virial in the output file.'
        )

        parser.add_argument(
            '--use-stress',
            action='store_true',
            default=False,
            help='Whether to include the stress tensor in the output file.'
        )

        parser.parse_mode()

    @classmethod
    def run(cls, args):
        """
        Runs the build_nep_traj function with the given arguments.

        Parameters
        ----------
        args : argparse.Namespace
            The parsed arguments.
        """
        writer = NEPWriter(filename=args.output, mode=args.mode)

        writer.write_from_files(
            file_prefixes=args.file_prefixes,
            test_ratio=args.test_ratio,
            total_ratios=args.total_ratios,
            use_forces=args.use_forces,
            use_virial=args.use_virial,
            use_stress=args.use_stress
        )



def main():
    """
    Main function of the build_nep_traj command line tool, which is basically 
    just a wrapper for the build_nep_traj function. For more information on 
    the build_nep_traj function please visit :py:func:`PQAnalysis.io.api.build_nep_traj`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    BuildNEPTrajCLI.add_arguments(parser)

    args = parser.parse_args()

    BuildNEPTrajCLI.run(args)
