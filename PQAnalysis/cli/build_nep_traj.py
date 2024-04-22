"""
.. _cli.build_nep_traj:

Command Line Tool for Building Neuroevolution Potential (NEP) training/test trajectories
----------------------------------------------------------------------------------------


"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io.nep.NEPWriter import NEPWriter


__outputdoc__ = """

This command line tool can be used to converts output of PQ of QMCFC simulations to training and test files for the Neuroevolution Potential (NEP) method. The output is written to a xyz file.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.build_nep_traj.html.
"""

__doc__ += __outputdoc__


def main():
    """
    Main function of the build_nep_traj command line tool, which is basically just a wrapper for the build_nep_traj function. For more information on the build_nep_traj function please visit :py:func:`PQAnalysis.io.api.build_nep_traj`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

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
        help="The ratio of testing frames to the total number of frames, by default 0.0. If the test_ratio is 0.0 no train and test files are created. If the test_ratio is larger not equal to 0.0, the test_ratio is used to determine the number of training and testing frames. The final ratio will be as close to the test_ratio as possible, but if it is not possible to have the exact ratio, always the higher next higher ratio is chosen. As output filenames the original filename is used with the suffix _train or _test appended and the same FileWritingMode as the original file is used."
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

    args = parser.parse_args()

    writer = NEPWriter(filename=args.output, mode=args.mode)

    writer.write_from_files(
        file_prefixes=args.file_prefixes,
        test_ratio=args.test_ratio,
        use_forces=args.use_forces,
        use_virial=args.use_virial,
        use_stress=args.use_stress
    )
