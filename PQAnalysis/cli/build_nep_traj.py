"""
.. _cli.build_nep_traj:

Command Line Tool for Building Neuroevolution Potential (NEP) training/test trajectories
----------------------------------------------------------------------------------------


"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io.nep.nep import NEPWriter


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

    parser.add_argument('file_prefixes', type=str, nargs='+',
                        help='The file prefixes for all input files.')

    parser.parse_output_file()

    parser.add_argument('--use-forces', type=bool, default=False,
                        help='Whether to include forces in the output file.')

    parser.add_argument('--use-virial', type=bool, default=False,
                        help='Whether to include the virial in the output file.')

    parser.add_argument('--use-stress', type=bool, default=False,
                        help='Whether to include the stress tensor in the output file.')

    parser.parse_mode()

    args = parser.parse_args()

    writer = NEPWriter(filename=args.output, mode=args.mode)

    writer.write_from_files(
        file_prefixes=args.file_prefixes,
        use_forces=args.use_forces,
        use_virial=args.use_virial,
        use_stress=args.use_stress
    )
