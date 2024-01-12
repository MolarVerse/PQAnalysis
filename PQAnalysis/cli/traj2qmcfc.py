"""
.. _cli.traj2qmcfc:

Command Line Tool for Converting PIMD-QMCF to QMCFC Trajectory Files
---------------------------------------------------------------------------

Converts a PIMD-QMCF trajectory to a QMCFC trajectory format output.

Both formats are adapted xyz file formats with the box dimensions and box angles
being placed in the same line after the number of atoms. The QMCFC contains an 
additional dummy 'X' particle as first entry of the coordinates section for
visibility and debugging reasons in conjunction with vmd.
"""

import argparse

from ..io import traj2qmcfc


def main():
    """
    Wrapper for the command line interface of traj2qmcfc.
    """
    parser = argparse.ArgumentParser(
        description='Converts multiple trajectory files to a qmcfc trajectory.')
    parser.add_argument('trajectory_file', type=str, nargs='+',
                        help='The trajectory file(s) to be converted.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='The output file. If not specified, the output is printed to stdout.')
    args = parser.parse_args()

    traj2qmcfc(args.trajectory_file, args.output)
