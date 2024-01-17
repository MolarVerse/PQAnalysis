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

from ._argumentParser import ArgumentParser
from PQAnalysis.io import traj2qmcfc


def main():
    """
    Wrapper for the command line interface of traj2qmcfc.
    """
    parser = ArgumentParser(
        description='Converts multiple trajectory files to a qmcfc trajectory.')

    parser.parse_trajectory_file()
    parser.parse_output_file()

    args = parser.parse_args()

    traj2qmcfc(args.trajectory_file, args.output)
