"""
Converts a PIMD-QMCF trajectory to a QMCFC trajectory format output.

Both formats are adapted xyz file formats with the box dimensions and box angles
being placed in the same line after the number of atoms. The QMCFC contains an 
additional dummy 'X' particle as first entry of the coordinates section for
visibility and debugging reasons in conjunction with vmd.
"""

import argparse

from beartype.typing import List

from PQAnalysis.io.trajectoryWriter import TrajectoryWriter
from PQAnalysis.io.trajectoryReader import TrajectoryReader


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


def traj2qmcfc(trajectory_files: List[str], output: str):
    """
    Converts multiple trajectory files to a qmcfc trajectory.

    Parameters
    ----------
    trajectory_file : list of str
        The trajectory file(s) to be converted.
    output : str, optional
        The output file. If not specified, the output is printed to stdout.
    """
    writer = TrajectoryWriter(filename=output, format="qmcfc")
    for filename in trajectory_files:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory)

    import os
    os.system("cat test_traj.qmcfc.xyz")
