"""
Converts a restart file to a xyz file.

If the box information from the restart file should not be included in the xyz file, 
please use the --nobox option.
"""

import argparse

from ..io.restartReader import RestartFileReader
from ..io.trajectoryWriter import TrajectoryWriter
from ..traj.frame import Frame
from ..core.cell import Cell


def main():
    """
    Wrapper for the command line interface of rst2xyz.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('restart_file', type=str,
                        help='The restart file to be converted.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='The output file. If not specified, the output is printed to stdout.')
    parser.add_argument('--nobox', action='store_true',
                        help='Do not print the box.')
    args = parser.parse_args()

    rst2xyz(args.restart_file, args.output, not args.nobox)


def rst2xyz(restart_file: str, output: str | None = None, print_box: bool = True):
    """
    Converts a restart file to a xyz file and prints it to stdout or writes it to a file.

    Parameters
    ----------
    restart_file : str
        The restart file to be converted.
    output : str | None
        The output file. If not specified, the output is printed to stdout.
    print_box : bool
        If True, the box is printed. If False, the box is not printed. Default is True.
    """
    reader = RestartFileReader(restart_file)
    system, _ = reader.read()

    if not print_box:
        system.cell = Cell()

    frame = Frame(system)

    writer = TrajectoryWriter(filename=output)
    writer.write(frame, type="xyz")
