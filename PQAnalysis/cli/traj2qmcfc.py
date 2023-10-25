import numpy as np
import argparse

from PQAnalysis.traj.trajectory import read


def main():
    parser = argparse.ArgumentParser(
        description='Converts multiple trajectory files to a qmcfc trajectory.')
    parser.add_argument('trajectory_file', type=str, nargs='+',
                        help='The trajectory file(s) to be converted.')
    args = parser.parse_args()

    for filename in args.trajectory_file:
        trajectory = read(filename)
        for frame in trajectory:
            print_qmcfc_trajectory(frame)


def print_qmcfc_trajectory(frame):
    frame.print_xyz_header()
    print()
    print("X 0.0 0.0 0.0")
    frame.print_coordinates()
