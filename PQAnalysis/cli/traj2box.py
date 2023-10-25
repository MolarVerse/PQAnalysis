import numpy as np
import argparse

from PQAnalysis.traj.trajectory import read
from PQAnalysis.utils import decorators


def main():
    parser = argparse.ArgumentParser(
        description='Converts multiple trajectory files to a box file.')
    parser.add_argument('trajectory_file', type=str, nargs='+',
                        help='The trajectory file(s) to be converted.')
    parser.add_argument('-v', '--vmd', action='store_true',
                        help='Output in VMD format.')
    args = parser.parse_args()

    for filename in args.trajectory_file:
        trajectory = read(filename)
        for frame in trajectory:
            if args.vmd:
                print_vmd_box(frame.cell)
            else:
                print_box(frame.cell)


def print_vmd_box(cell):
    print("8")
    print(f"Box {cell.x} {cell.y} {cell.z}")

    for x in [-0.5, 0.5]:
        for y in [-0.5, 0.5]:
            for z in [-0.5, 0.5]:
                vec = cell.box_matrix @ np.array([x, y, z])
                print(f"X {vec[0]} {vec[1]} {vec[2]}")


@decorators.count_decorator
def print_box(cell):
    print(f"{print_box.counter} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}")
