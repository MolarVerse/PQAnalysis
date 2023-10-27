import argparse

from PQAnalysis.io.boxWriter import BoxWriter
from PQAnalysis.io.trajectoryReader import TrajectoryReader

"""


"""


def main():
    parser = argparse.ArgumentParser(
        description='Converts multiple trajectory files to a box file.')
    parser.add_argument('trajectory_file', type=str, nargs='+',
                        help='The trajectory file(s) to be converted.')
    parser.add_argument('-v', '--vmd', action='store_true',
                        help='Output in VMD format.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='The output file. If not specified, the output is printed to stdout.')
    args = parser.parse_args()

    if args.vmd:
        output_format = "vmd"
    else:
        output_format = None

    writer = BoxWriter(filename=args.output, format=output_format)
    for filename in args.trajectory_file:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory)
