import argparse


def main():
    """
    Wrapper for the command line interface of xyz2rst.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('xyz_file', type=str,
                        help='The xyz file to be converted.')
    parser.add_argument('vel_file', type=str, default=None,
                        help='The velocity file to be converted.')
    parser.add_argument('force_file', type=str, default=None,
                        help='The force file to be converted.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='The output file. If not specified, the output is printed to stdout.')

    args = parser.parse_args()

    # xyz2rst(args.xyz_file, args.vel_file, args.force_file, args.output)

# def xyz2rst(xyz_file: str, vel_file: str, force_file: str, output: str | None = None):
#     """
#     Converts a xyz file to a restart file and prints it to stdout or writes it to a file.

#     Parameters
#     ----------
#     xyz_file : str
#         The xyz file to be converted.
#     output : str | None
#         The output file. If not specified, the output is printed to stdout.
#     """
#     writer = RestartFileWriter(filename=output)
#     writer.write(frame, type="rst"
