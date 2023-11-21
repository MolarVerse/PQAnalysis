import argparse


def main():
    """
    Wrapper for the command line interface of xyz2rst.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('xyz_file', type=str,
                        help='The xyz file to be converted.')
    parser.add_argument('vel_file', type=str,
                        help='The velocity file to be converted.')
    parser.add_argument('force_file', type=str,
                        help='The force file to be converted.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='The output file. If not specified, the output is printed to stdout.')
