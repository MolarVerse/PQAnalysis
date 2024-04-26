"""
.. _cli.traj2box:

Command Line Tool for Converting Trajectory Files to Box Files
--------------------------------------------------------------

"""

import PQAnalysis.config as config

from ._argument_parser import _ArgumentParser
from PQAnalysis.io import traj2box


__outputdoc__ = """

Converts multiple trajectory files to a box file.

Without the --vmd option the output is printed in a data file format.
The first column represents the step starting from 1, the second to fourth column
represent the box vectors a, b, c, the fifth to seventh column represent the box angles.

With the --vmd option the output is printed in a VMD file format. Meaning the output is
in xyz format with 8 particle entries representing the vertices of the box. The comment
line contains the information about the box dimensions a, b and c and the box angles.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.traj2box.html.
"""

__doc__ += __outputdoc__


def main():
    """
    Main function of the traj2box command line tool, which is basically just a wrapper for the traj2box function. For more information on the traj2box function please visit :py:func:`PQAnalysis.io.api.traj2box`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    parser.parse_output_file()

    parser.add_argument('trajectory_file', type=str, nargs='+',
                        help='The trajectory file(s) to be converted.')
    parser.add_argument('--vmd', action='store_true',
                        help='Output in VMD format.')

    args = parser.parse_args()

    traj2box(args.trajectory_file, args.vmd, args.output)
