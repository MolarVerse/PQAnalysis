"""
.. _cli.traj2qmcfc:

Command Line Tool for Converting PIMD-QMCF to QMCFC Trajectory Files
---------------------------------------------------------------------------

"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.io import traj2qmcfc


__outputdoc__ = """
Converts a PIMD-QMCF trajectory to a QMCFC trajectory format output.

Both formats are adapted xyz file formats with the box dimensions and box angles
being placed in the same line after the number of atoms. The QMCFC contains an 
additional dummy 'X' particle as first entry of the coordinates section for
visibility and debugging reasons in conjunction with vmd.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.traj2qmcfc.html.
"""

__doc__ += __outputdoc__


def main():
    """
    Main function of the traj2qmcfc command line tool, which is basically just a wrapper for the traj2qmcfc function. For more information on the traj2qmcfc function please visit :py:func:`PQAnalysis.io.api.traj2qmcfc`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    parser.parse_trajectory_file()
    parser.parse_output_file()

    args = parser.parse_args()

    traj2qmcfc(args.trajectory_file, args.output)
