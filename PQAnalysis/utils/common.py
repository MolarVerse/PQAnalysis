"""
Common things needed by command line scripts.
"""

import sys
import warnings

from .._version import __version__

line = "                                              *"
header = r"""
**************************************************************
*                                                            *
*      ____  ____    ___                __           _       *
*     / __ \/ __ \  /   |  ____  ____ _/ /_  _______(_)____  *
*    / /_/ / / / / / /| | / __ \/ __ `/ / / / / ___/ / ___/  *
*   / ____/ /_/ / / ___ |/ / / / /_/ / / /_/ (__  ) (__  )   *
*  /_/    \___\_\/_/  |_/_/ /_/\__,_/_/\__, /____/_/____/    *
*                                     /____/                 *
*                                                            *
*                                                            *"""

header += f"""
*  authors:    Jakob Gamper, Josef M. Gallmetzer             *
*  version:    {__version__}{line[len(__version__):]}
*                                                            *
**************************************************************
"""


def print_header(file: str = None) -> None:
    """
    A function to print the header of the program.

    This header is printed to stderr so that it can be used
    by all cli scripts without interfering with the output.
    """

    if file is None:
        print(header, file=sys.stderr)
    else:
        print(header, file=file)
