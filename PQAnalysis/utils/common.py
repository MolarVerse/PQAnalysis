"""
Common things needed by command line scripts.
"""

import sys

from .._version import __version__


def print_header() -> None:
    """
    A function to print the header of the program.

    This header is printed to stderr so that it can be used
    by all cli scripts without interfering with the output.
    """

    line = "                                              *"

    header = f"""
**************************************************************
*                                                            *
*      ____  ____    ___                __           _       *
*     / __ \/ __ \  /   |  ____  ____ _/ /_  _______(_)____  *
*    / /_/ / / / / / /| | / __ \/ __ `/ / / / / ___/ / ___/  *
*   / ____/ /_/ / / ___ |/ / / / /_/ / / /_/ (__  ) (__  )   *
*  /_/    \___\_\/_/  |_/_/ /_/\__,_/_/\__, /____/_/____/    *
*                                     /____/                 *
*                                                            *
*                                                            *
*  authors:    Jakob Gamper, Josef M. Gallmetzer             *
*  version:    {__version__}{line[len(__version__):]}
*                                                            *
**************************************************************
"""
    print(header, file=sys.stderr)
