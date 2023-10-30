import pytest
from _pytest.capture import CaptureFixture

from PQAnalysis.utils.common import print_header
from PQAnalysis._version import __version__


def test_print_header(capsys: CaptureFixture):
    print_header()

    captured = capsys.readouterr()

    line = "                                              *"

    assert captured.err == f"""
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
