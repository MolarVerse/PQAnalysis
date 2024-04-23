"""
.. _cli.diffcalc:

Command Line Tool for Self-Diffusion Coefficient Analysis
==================================

"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.analysis.diffcalc import diffcalc
from PQAnalysis.analysis.diffcalc.diffcalcInputFileReader import input_keys_documentation
from PQAnalysis.traj import MDEngineFormat

__outputdoc__ = """

This command line tool can be used to calculate the self-diffusion coefficient (Diffcalc) of given trajectory file(s). This is an input file based tool, so that the input file can be used to specify the parameters of the self-diffusion coefficient calculation.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.diffcalc.html.
"""

__doc__ += __outputdoc__
__doc__ += "For more information on the general the self-diffusion coeffficent analysis and its input file options please visit :py:class:`PQAnalysis.analysis.diffcalc.diffcalc.Diffcalc` and :py:mod:`PQAnalysis.analysis.diffcalc.diffcalcInputFileReader`\n"
__doc__ += input_keys_documentation


def main():
    """
    The main function of the diffcalc analysis command line tool, which is basically just a wrapper for the diffcalc function. For more information on the diffcalc function please visit :py:func:`PQAnalysis.analysis.diffcalc.api.diffcalc`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)
    parser.parse_engine_format()
    parser.parse_input_file()

    args = parser.parse_args()

    diffcalc(args.input_file, args.engine_format)
