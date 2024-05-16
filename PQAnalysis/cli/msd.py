"""
.. _cli.msd:

Command Line Tool for Self-Diffusion Coefficient Analysis
==================================

"""

from PQAnalysis.analysis.msd import msd
from PQAnalysis.analysis.msd.msd_input_file_reader import input_keys_documentation
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser

__outputdoc__ = """

This command line tool can be used to calculate the self-diffusion coefficient (Diffcalc) of given trajectory file(s). This is an input file based tool, so that the input file can be used to specify the parameters of the self-diffusion coefficient calculation.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {code_base_url}PQAnalysis.cli.msd.html.
"""

__doc__ += __outputdoc__
__doc__ += "For more information on the general the self-diffusion coeffficent analysis and its input file options please visit :py:class:`PQAnalysis.analysis.msd.msd.Diffcalc` and :py:mod:`PQAnalysis.analysis.msd.msdInputFileReader`\n"
__doc__ += input_keys_documentation


def main():
    """
    The main function of the msd analysis command line tool, which is basically just a wrapper for the msd function. For more information on the msd function please visit :py:func:`PQAnalysis.analysis.msd.api.msd`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)
    parser.parse_engine_format()
    parser.parse_input_file()

    args = parser.parse_args()

    msd(args.input_file, args.engine_format)
