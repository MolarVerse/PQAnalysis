"""
.. _cli.rdf:

Command Line Tool for RDF Analysis
==================================

"""

import PQAnalysis.config as config

from ._argumentParser import _ArgumentParser
from PQAnalysis.analysis.rdf import rdf
from PQAnalysis.analysis.rdf.rdfInputFileReader import input_keys_documentation
from PQAnalysis.traj import MDEngineFormat

__outputdoc__ = """

This command line tool can be used to calculate the radial distribution function (RDF) of given trajectory file(s). This is an input file based tool, so that the input file can be used to specify the parameters of the RDF calculation.
"""

__epilog__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.rdf.html.
"""

__doc__ += __outputdoc__
__doc__ += "For more information on the general the RDF analysis and its input file options please visit :py:class:`PQAnalysis.analysis.rdf.rdf.RDF` and :py:mod:`PQAnalysis.analysis.rdf.rdfInputFileReader`\n"
__doc__ += input_keys_documentation


def main():
    """
    The main function of the RDF analysis command line tool, which is basically just a wrapper for the rdf function. For more information on the rdf function please visit :py:func:`PQAnalysis.analysis.rdf.api.rdf`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)
    parser.parse_md_format()
    parser.parse_input_file()
    parser._parse_progress()

    args = parser.parse_args()

    rdf(args.input_file, args.engine_format)
