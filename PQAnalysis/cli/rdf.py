"""
.. _cli.rdf:

Command Line Tool for RDF Analysis
==================================

"""

import PQAnalysis.config as config

from ._argumentParser import Argumentparser
from PQAnalysis.analysis.rdf import rdf
from PQAnalysis.analysis.rdf.rdfInputFileReader import input_keys_documentation
from PQAnalysis.traj import MDEngineFormat

__outputdoc__ = """

This command line tool can be used to calculate the radial distribution function (RDF) of given trajectory file(s). This is an input file based tool, so that the input file can be used to specify the parameters of the RDF calculation.
"""

__reference__ = f"""
For more information on required and optional input file keys please visit {config.code_base_url}PQAnalysis.cli.rdf.html.
"""

__doc__ += __outputdoc__
__doc__ += input_keys_documentation

__outputdoc__ += __reference__


def main():
    """
    The main function of the RDF analysis command line tool.
    """
    parser = Argumentparser(description=__outputdoc__)
    parser.parse_md_format()
    parser.parse_input_file()
    parser._parse_progress()

    args = parser.parse_args()

    rdf(args.input_file, args.engine_format)
