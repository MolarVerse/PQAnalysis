"""
.. _cli.thermal_expansion:

Command Line Tool for linear thermal expansion coefficient Analysis
==================================

"""

from PQAnalysis.analysis.thermal_expansion import thermal_expansion
from PQAnalysis.analysis.thermal_expansion.thermal_expansion_input_file_reader import input_keys_documentation
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to calculate the 
linear thermal expansion coefficient of given lattice parameters or 
volumetric thermal expansion coefficient of given volume data.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.thermal_expansion.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__
__doc__ += "For more information on the general the "
__doc__ += "linear thermal expansion coefficient analysis and its input file options "
__doc__ += "please visit "
__doc__ += ":py:class:`PQAnalysis.analysis.thermal_expansion.thermal_expansion.thermal_expansion.ThermalExpansion` "
__doc__ += "and :py:mod:`PQAnalysis.analysis.thermal_expansion.thermal_input_file_reader`\n"
__doc__ += input_keys_documentation


class ThermalExpansionCLI(CLIBase):

    """
    Command Line Tool for linear or volumetric thermal expansion coefficient Analysis
    """

    @classmethod
    def program_name(cls) -> str:
        """
        Returns the name of the program.

        Returns
        -------
        str
            The name of the program.
        """
        return 'thermal_expansion'

    @classmethod
    def add_arguments(cls, parser: _ArgumentParser) -> None:
        """
        Adds the arguments to the parser.

        Parameters
        ----------
        parser : _ArgumentParser
            The parser to which the arguments should be added.
        """
        parser.parse_input_file()
        parser.parse_engine()

    @classmethod
    def run(cls, args):
        """
        Runs the command line tool.

        Parameters
        ----------
        args : argparse.Namespace
            The arguments parsed by the parser.
        """
        thermal_expansion(args.input_file, args.engine)


def main():
    """
    The main function of the linear thermal expansion coefficient analysis command line tool,
    which is basically just a wrapper for the thermal_expansion function. 
    For more information on the thermal_expansion function please 
    visit :py:func:`PQAnalysis.analysis.thermal_expansion.api.thermal_expansion`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    ThermalExpansionCLI.add_arguments(parser)

    args = parser.parse_args()

    ThermalExpansionCLI.run(args)
