"""
.. _cli.green_kubo:

Command Line Tool for Green-Kubo Analysis
=========================================

"""

from PQAnalysis.analysis.green_kubo import green_kubo
from PQAnalysis.analysis.green_kubo.green_kubo_input_file_reader import (
    input_keys_documentation,
)
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to calculate a
Green-Kubo transport coefficient (currently the
self-diffusion coefficient) of given velocity
trajectory file(s). The diffusion coefficient is
obtained from the time integral of the un-normalized
velocity auto-correlation function. This is an input
file based tool, so that the input file can be used to
specify the parameters of the Green-Kubo calculation.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.green_kubo.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__
__doc__ += "For more information on the general the "
__doc__ += "Green-Kubo analysis and its input file options "
__doc__ += "please visit "
__doc__ += ":py:class:`PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo` "
__doc__ += "and :py:mod:`PQAnalysis.analysis.green_kubo.green_kubo_input_file_reader`\n"
__doc__ += input_keys_documentation



class GreenKuboCLI(CLIBase):

    """
    Command Line Tool for Green-Kubo Analysis
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
        return 'green_kubo'

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
        green_kubo(args.input_file, args.engine)



def main():
    """
    The main function of the Green-Kubo analysis command line tool,
    which is basically just a wrapper for the green_kubo function.
    For more information on the green_kubo function please
    visit :py:func:`PQAnalysis.analysis.green_kubo.api.green_kubo`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    GreenKuboCLI.add_arguments(parser)

    args = parser.parse_args()

    GreenKuboCLI.run(args)
