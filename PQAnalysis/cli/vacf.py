"""
.. _cli.vacf:

Command Line Tool for VACF Analysis
===================================

"""

from PQAnalysis.analysis.vacf import vacf
from PQAnalysis.analysis.vacf.vacf_input_file_reader import input_keys_documentation
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to calculate the
normalized velocity auto-correlation function (VACF)
of given velocity trajectory file(s) using multiple
time origins on a sliding window. In the charge-flux
mode the velocities are weighted with static or
time-dependent atomic partial charges. Optionally,
the legacy cosine-transform spectrum of the VACF is
calculated. This is an input file based tool, so that
the input file can be used to specify the parameters
of the VACF calculation.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.vacf.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__
__doc__ += "For more information on the general the "
__doc__ += "VACF analysis and its input file options "
__doc__ += "please visit "
__doc__ += ":py:class:`PQAnalysis.analysis.vacf.vacf.VACF` "
__doc__ += "and :py:mod:`PQAnalysis.analysis.vacf.vacf_input_file_reader`\n"
__doc__ += input_keys_documentation



class VACFCLI(CLIBase):

    """
    Command Line Tool for VACF Analysis
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
        return 'vacf'

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
        vacf(args.input_file, args.engine)



def main():
    """
    The main function of the VACF analysis command line tool,
    which is basically just a wrapper for the vacf function.
    For more information on the vacf function please
    visit :py:func:`PQAnalysis.analysis.vacf.api.vacf`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    VACFCLI.add_arguments(parser)

    args = parser.parse_args()

    VACFCLI.run(args)
