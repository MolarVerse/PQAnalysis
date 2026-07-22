"""
.. _cli.msd:

Command Line Tool for MSD Analysis
==================================

"""

from PQAnalysis.analysis.msd import msd
from PQAnalysis.analysis.msd.msd_input_file_reader import input_keys_documentation
from PQAnalysis.config import code_base_url
from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to calculate the
mean square displacement (MSD) of given trajectory
file(s) using multiple time origins on a sliding
window. This is an input file based tool, so that
the input file can be used to specify the parameters
of the MSD calculation.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.msd.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__
__doc__ += "For more information on the general the "
__doc__ += "MSD analysis and its input file options "
__doc__ += "please visit "
__doc__ += ":py:class:`PQAnalysis.analysis.msd.msd.MSD` "
__doc__ += "and :py:mod:`PQAnalysis.analysis.msd.msd_input_file_reader`\n"
__doc__ += input_keys_documentation



class MSDCLI(CLIBase):

    """
    Command Line Tool for MSD Analysis
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
        return 'msd'

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
        msd(args.input_file, args.engine)



def main():
    """
    The main function of the MSD analysis command line tool,
    which is basically just a wrapper for the msd function.
    For more information on the msd function please
    visit :py:func:`PQAnalysis.analysis.msd.api.msd`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    MSDCLI.add_arguments(parser)

    args = parser.parse_args()

    MSDCLI.run(args)
