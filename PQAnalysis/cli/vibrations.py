"""
.. _cli.vibrations:

Command Line Tool for Vibrational Analysis
==========================================
"""

from PQAnalysis.analysis.vibrational import vibrations
from PQAnalysis.analysis.vibrational.vibrational_input_file_reader import (
    input_keys_documentation,
)
from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool calculates vibrational frequencies, force constants,
reduced masses, normal modes and optional IR intensities from a structure file
and a Hessian file.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.vibrations.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__
__doc__ += input_keys_documentation



class VibrationsCLI(CLIBase):

    """
    Command Line Tool for Vibrational Analysis
    """

    @classmethod
    def program_name(cls) -> str:
        """
        Returns the name of the program.
        """
        return "vibrations"

    @classmethod
    def add_arguments(cls, parser: _ArgumentParser) -> None:
        """
        Adds the arguments to the parser.
        """
        parser.parse_input_file()

    @classmethod
    def run(cls, args) -> None:
        """
        Runs the command line tool.
        """
        vibrations(args.input_file)



def main():
    """
    Main function for the standalone vibrational analysis CLI.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    VibrationsCLI.add_arguments(parser)

    args = parser.parse_args()

    VibrationsCLI.run(args)
