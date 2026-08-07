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

The main table contains signed wavenumbers, optional IR intensities,
force constants and reduced masses. Optional files contain the
normal-mode matrix or XYZ mode representations.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.vibrations.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__
__doc__ += (
    "For column and mode-file definitions see "
    ":ref:`vibrational-analysis output <analysis-output-vibrations>`.\n"
)
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
        parser.parse_export_files()

    @classmethod
    def run(cls, args) -> None:
        """
        Runs the command line tool.
        """
        export_kwargs = {}
        if args.export_files is not None:
            export_kwargs['export_files'] = args.export_files

        vibrations(args.input_file, **export_kwargs)



def main():
    """
    Main function for the standalone vibrational analysis CLI.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    VibrationsCLI.add_arguments(parser)

    args = parser.parse_args()

    VibrationsCLI.run(args)
