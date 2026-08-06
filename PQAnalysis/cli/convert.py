"""
Command line interface for converting analysis-table output.
"""

from PQAnalysis.analysis.output import convert_analysis_output
from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

Convert a native PQAnalysis, CSV or TSV analysis table into one or
more output files. The output extension selects CSV, TSV, XVG or the
native self-describing text format.
"""

__epilog__ = "\n"
__epilog__ += "For analysis output column definitions please visit "
__epilog__ += f"{code_base_url}userGuide/analysisOutputFiles.html."
__epilog__ += "\n\n"



class ConvertCLI(CLIBase):

    """
    Convert analysis-table output files.
    """

    @classmethod
    def program_name(cls) -> str:
        """
        Return the pqanalysis subcommand name.
        """
        return "convert"

    @classmethod
    def add_arguments(cls, parser: _ArgumentParser) -> None:
        """
        Add converter arguments to a parser.
        """
        parser.parse_input_file()
        parser.add_argument(
            "-o",
            "--output",
            dest="output_files",
            action="append",
            required=True,
            metavar="FILE",
            help=(
                "Output file. Repeat for multiple outputs. The .csv, .tsv "
                "and .xvg extensions select those formats; other extensions "
                "use native PQAnalysis text."
            ),
        )
        parser.add_argument(
            "--x",
            dest="x_field",
            default=None,
            metavar="FIELD",
            help="Override the x-axis field used for XVG output.",
        )
        parser.add_argument(
            "--y",
            dest="y_fields",
            action="append",
            default=None,
            metavar="FIELD",
            help="Override an XVG y-axis field. Repeat for multiple sets.",
        )
        parser.parse_mode()

    @classmethod
    def run(cls, args) -> None:
        """
        Run the requested conversions.
        """
        convert_analysis_output(
            input_file=args.input_file,
            output_files=args.output_files,
            x_field=args.x_field,
            y_fields=args.y_fields,
            mode=args.mode,
        )



def main():
    """
    Run the standalone analysis-table converter CLI.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)
    ConvertCLI.add_arguments(parser)
    args = parser.parse_args()
    ConvertCLI.run(args)
