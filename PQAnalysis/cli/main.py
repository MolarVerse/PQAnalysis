"""
A command line interface for the PQAnalysis package.
"""

import sys

from importlib import import_module

from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser

__outputdoc__ = """

This is the command line interface for the PQAnalysis package.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.html."
__epilog__ += "\n"
__epilog__ += "\n"

_COMMANDS = {  # pylint: disable=consider-using-namedtuple-or-dataclass
    "add_molecules": (
        ".add_molecules",
        "AddMoleculesCLI",
        "Add molecules to restart files.",
    ),
    "build_nep_traj": (
        ".build_nep_traj",
        "BuildNEPTrajCLI",
        "Build Neuroevolution Potential training and test trajectories.",
    ),
    "build_spectrum": (
        ".build_spectrum",
        "BuildSpectrumCLI",
        "Broaden stick spectra.",
    ),
    "check_momentum": (
        ".check_momentum",
        "CheckMomentumCLI",
        "Check the total linear momentum.",
    ),
    "continue_input": (
        ".continue_input",
        "ContinueInputCLI",
        "Extend PQ molecular-dynamics input files.",
    ),
    "convert": (
        ".convert",
        "ConvertCLI",
        "Convert analysis-table output files.",
    ),
    "gen2xyz": (
        ".gen2xyz",
        "GEN2XYZCLI",
        "Convert GEN files to XYZ files.",
    ),
    "msd": (".msd", "MSDCLI", "Calculate mean square displacements."),
    "rdf": (".rdf", "RDFCLI", "Calculate radial distribution functions."),
    "rst2xyz": (
        ".rst2xyz",
        "Rst2XYZCLI",
        "Convert restart files to XYZ files.",
    ),
    "xyz2rst": (
        ".xyz2rst",
        "XYZ2RstCLI",
        "Convert XYZ files to restart files.",
    ),
    "traj2box": (
        ".traj2box",
        "Traj2BoxCLI",
        "Convert trajectory files to box files.",
    ),
    "traj2extxyz": (
        ".traj2extxyz",
        "Traj2ExtXYZCLI",
        "Convert PQ trajectory files to extended XYZ.",
    ),
    "traj2qmcfc": (
        ".traj2qmcfc",
        "Traj2QMCFCCLI",
        "Convert PQ trajectory files to QMCFC format.",
    ),
    "vacf": (
        ".vacf",
        "VACFCLI",
        "Calculate velocity autocorrelation functions.",
    ),
    "vibrations": (
        ".vibrations",
        "VibrationsCLI",
        "Calculate molecular vibrations.",
    ),
    "xyz2gen": (
        ".xyz2gen",
        "XYZ2GENCLI",
        "Convert XYZ files to GEN files.",
    ),
}



def _detect_command(arguments: list[str]) -> str | None:
    """Scan root options to find the first positional CLI command."""
    index = 0
    while index < len(arguments):
        argument = arguments[index]

        if argument in {"-h", "--help", "--version"}:
            return None

        if argument == "--logging-level":
            index += 2
            continue

        if argument == "--log-file":
            index += 1
            if (
                index < len(arguments) and
                not arguments[index].startswith("-")
            ):
                index += 1
            continue

        if argument.startswith(("--logging-level=", "--log-file=")):
            index += 1
            continue

        if argument.startswith("-"):
            index += 1
            continue

        return argument if argument in _COMMANDS else None

    return None



def _load_command(command: str):
    module_name, class_name, _ = _COMMANDS[command]
    return getattr(import_module(module_name, __package__), class_name)



def main():
    """
    The main function of the PQAnalysis command line interface.
    """
    selected_command = _detect_command(sys.argv[1:])
    selected_cli = (
        _load_command(selected_command)
        if selected_command in _COMMANDS else None
    )

    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)
    subparsers = parser.add_subparsers(dest='cli_command', )

    for command, (_, _, help_text) in _COMMANDS.items():
        sub_parser = subparsers.add_parser(command, help=help_text)
        if command == selected_command:
            selected_cli.add_arguments(sub_parser)

    args = parser.parse_args()

    if args.cli_command in _COMMANDS:
        command_cli = (
            selected_cli if args.cli_command == selected_command else
            _load_command(args.cli_command)
        )
        command_cli.run(args)
    else:
        parser.print_help()
