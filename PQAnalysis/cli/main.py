"""
A command line interface for the PQAnalysis package.
"""

from PQAnalysis.config import code_base_url

from .xyz2gen import XYZ2GENCLI
from .traj2qmcfc import Traj2QMCFCCLI
from .traj2box import Traj2BoxCLI
from .rst2xyz import Rst2XYZCLI
from .rdf import RDFCLI
from .gen2xyz import GEN2XYZCLI
from .continue_input import ContinueInputCLI
from .add_molecules import AddMoleculesCLI
from .build_nep_traj import BuildNEPTrajCLI
from ._argument_parser import _ArgumentParser

__outputdoc__ = """

This is the command line interface for the PQAnalysis package.
"""

__epilog__ = "\n"
__epilog__ += "For more information on required and optional input file keys please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.html."
__epilog__ += "\n"
__epilog__ += "\n"



def main():
    """
    The main function of the PQAnalysis command line interface.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    subparsers = parser.add_subparsers(dest='cli_command', )

    sub_parser_dict = {
        AddMoleculesCLI.program_name(): AddMoleculesCLI,
        BuildNEPTrajCLI.program_name(): BuildNEPTrajCLI,
        ContinueInputCLI.program_name(): ContinueInputCLI,
        GEN2XYZCLI.program_name(): GEN2XYZCLI,
        RDFCLI.program_name(): RDFCLI,
        Rst2XYZCLI.program_name(): Rst2XYZCLI,
        Traj2BoxCLI.program_name(): Traj2BoxCLI,
        Traj2QMCFCCLI.program_name(): Traj2QMCFCCLI,
        XYZ2GENCLI.program_name(): XYZ2GENCLI,
    }

    for key, value in sub_parser_dict.items():
        sub_parser = subparsers.add_parser(key, help=value.__doc__)
        value.add_arguments(sub_parser)

    args = parser.parse_args()

    if args.cli_command in sub_parser_dict:
        sub_parser_dict[args.cli_command].run(args)
    else:
        parser.print_help()
