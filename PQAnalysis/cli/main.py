"""
A command line interface for the PQAnalysis package.
"""

from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser
from .add_molecules import AddMoleculesCLI
from .build_nep_traj import BuildNEPTrajCLI
from .build_spectrum import BuildSpectrumCLI
from .check_momentum import CheckMomentumCLI
from .continue_input import ContinueInputCLI
from .gen2xyz import GEN2XYZCLI
from .green_kubo import GreenKuboCLI
from .msd import MSDCLI
from .rdf import RDFCLI
from .rst2xyz import Rst2XYZCLI
from .traj2box import Traj2BoxCLI
from .traj2extxyz import Traj2ExtXYZCLI
from .traj2qmcfc import Traj2QMCFCCLI
from .vacf import VACFCLI
from .vibrations import VibrationsCLI
from .xyz2gen import XYZ2GENCLI
from .xyz2rst import XYZ2RstCLI

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
        BuildSpectrumCLI.program_name(): BuildSpectrumCLI,
        CheckMomentumCLI.program_name(): CheckMomentumCLI,
        ContinueInputCLI.program_name(): ContinueInputCLI,
        GEN2XYZCLI.program_name(): GEN2XYZCLI,
        GreenKuboCLI.program_name(): GreenKuboCLI,
        MSDCLI.program_name(): MSDCLI,
        RDFCLI.program_name(): RDFCLI,
        Rst2XYZCLI.program_name(): Rst2XYZCLI,
        XYZ2RstCLI.program_name(): XYZ2RstCLI,
        Traj2BoxCLI.program_name(): Traj2BoxCLI,
        Traj2ExtXYZCLI.program_name(): Traj2ExtXYZCLI,
        Traj2QMCFCCLI.program_name(): Traj2QMCFCCLI,
        VACFCLI.program_name(): VACFCLI,
        VibrationsCLI.program_name(): VibrationsCLI,
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
