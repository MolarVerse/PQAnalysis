from ._argument_parser import _ArgumentParser

__outputdoc__ = """

This is the command line interface for the PQAnalysis package.
"""

from .build_nep_traj import BuildNEPTrajCLI


def main():
    parser = _ArgumentParser(description=__outputdoc__, epilog=None)

    subparsers = parser.add_subparsers(dest='command')

    build_nep_traj_parser = subparsers.add_parser(
        BuildNEPTrajCLI.program_name(), help=BuildNEPTrajCLI.__doc__
    )
    BuildNEPTrajCLI.add_arguments(build_nep_traj_parser)

    args = parser.parse_args()

    if args.command == BuildNEPTrajCLI.program_name():
        BuildNEPTrajCLI.run(args)
