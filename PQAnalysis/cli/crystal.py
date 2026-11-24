"""Convert an ordered CIF structure to a PQ restart.

.. _cli.crystal:

Crystal preparation
-------------------
"""

from PQAnalysis.io.crystal import prepare_crystal
from PQAnalysis.io.restart_file.api import write_restart_file
from ._cli_base import CLIBase



class CrystalCLI(CLIBase):

    """Expose crystal preparation on the command line."""

    @classmethod
    def program_name(cls):
        """Return the subcommand name."""
        return "crystal"

    @classmethod
    def add_arguments(cls, parser):
        """Register the CIF path, output path and replication options."""
        parser.add_argument("cif")
        parser.add_argument("output")
        parser.add_argument("--repeat", nargs=3, type=int, default=(1, 1, 1))
        parser.add_argument("--block")

    @classmethod
    def run(cls, args):
        """Prepare the geometry and write a new restart file."""
        try:
            system = prepare_crystal(
                args.cif, tuple(args.repeat), block=args.block
            )
            write_restart_file(system, args.output)
        except (ImportError, OSError, ValueError) as error:
            raise SystemExit(str(error)) from error
