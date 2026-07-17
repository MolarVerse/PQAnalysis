"""
.. _cli.build_spectrum:

Command Line Tool for Broadening Stick Spectra
==============================================


"""

from PQAnalysis.analysis.spectrum_broadening import build_spectrum
from PQAnalysis.config import code_base_url

from ._argument_parser import _ArgumentParser
from ._cli_base import CLIBase

__outputdoc__ = """

This command line tool can be used to broaden a stick spectrum
(a two-column file with wavenumbers in cm^-1 and intensities) with
a Gaussian (or Lorentzian) kernel on a regular wavenumber grid.

The broadening uses the peak-height convention, i.e. the broadened
profile of a single stick reaches exactly the stick intensity at the
stick position and no area normalization is applied. The default
Gaussian exponent alpha of 0.0025 cm^-2 corresponds to a full width
at half maximum of about 33.3 cm^-1. Alternatively, the width can be
given directly as a full width at half maximum via --fwhm.
"""

__epilog__ = "\n"
__epilog__ += "For more information on the command line options of this tool please visit "
__epilog__ += f"{code_base_url}PQAnalysis.cli.build_spectrum.html."
__epilog__ += "\n"
__epilog__ += "\n"

__doc__ += __outputdoc__



class BuildSpectrumCLI(CLIBase):

    """
    Command Line Tool for Broadening Stick Spectra
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
        return 'build_spectrum'

    @classmethod
    def add_arguments(cls, parser: _ArgumentParser) -> None:
        """
        Adds the arguments to the parser.

        Parameters
        ----------
        parser : _ArgumentParser
            The parser to which the arguments should be added.
        """
        parser.parse_output_file()

        parser.add_argument(
            'input_file',
            type=str,
            help=(
                'The stick spectrum file to broaden. It must contain '
                'one line per stick with the wavenumber in cm^-1 in '
                'the first column and the intensity in the second column.'
            )
        )

        width_group = parser.add_mutually_exclusive_group()

        width_group.add_argument(
            '--alpha',
            type=float,
            default=None,
            help=(
                'The Gaussian exponent alpha in cm^-2. If neither '
                '--alpha nor --fwhm is given, 0.0025 cm^-2 is used.'
            )
        )

        width_group.add_argument(
            '--fwhm',
            type=float,
            default=None,
            help=(
                'The full width at half maximum in cm^-1 as an '
                'alternative way to specify the broadening width.'
            )
        )

        parser.add_argument(
            '--min',
            dest='wavenumber_min',
            type=float,
            default=10.0,
            help='The first grid point in cm^-1.'
        )

        parser.add_argument(
            '--max',
            dest='wavenumber_max',
            type=float,
            default=4000.0,
            help='The exclusive upper bound of the grid in cm^-1.'
        )

        parser.add_argument(
            '--step',
            dest='wavenumber_step',
            type=float,
            default=0.25,
            help='The grid spacing in cm^-1.'
        )

        parser.add_argument(
            '--lorentzian',
            action='store_true',
            default=False,
            help=(
                'Use a Lorentzian kernel instead of a Gaussian kernel. '
                'The Lorentzian width is chosen such that it has the '
                'same full width at half maximum as the corresponding '
                'Gaussian kernel.'
            )
        )

        parser.parse_mode()

    @classmethod
    def run(cls, args):
        """
        Runs the command line tool.

        Parameters
        ----------
        args : argparse.Namespace
            The arguments parsed by the parser.
        """
        build_spectrum(
            input_file=args.input_file,
            output=args.output,
            alpha=args.alpha,
            fwhm=args.fwhm,
            wavenumber_min=args.wavenumber_min,
            wavenumber_max=args.wavenumber_max,
            wavenumber_step=args.wavenumber_step,
            kernel='lorentzian' if args.lorentzian else 'gaussian',
            mode=args.mode,
        )



def main():
    """
    Main function of the build_spectrum command line tool, which is
    basically just a wrapper for the build_spectrum function. For more
    information on the build_spectrum function please visit
    :py:func:`PQAnalysis.analysis.spectrum_broadening.api.build_spectrum`.
    """
    parser = _ArgumentParser(description=__outputdoc__, epilog=__epilog__)

    BuildSpectrumCLI.add_arguments(parser)

    args = parser.parse_args()

    BuildSpectrumCLI.run(args)
