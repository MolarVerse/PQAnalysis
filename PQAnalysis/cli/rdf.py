import argparse

from ..analysis import RadialDistributionFunction, RDFInputFileReader, RDFOutputFileWriter
from ..io import TrajectoryReader
from ..traj import MDEngineFormat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file', type=str, help='The input file.')
    parser.add_argument('-f', '--format', type=str, default="pimd-qmcf",
                        help='The format of the input trajectory. Default is "pimd-qmcf".')
    parser.add_argument('--progress', action='store_true',
                        default=False, help='Show progress bar.')
    args = parser.parse_args()

    format = MDEngineFormat(args.format)

    _rdf(args.input_file, format, args.progress)


def _rdf(input_file: str, format: MDEngineFormat, with_progress_bar: bool):
    reader = RDFInputFileReader(input_file)
    reader.read()

    reader = TrajectoryReader(reader.traj_files)
    traj = reader.read(format=format)

    rdf = RadialDistributionFunction(
        traj=traj,
        reference_species=reader.reference_selection,
        reference_indices=reader.reference_selection,
        use_full_atom_info=reader.use_full_atom_info_for_selection,
        n_bins=reader.n_bins,
        delta_r=reader.delta_r,
        r_max=reader.r_max,
        r_min=reader.r_min,
    )

    rdf_data = rdf.run(with_progress_bar=with_progress_bar)

    writer = RDFOutputFileWriter(
        reader.out_file, reader.log_file, rdf_data, rdf)
    writer.write()
