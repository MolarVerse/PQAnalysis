"""

"""

import argparse

import PQAnalysis.config as config

from ..analysis import RDF, RDFInputFileReader, RDFDataWriter, RDFLogWriter
from ..analysis.rdf.rdfInputFileReader import input_keys_documentation
from ..io import TrajectoryReader, RestartFileReader, MoldescriptorReader
from ..traj import MDEngineFormat
from ..topology import Topology

import cProfile

__outputdoc__ = """"""  # TODO: include here url of the documentation of the output file dynamically
__doc__ += input_keys_documentation


def main():
    parser = argparse.ArgumentParser(description=__outputdoc__)
    parser.add_argument('input_file', type=str, help='The input file.')
    parser.add_argument('-f', '--format', type=str, default="pimd-qmcf",
                        help='The format of the input trajectory. Default is "pimd-qmcf".')
    parser.add_argument('--progress', action='store_true',
                        default=False, help='Show progress bar.')
    args = parser.parse_args()

    engine_format = MDEngineFormat(args.format)

    config.with_progress_bar = args.progress
    _rdf(args.input_file, engine_format)


def _rdf(input_file: str, format: MDEngineFormat):
    input_reader = RDFInputFileReader(input_file)
    input_reader.read()

    restart_reader = RestartFileReader(input_reader.restart_file)
    restart_frame = restart_reader.read()

    moldescriptor_reader = MoldescriptorReader(
        input_reader.moldescriptor_file)
    reference_residues = moldescriptor_reader.read()

    topology = restart_frame.topology
    topology = Topology(atoms=topology.atoms, residue_ids=topology.residue_ids,
                        reference_residues=reference_residues)

    traj_reader = TrajectoryReader(
        input_reader.traj_files, md_format=format, topology=topology)

    rdf = RDF(
        traj=traj_reader,
        reference_species=input_reader.reference_selection,
        target_species=input_reader.target_selection,
        use_full_atom_info=input_reader.use_full_atom_info,
        no_intra_molecular=input_reader.no_intra_molecular,
        n_bins=input_reader.n_bins,
        delta_r=input_reader.delta_r,
        r_max=input_reader.r_max,
        r_min=input_reader.r_min,
    )

    data_writer = RDFDataWriter(input_reader.out_file)
    log_writer = RDFLogWriter(input_reader.log_file)
    log_writer.write_before_run(rdf)

    rdf_data = rdf.run()

    data_writer.write(rdf_data)
    log_writer.write_after_run(rdf)
