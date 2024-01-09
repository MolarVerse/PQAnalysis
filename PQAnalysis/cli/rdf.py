import argparse

from ..analysis import RDF, RDFInputFileReader, RDFDataWriter, RDFLogWriter
from ..io import TrajectoryReader, RestartFileReader, MoldescriptorReader
from ..traj import MDEngineFormat
from ..topology import Topology


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file', type=str, help='The input file.')
    parser.add_argument('-f', '--format', type=str, default="pimd-qmcf",
                        help='The format of the input trajectory. Default is "pimd-qmcf".')
    parser.add_argument('--progress', action='store_true',
                        default=False, help='Show progress bar.')
    args = parser.parse_args()

    engine_format = MDEngineFormat(args.format)

    _rdf(args.input_file, engine_format, args.progress)


def _rdf(input_file: str, format: MDEngineFormat, with_progress_bar: bool):
    input_reader = RDFInputFileReader(input_file)
    input_reader.read()

    traj_reader = TrajectoryReader(input_reader.traj_files)
    traj = traj_reader.read(md_format=format)

    restart_reader = RestartFileReader(input_reader.restart_file)
    restart_frame = restart_reader.read()

    moldescriptor_reader = MoldescriptorReader(input_reader.moldescriptor_file)
    reference_residues = moldescriptor_reader.read()

    topology = restart_frame.topology
    topology = Topology(atoms=topology.atoms, residue_ids=topology.residue_ids,
                        reference_residues=reference_residues)

    print(topology.n_atoms)
    print(traj[0].n_atoms)
    print(traj.topology.n_atoms)

    traj.topology = topology

    rdf = RDF(
        traj=traj,
        reference_species=input_reader.reference_selection,
        target_species=input_reader.target_selection,
        use_full_atom_info=input_reader.use_full_atom_info,
        no_intra_molecular=input_reader.no_intra_molecular,
        n_bins=input_reader.n_bins,
        delta_r=input_reader.delta_r,
        r_max=input_reader.r_max,
        r_min=input_reader.r_min,
    )

    RDFLogWriter(input_reader.log_file, rdf).write_before_run()

    rdf_data = rdf.run(with_progress_bar=with_progress_bar)

    RDFLogWriter(input_reader.log_file, rdf).write_after_run()
    RDFDataWriter(input_reader.out_file, rdf_data).write()
