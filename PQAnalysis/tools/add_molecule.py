from PQAnalysis.io import read_restart_file, TrajectoryReader
from PQAnalysis.io.formats import OutputFileFormat


def add_molecule_to_restart_file(restart_file: str,
                                 molecule_file: str,
                                 molecule_file_type: OutputFileFormat | str = OutputFileFormat.AUTO,
                                 restart_moldescriptor_file: str | None = None,
                                 molecule_moldescriptor_file: str | None = None,
                                 ):

    restart_frame = read_restart_file(
        restart_file,
        restart_moldescriptor_file,
        md_engine_format="PIMD_QMCF"  # TODO: add kwargs
    )

    molecule_file_type = OutputFileFormat(molecule_file_type)
    molecule_file_type = OutputFileFormat.infer_format_from_extension(
        molecule_file
    )

    if molecule_file_type != OutputFileFormat.RESTART and molecule_moldescriptor_file is not None:
        raise ValueError(
            "A moldescriptor file can only be specified for restart files."
        )

    if molecule_file_type == OutputFileFormat.RESTART:
        reader = RestartFileReader(
            molecule_file,
            molecule_moldescriptor_file,
        )
        molecule_frame = reader.read_frame()
    else:
        frame_generator = TrajectoryReader(
            molecule_file,
            molecule_file_type
        ).frame_generator()

        molecule_frame = next(frame_generator)

    restart_frame.fit_atomic_system(
        molecule_frame.system
        # TODO: add kwargs
    )
