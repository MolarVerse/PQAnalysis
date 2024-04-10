import numpy as np

from beartype.typing import List

from PQAnalysis.io import read_restart_file, TrajectoryReader, RestartFileWriter
from PQAnalysis.io.formats import OutputFileFormat
from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.types import PositiveInt, PositiveReal, Np1DNumberArray


class AddMolecule:
    """
    A class for adding a molecule to a restart file.

    The class reads a restart file and a molecule file and adds the molecule to the restart file. The molecule is added by fitting the molecule to the restart file. The fitting is done randomly by rotating the molecule and translating it to a random position. After the fitting, the molecule is added to the restart file. The class can add multiple molecules to the restart file. The class can also add a moldescriptor file to the restart file to keep track of the fitting.
    """

    def __init__(self,
                 restart_file: str,
                 molecule_file: str,
                 output_file: str | None = None,
                 molecule_file_type: OutputFileFormat | str = OutputFileFormat.AUTO,
                 restart_moldescriptor_file: str | None = None,
                 molecule_moldescriptor_file: str | None = None,
                 number_of_additions: PositiveInt = 1,
                 max_iterations: PositiveInt = 100,
                 distance_cutoff: PositiveReal = 1.0,
                 max_displacement: PositiveReal | Np1DNumberArray = 0.1,
                 rotation_angle_step: PositiveInt = 10,
                 md_engine_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                 ) -> None:
        """
        Parameters
        ----------
        restart_file : str
            The filename of the restart file.
        molecule_file : str
            The filename of the molecule file.
        output_file : str | None, optional
            The filename of the output file, by default None
        molecule_file_type : OutputFileFormat | str, optional
            The type of the molecule file, by default OutputFileFormat.AUTO. If the type is AUTO, the type will be inferred from the file extension. If the type is RESTART, a moldescriptor file can be specified.
        restart_moldescriptor_file : str | None, optional
            The filename of the moldescriptor file of the restart file, by default None.
        molecule_moldescriptor_file : str | None, optional
            The filename of the moldescriptor file of the molecule file, by default None. A moldescriptor file can only be specified for restart file types.
        number_of_additions : PositiveInt, optional
            The number of times the molecule should be added to the restart file, by default 1
        max_iterations : PositiveInt, optional
            The maximum number of iterations to try to fit the molecule into the restart file, by default 100
        distance_cutoff : PositiveReal, optional
            The distance cutoff for the fitting, by default 1.0
        max_displacement : PositiveReal | Np1DNumberArray, optional
            The maximum displacement for the fitting, by default 0.1
        rotation_angle_step : PositiveInt, optional
            The angle step for the fitting, by default 10
        md_engine_format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PIMD_QMCF

        Raises
        ------
        ValueError
            If the molecule file type is not RESTART and a moldescriptor file is specified.
        """
        self.restart_file = restart_file
        self.molecule_file = molecule_file
        self.output_file = output_file
        self.molecule_file_type = OutputFileFormat(molecule_file_type)
        self.restart_moldescriptor_file = restart_moldescriptor_file
        self.molecule_moldescriptor_file = molecule_moldescriptor_file
        self.md_engine_format = MDEngineFormat(md_engine_format)
        self.number_of_additions = number_of_additions
        self.max_iterations = max_iterations
        self.distance_cutoff = distance_cutoff
        self.max_displacement = max_displacement
        self.rotation_angle_step = rotation_angle_step

        self.molecule_file_type = OutputFileFormat.infer_format_from_extension(
            self.molecule_file
        )

        if self.molecule_file_type != OutputFileFormat.RESTART and self.molecule_moldescriptor_file is not None:
            raise ValueError(
                "A moldescriptor file can only be specified for restart files."
            )

    def write_restart_file(self) -> None:
        """
        Write the restart file with the added molecule.
        """
        writer = RestartFileWriter(
            self.output_file,
            md_engine_format=self.md_engine_format
        )

        fitted_systems = self.add_molecules()

        self.restart_frame.image()
        lines = writer.get_lines(self.restart_frame, atom_counter=0)
        for i, system in enumerate(fitted_systems):
            system.image()
            lines += writer.get_lines(system, atom_counter=i+1)

        writer.write_lines_to_file(lines)

    def add_molecules(self) -> List[AtomicSystem]:
        """
        Calculates all newly fitted atomic systems

        Returns
        -------
        List[AtomicSystem]
            The fitted atomic systems.
        """
        self.read_files()

        fitted_systems = self.restart_frame.fit_atomic_system(
            system=self.molecule_frame.system,
            number_of_additions=self.number_of_additions,
            max_iterations=self.max_iterations,
            distance_cutoff=self.distance_cutoff,
            max_displacement=self.max_displacement,
            rotation_angle_step=self.rotation_angle_step
        )

        return list(np.atleast_1d(fitted_systems))

    def read_files(self):
        """
        Read the restart and molecule files.
        """
        self.restart_frame = read_restart_file(
            self.restart_file,
            self.restart_moldescriptor_file,
            md_engine_format=self.md_engine_format
        )

        self.molecule_frame = self.read_molecule_file()

    def read_molecule_file(self):
        """
        Read the molecule file.

        Returns
        -------
        Frame
            The frame of the molecule file.
        """
        if self.molecule_file_type == OutputFileFormat.RESTART:

            molecule_frame = read_restart_file(
                self.molecule_file,
                self.molecule_moldescriptor_file,
                md_engine_format=self.md_engine_format
            )

        else:

            frame_generator = TrajectoryReader(
                self.molecule_file,
            ).frame_generator()

            molecule_frame = next(frame_generator)

        return molecule_frame
