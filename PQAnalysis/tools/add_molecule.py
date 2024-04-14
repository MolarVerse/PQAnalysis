import numpy as np
import warnings

from beartype.typing import List

from PQAnalysis.io.formats import OutputFileFormat, FileWritingMode
from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.types import (
    PositiveInt,
    PositiveReal,
    Np1DNumberArray,
)
from PQAnalysis.io import (
    read_restart_file,
    TrajectoryReader,
    RestartFileWriter,
    read_topology_file,
    write_topology_file,
)


def add_molecule(restart_file: str,
                 molecule_file: str,
                 output_file: str | None = None,
                 mode: FileWritingMode | str = "w",
                 molecule_file_type: OutputFileFormat | str = OutputFileFormat.AUTO,
                 restart_moldescriptor_file: str | None = None,
                 molecule_moldescriptor_file: str | None = None,
                 number_of_additions: PositiveInt = 1,
                 max_iterations: PositiveInt = 100,
                 distance_cutoff: PositiveReal = 1.0,
                 max_displacement: PositiveReal | Np1DNumberArray = 0.1,
                 rotation_angle_step: PositiveInt = 10,
                 md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
                 topology_file: str | None = None,
                 topology_file_to_add: str | None = None,
                 topology_file_output: str | None = None,
                 ) -> None:
    """
    Add a molecule to a restart file.

    The function reads a restart file and a molecule file and adds the molecule to the restart file. The molecule is added by fitting the molecule to the restart file. The fitting is done randomly by rotating the molecule and translating it to a random position. After the fitting, the molecule is added to the restart file. The function can add multiple molecules to the restart file. The function can also add a moldescriptor file to the restart file to keep track of the fitting.

    Parameters
    ----------
    restart_file : str
        The filename of the restart file.
    molecule_file : str
        The filename of the molecule file.
    output_file : str | None, optional
        The filename of the output file, by default None
    mode : FileWritingMode | str, optional
        The writing mode, by default "w" The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
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
        The format of the restart file, by default MDEngineFormat.PQ
    topology_file : str | None, optional
        The filename of the topology file, by default None
    topology_file_to_add : str | None, optional
        The filename of the topology file to add, by default None
    topology_file_output : str | None, optional
        The filename of the output topology file, by default None

    Raises
    ------
    ValueError
        If the molecule file type is not RESTART and a moldescriptor file is specified.
    """

    check_topology_args(topology_file, topology_file_to_add,
                        topology_file_output)

    add_molecule = AddMolecule(
        restart_file=restart_file,
        molecule_file=molecule_file,
        output_file=output_file,
        mode=mode,
        molecule_file_type=molecule_file_type,
        restart_moldescriptor_file=restart_moldescriptor_file,
        molecule_moldescriptor_file=molecule_moldescriptor_file,
        number_of_additions=number_of_additions,
        max_iterations=max_iterations,
        distance_cutoff=distance_cutoff,
        max_displacement=max_displacement,
        rotation_angle_step=rotation_angle_step,
        md_engine_format=md_engine_format
    )

    add_molecule.write_restart_file()

    if topology_file is not None:
        add_molecule.extend_topology_file(
            original_shake_file=topology_file,
            extension_shake_file=topology_file_to_add,
            output=topology_file_output
        )


def check_topology_args(topology_file: str | None,
                        topology_file_to_add: str | None,
                        topology_file_output: str | None
                        ) -> None:

    if topology_file is None and topology_file_to_add is None:
        if topology_file_output is not None:
            warnings.warn(
                "The output topology file is specified, but no topology files are given to add. The output topology file will be ignored."
            )
        return
    elif topology_file is None and topology_file_to_add is not None:
        raise ValueError(
            "The topology file must be specified if a topology file to add is given."
        )
    elif topology_file is not None and topology_file_to_add is None:
        raise ValueError(
            "The topology file to add must be specified if a topology file is given."
        )

    if topology_file_output is None:
        raise ValueError(
            "The output topology file must be specified if topology files are given to add. This is a special case where None cannot be treated as stdout as it is already used for the restart file."
        )


class AddMolecule:
    """
    A class for adding a molecule to a restart file.

    The class reads a restart file and a molecule file and adds the molecule to the restart file. The molecule is added by fitting the molecule to the restart file. The fitting is done randomly by rotating the molecule and translating it to a random position. After the fitting, the molecule is added to the restart file. The class can add multiple molecules to the restart file. The class can also add a moldescriptor file to the restart file to keep track of the fitting.
    """

    def __init__(self,
                 restart_file: str,
                 molecule_file: str,
                 output_file: str | None = None,
                 mode: FileWritingMode | str = "w",
                 molecule_file_type: OutputFileFormat | str = OutputFileFormat.AUTO,
                 restart_moldescriptor_file: str | None = None,
                 molecule_moldescriptor_file: str | None = None,
                 number_of_additions: PositiveInt = 1,
                 max_iterations: PositiveInt = 100,
                 distance_cutoff: PositiveReal = 1.0,
                 max_displacement: PositiveReal | Np1DNumberArray = 0.1,
                 rotation_angle_step: PositiveInt = 10,
                 md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
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
        mode : FileWritingMode | str, optional
            The writing mode, by default "w" The following modes are available:
            - "w": write
            - "a": append
            - "o": overwrite
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
            The format of the restart file, by default MDEngineFormat.PQ

        Raises
        ------
        ValueError
            If the molecule file type is not RESTART and a moldescriptor file is specified.
        ValueError
            If the molecule file type is not RESTART or XYZ.
        """
        self.restart_file = restart_file
        self.molecule_file = molecule_file
        self.output_file = output_file
        self.mode = FileWritingMode(mode)
        self.restart_moldescriptor_file = restart_moldescriptor_file
        self.molecule_moldescriptor_file = molecule_moldescriptor_file
        self.md_engine_format = MDEngineFormat(md_engine_format)
        self.number_of_additions = number_of_additions
        self.max_iterations = max_iterations
        self.distance_cutoff = distance_cutoff
        self.max_displacement = max_displacement
        self.rotation_angle_step = rotation_angle_step

        # Initialize the variables, just to make sure they are defined
        self.molecule = None
        self.restart_system = None

        self.molecule_file_type = OutputFileFormat(
            molecule_file_type,
            self.molecule_file
        )

        if self.molecule_file_type != OutputFileFormat.RESTART and self.molecule_moldescriptor_file is not None:
            raise ValueError(
                "A moldescriptor file can only be specified for restart files."
            )

        if self.molecule_file_type not in [OutputFileFormat.RESTART, OutputFileFormat.XYZ]:
            raise ValueError(
                "The molecule file type must be either RESTART or XYZ."
            )

    def write_restart_file(self) -> None:
        """
        Write the restart file with the added molecule.
        """
        writer = RestartFileWriter(
            self.output_file,
            md_engine_format=self.md_engine_format,
            mode=self.mode
        )

        fitted_systems = self.add_molecules()

        self.restart_system.image()
        lines = writer.get_lines(self.restart_system, atom_counter=0)
        for i, system in enumerate(fitted_systems):
            system.image()
            lines += writer._get_atom_lines(system, atom_counter=i+1)

        writer.write_lines_to_file(lines)

    def extend_topology_file(self,
                             original_shake_file: str,
                             extension_shake_file: str,
                             output: str | None = None,
                             ) -> None:
        warnings.warn(
            "Extension of the topology file is only implemented for shake bonds. The extension of general bonded topologies is not implemented yet.", UserWarning)

        original_topology = read_topology_file(original_shake_file)
        new_topology = read_topology_file(extension_shake_file)

        if self.restart_system is None or self.molecule is None:
            raise ValueError(
                "The restart frame and the molecule must be read before extending the topology file. Either call the read_files method or the add_molecules method first."
            )

        original_topology.extend_shake_bonds(
            new_topology.shake_bonds,
            n_atoms=self.restart_system.n_atoms,
            n_extensions=self.number_of_additions,
            n_atoms_per_extension=self.molecule.n_atoms
        )

        write_topology_file(original_topology, output, mode=self.mode)

    def add_molecules(self) -> List[AtomicSystem]:
        """
        Calculates all newly fitted atomic systems

        Returns
        -------
        List[AtomicSystem]
            The fitted atomic systems.
        """
        self.read_files()

        fitted_systems = self.restart_system.fit_atomic_system(
            system=self.molecule,
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
        self.restart_system = read_restart_file(
            self.restart_file,
            self.restart_moldescriptor_file,
            md_engine_format=self.md_engine_format
        )

        self.molecule = self.read_molecule_file()

    def read_molecule_file(self) -> AtomicSystem:
        """
        Read the molecule file.

        Returns
        -------
        AtomicSystem
            The AtomicSystem of the molecule file.
        """
        if self.molecule_file_type == OutputFileFormat.RESTART:

            molecule = read_restart_file(
                self.molecule_file,
                self.molecule_moldescriptor_file,
                md_engine_format=self.md_engine_format
            )

        else:

            frame_generator = TrajectoryReader(
                self.molecule_file,
            ).frame_generator()

            molecule = next(frame_generator)

        return molecule
