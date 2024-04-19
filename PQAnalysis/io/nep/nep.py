import numpy as np
import glob

from beartype.typing import List

from PQAnalysis.utils.units import *
from PQAnalysis.io import BaseWriter, FileWritingMode, OutputFileFormat, read_trajectory_generator, EnergyFileReader
from PQAnalysis.io.virial.api import read_stress_file, read_virial_file
from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.traj import Trajectory, TrajectoryFormat


class NEPWriter(BaseWriter):
    """
    A class to write NEP training and testing files.
    """

    def __init__(self,
                 filename: str | None,
                 mode: FileWritingMode | str = "w",
                 ) -> None:
        """
        Parameters
        ----------
        filename : str
            _description_
        mode : FileWritingMode | str, optional
            _description_, by default "w"
        """
        super().__init__(filename, mode)

    def write_from_files(self,
                         file_prefixes: List[str] | str,
                         use_forces: bool = False,
                         use_stress: bool = False,
                         use_virial: bool = False,
                         xyz_file_extension: str = None,
                         energy_file_extension: str = None,
                         force_file_extension: str = None,
                         stress_file_extension: str = None,
                         virial_file_extension: str = None,
                         ) -> None:

        file_prefixes = list(np.atleast_1d(file_prefixes))
        files = [glob.glob(prefix + ".*") for prefix in file_prefixes]
        files = [file for sublist in files for file in sublist]

        print(files)
        print(OutputFileFormat.get_file_extensions(OutputFileFormat.XYZ))

        # filter all possible xyz files
        xyz_files = [
            file
            for file in files
            if "." + file.split(".")[-1] in OutputFileFormat.get_file_extensions(OutputFileFormat.XYZ)
        ]
        sorted(xyz_files)

        en_files = [
            file
            for file in files
            if "." + file.split(".")[-1] in OutputFileFormat.get_file_extensions(OutputFileFormat.ENERGY)
        ]
        sorted(en_files)

        info_files = [
            file
            for file in files
            if "." + file.split(".")[-1] in OutputFileFormat.get_file_extensions(OutputFileFormat.INFO)
        ]

        if use_forces:
            force_files = [
                file
                for file in files
                if "." + file.split(".")[-1] in OutputFileFormat.get_file_extensions(OutputFileFormat.FORCE)
            ]
            sorted(force_files)

        if use_stress:
            stress_files = [
                file
                for file in files
                if "." + file.split(".")[-1] in OutputFileFormat.get_file_extensions(OutputFileFormat.STRESS)
            ]
            sorted(stress_files)

        if use_virial:
            virial_files = [
                file
                for file in files
                if "." + file.split(".")[-1] in OutputFileFormat.get_file_extensions(OutputFileFormat.VIRIAL)
            ]
            sorted(virial_files)

        if len(xyz_files) == 0:
            raise ValueError(
                "No coordinate files found with the specified file prefixes.")

        if len(en_files) != len(xyz_files):
            raise ValueError(
                "The number of energy files does not match the number of coordinate files.")

        if len(info_files) != len(xyz_files):
            raise ValueError(
                "The number of info files does not match the number of coordinate files.")

        if use_forces:
            if len(force_files) != len(xyz_files):
                raise ValueError(
                    "The number of force files does not match the number of coordinate files.")

        if use_stress:
            if len(stress_files) != len(xyz_files):
                raise ValueError(
                    "The number of stress files does not match the number of coordinate files.")

        if use_virial:
            if len(virial_files) != len(xyz_files):
                raise ValueError(
                    "The number of virial files does not match the number of coordinate files.")

        self.open()

        for i in range(len(xyz_files)):

            # read stress file in format step stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz as n numpy 2D 3x3 arrays
            stress = read_stress_file(stress_files[i]) if use_stress else None
            virial = read_virial_file(virial_files[i]) if use_virial else None
            energy = EnergyFileReader(en_files[i], info_files[i]).read()

            xyz_generator = read_trajectory_generator(
                xyz_files[i], traj_format=TrajectoryFormat.XYZ)
            force_generator = read_trajectory_generator(
                force_files[i], traj_format=TrajectoryFormat.FORCE) if use_forces else None

            for j in range(len(energy.data)):
                system = next(xyz_generator)

                if use_forces:
                    force_system = next(force_generator)
                    system.forces = force_system.forces

                system.energy = energy.qm_energy[j]
                # system.stress = stress if stress is not None else None
                # system.virial = virial if virial is not None else None

                self.write_from_atomic_system(
                    system,
                    use_forces,
                    use_stress,
                    use_virial
                )

        self.close()

    def write_from_trajectory(self,
                              trajectory: Trajectory,
                              use_forces: bool = False,
                              use_stress: bool = False,
                              use_virial: bool = False,
                              ) -> None:

        self.open()
        for frame in trajectory:
            self.write_from_atomic_system(
                frame,
                use_forces,
                use_stress,
                use_virial
            )

        self.close()

    def write_from_atomic_system(self,
                                 system: AtomicSystem,
                                 use_forces: bool = False,
                                 use_stress: bool = False,
                                 use_virial: bool = False,
                                 ) -> None:

        if not system.has_pos:
            raise ValueError(
                "The system does not have coordinates, which are required for NEP trajectory files.")

        if not system.has_energy:
            raise ValueError(
                "The system does not have an energy, which is required for NEP trajectory files.")

        if use_forces and not system.has_forces:
            raise ValueError(
                "The system does not have forces, and they were specified to be written to the NEP trajectory file.")

        if use_stress and use_virial:
            raise ValueError(
                "Both the stress and the virial tensor were specified to be written to the NEP trajectory file. "
                "Only one of them can be written at a time.")

        if use_stress and not system.has_stress:
            raise ValueError(
                "The system does not have a stress tensor, and it was specified to be written to the NEP trajectory file.")

        if use_virial and not system.has_virial:
            raise ValueError(
                "The system does not have a virial tensor, and it was specified to be written to the NEP trajectory file.")

        self.write_header(system, use_forces, use_stress, use_virial)
        self.write_body(system, use_forces)

    def write_header(self,
                     system: AtomicSystem,
                     use_forces: bool = False,
                     use_stress: bool = False,
                     use_virial: bool = False,
                     ) -> None:

        box_matrix = np.transpose(system.cell.box_matrix)

        print(system.n_atoms, file=self.file)

        energy_unit = kcal_per_mole
        energy_conversion = energy_unit.asUnit(eV).asNumber()
        energy *= energy_conversion

        self.file.write(f"energy={energy} ")

        self.file.write("config_type=nep2xyz ")

        self.file.write("lattice=\"")
        for i in range(3):
            for j in range(3):
                self.file.write(f"{box_matrix[i][j]} ")
        self.file.write("\" ")

        virial_unit = kcal_per_mole
        virial_conversion = virial_unit.asUnit(eV).asNumber()
        virial = system.virial * virial_conversion

        if use_virial:
            self.file.write("virial=\"")
            for i in range(3):
                for j in range(3):
                    self.file.write(f"{virial[i][j]} ")
            self.file.write("\" ")

        stress_unit = kcal_per_mole / angstrom**3
        stress_conversion = stress_unit.asUnit(eV / angstrom**3).asNumber()
        stress = system.stress * stress_conversion

        if use_stress:
            self.file.write("stress=\"")
            for i in range(3):
                for j in range(3):
                    self.file.write(f"{stress[i][j]} ")
            self.file.write("\" ")

        self.file.write("properties=species:S:1:pos:R:3")
        if use_forces:
            self.file.write("forces:R:3")
        self.file.write("\n")

    def write_body(self,
                   system: AtomicSystem,
                   use_forces: bool = False,
                   ) -> None:
        """
        Writes the body of the NEP trajectory file.

        Parameters
        ----------
        system : AtomicSystem
            The system to be written to the NEP trajectory file.
        use_forces : bool, optional
            Whether to write the forces to the NEP trajectory file, by default False
        """

        for i in range(system.n_atoms):
            atom = system.atoms[i]

            print(
                f"{atom.symbol} {system.pos[i][0]} {system.pos[i][1]} {system.pos[i][2]}", file=self.file, end=" "
            )

            force_unit = kcal_per_mole / angstrom
            force_conversion = force_unit.asUnit(eV / angstrom).asNumber()
            forces = system.forces * force_conversion

            if use_forces:
                print(
                    f"{forces[i][0]} {forces[i][1]} {forces[i][2]}", file=self.file, end=" "
                )

            print(file=self.file)
