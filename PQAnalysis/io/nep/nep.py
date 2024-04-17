import numpy as np

from beartype.typing import List

from PQAnalysis.io import BaseWriter, FileWritingMode
from PQAnalysis.atomicSystem import AtomicSystem


class NEPWriter(BaseWriter):
    """
    A class to write NEP training and testing files.
    """

    def __init__(self,
                 filename: str,
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

        self.file.write(f"energy={system.energy} ")

        self.file.write("config_type=nep2xyz ")

        self.file.write("lattice=\"")
        for i in range(3):
            for j in range(3):
                self.file.write(f"{box_matrix[i][j]} ")
        self.file.write("\" ")

        if use_virial:
            self.file.write("virial=\"")
            for i in range(3):
                for j in range(3):
                    self.file.write(f"{system.virial[i][j]} ")
            self.file.write("\" ")

        if use_stress:
            self.file.write("stress=\"")
            for i in range(3):
                for j in range(3):
                    self.file.write(f"{system.stress[i][j]} ")
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

            if use_forces:
                print(
                    f"{system.forces[i][0]} {system.forces[i][1]} {system.forces[i][2]}", file=self.file, end=" "
                )

            print(file=self.file)
