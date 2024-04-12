"""
A module containing the GenFileWriter class
"""

import numpy as np

from .. import BaseWriter, FileWritingMode
from PQAnalysis.atomicSystem import AtomicSystem


class GenFileWriter(BaseWriter):
    """
    A class for writing gen files.
    """

    def __init__(self,
                 filename: str,
                 mode: FileWritingMode | str = "w",
                 ) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the gen file.
        mode : FileWritingMode | str, optional
            The writing mode, by default "w". The following modes are available:
            - "w": write
            - "a": append
            - "o": overwrite
        """

        super().__init__(filename, mode)

    def write(self, system: AtomicSystem) -> None:
        """
        Writes the system to the file.

        Parameters
        ----------
        system : AtomicSystem
            The system to write.
        """

        self.system = system

        self.open()
        self.write_header()
        self.write_coords()
        if not self.system.cell.is_vacuum:
            self.write_box_matrix()

        self.close()

    def write_header(self) -> None:
        """
        Writes the header of the gen file.

        The header consists of two lines. The first line contains the number of atoms and the periodicity of the system. The second line contains the atom names from which the indices of the atoms in the coordinates block can be derived.
        """
        periodicity = "C" if self.system.cell.is_vacuum else "S"

        print({self.system.n_atoms}, {periodicity}, file=self.file)

        element_names = self.system.unique_element_names

        print(" ".join(element_names), file=self.file)

    def write_coords(self) -> None:
        """
        Writes the coordinates of the system.

        The coordinates are written in the following format:
        - The atom index (starting from 1)
        - The element index (starting from 1)
        - The x-coordinate
        - The y-coordinate
        - The z-coordinate
        """
        element_names = self.system.unique_element_names
        for i in range(self.system.n_atoms):
            print(
                "\t",
                i + 1,
                element_names.index(self.system.atoms[i].element_name) + 1,
                self.system.pos[i, 0],
                self.system.pos[i, 1],
                self.system.pos[i, 2],
                file=self.file
            )

    def write_box_matrix(self) -> None:
        """
        Writes the box matrix of the system.

        Before writing the box matrix, the box vector (0, 0, 0) is written.

        The box matrix is written in the following format:
        - The x-coordinate of the first box vector
        - The y-coordinate of the first box vector
        - The z-coordinate of the first box vector
        """
        box_matrix = np.transpose(self.system.cell.box_matrix)

        print(
            0.0,
            0.0,
            0.0,
            file=self.file
        )
        for i in range(3):
            print(
                box_matrix[i, 0],
                box_matrix[i, 1],
                box_matrix[i, 2],
                file=self.file
            )
