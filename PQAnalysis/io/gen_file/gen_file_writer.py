"""
A module containing the GenFileWriter class
"""

import logging

import numpy as np

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.io.base import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.type_checking import runtime_type_checking



class GenFileWriter(BaseWriter):

    """
    A class for writing gen files.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        mode: FileWritingMode | str = "w",
    ) -> None:
        """
        Parameters
        ----------
        filename : str | None, optional
            The filename of the gen file. If None, the file is not opened,
            by default None.
        mode : FileWritingMode | str, optional
            The writing mode, by default "w". The following modes are available:
            - "w": write
            - "a": append
            - "o": overwrite
        """

        super().__init__(filename, mode)

        self.system = None
        self.periodic = None

    @runtime_type_checking
    def write(
        self,
        system: AtomicSystem,
        periodic: bool | None = None,
    ) -> None:
        """
        Writes the system to the file.

        Parameters
        ----------
        system : AtomicSystem
            The system to write.
        periodic : bool, optional
            The periodicity of the system. If True, the system is considered periodic. 
            If False, the system is considered non-periodic. If None, the periodicity 
            is inferred from the system, by default None.
            
        Raises
        ------
        PQValueError
            If the system is non-periodic and periodic is set to True.
        """

        self.system = system

        if periodic is not None:
            if periodic and self.system.cell.is_vacuum:
                self.logger.error(
                    "Invalid periodicity. The system is non-periodic.",
                    exception=PQValueError,
                )

            if periodic:
                self.periodic = "S"
            else:
                self.periodic = "C"
        else:
            self.periodic = "C" if self.system.cell.is_vacuum else "S"

        self.open()
        self._write_header(self.periodic)
        self._write_coords()

        if self.periodic == "S":
            self._write_box_matrix()

        self.close()

    def _write_header(self, periodicity: str) -> None:
        """
        Writes the header of the gen file.

        The header consists of two lines. The first line contains the number of atoms and 
        the periodicity of the system. The second line contains the atom names from which 
        the indices of the atoms in the coordinates block can be derived.
        """
        print(self.system.n_atoms, periodicity, file=self.file)

        element_names = self.system.unique_element_names

        print(" ".join(element_names), file=self.file)

    def _write_coords(self) -> None:
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

    def _write_box_matrix(self) -> None:
        """
        Writes the box matrix of the system.

        Before writing the box matrix, the box vector (0, 0, 0) is written.

        The box matrix is written in the following format:
        - The x-coordinate of the first box vector
        - The y-coordinate of the first box vector
        - The z-coordinate of the first box vector
        """
        box_matrix = np.transpose(self.system.cell.box_matrix)

        print(0.0, 0.0, 0.0, file=self.file)
        for i in range(3):
            print(
                box_matrix[i, 0],
                box_matrix[i, 1],
                box_matrix[i, 2],
                file=self.file
            )
