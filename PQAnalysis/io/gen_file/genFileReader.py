"""
A module containing the GenFileReader class
"""

import numpy as np

from beartype.typing import List, Tuple

from .exceptions import GenFileReaderError
from PQAnalysis.io import BaseReader
from PQAnalysis.types import PositiveInt, Np2DNumberArray, Np1DIntArray
from PQAnalysis.core import Cell, Atom
from PQAnalysis.traj import Frame
from PQAnalysis.atomicSystem import AtomicSystem


class GenFileReader(BaseReader):
    """
    A class for reading gen files.
    """

    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the gen file.
        """
        super().__init__(filename)

    def read(self) -> Frame:
        """
        Reads the gen file and returns a Frame object.

        Returns
        -------
        Frame
            The Frame object including the AtomicSystem and the Cell.
        """
        with open(self.filename, 'r') as file:
            lines = file.read_lines()

            self.n_atoms, is_periodic, atom_names = self.read_header(lines[:2])

            coords, ids = self.read_coords(lines[2:2+self.n_atoms])

            if is_periodic:
                cell = self.read_cell(lines[2+self.n_atoms:2+self.n_atoms+3])
            else:
                cell = Cell()

            atoms = [Atom(atom_names[id - 1]) for id in ids]

            return Frame(AtomicSystem(atoms=atoms, pos=coords, cell=cell))

    def read_header(self, header: List[str]) -> Tuple[PositiveInt, bool, List[str]]:
        """
        Reads the header of the gen file.

        The header consists of two lines. The first line contains the number of atoms and the periodicity of the system. The second line contains the atom names from which the indices of the atoms in the coordinates block can be derived.

        Parameters
        ----------
        header : List[str]
            The header of the gen file.

        Returns
        -------
        Tuple[PositiveInt, bool, List[str]]
            The number of atoms, the periodicity of the system and the atom names.

        Raises
        ------
        GenFileReaderError
            If the periodicity is not "c" or "s".
        """
        n_atoms, periodicity = header[0].split()
        n_atoms = int(n_atoms)

        if periodicity.lower() == "c":
            is_periodic = False
        elif periodicity.lower() == "s":
            is_periodic = True
        else:
            raise GenFileReaderError(
                f"Unknown periodicity: {periodicity} in line {header[0:-1]}"
            )

        atom_names = [name.lower() for name in header[1].split()]

        return n_atoms, is_periodic, atom_names

    def read_coords(self, lines: List[str]) -> Tuple[Np2DNumberArray, Np1DIntArray]:
        """
        Reads the coordinates block of the gen file.

        One line corresponds to one atom and contains a counter index, the atom index and the x, y and z coordinates.

        Parameters
        ----------
        lines : List[str]
            The coordinates block of the gen file.

        Returns
        -------
        Tuple[Np2DNumberArray, Np1DIntArray]
            The coordinates and atom indices.
        """
        coords = np.zeros((self.n_atoms, 3))
        ids = np.zeros(self.n_atoms, dtype=int)

        for i, line in enumerate(lines):
            coords[i] = np.array(line.split()[2:], dtype=float)
            ids[i] = int(line.split()[1])

        return coords, ids

    def read_cell(self, lines: List[str]) -> Cell:
        """
        Reads the cell block of the gen file.

        The cell block contains the box matrix in row-major order.

        Parameters
        ----------
        lines : List[str]
            The cell block of the gen file.

        Returns
        -------
        Cell
            The cell object.
        """
        box_matrix = np.zeros((3, 3))
        for i, line in enumerate(lines):
            box_matrix[i] = np.array(line.split(), dtype=float)

        # NOTE: The box matrix is stored in row-major order in the gen file, so we need to transpose it.
        return Cell(box_matrix=np.transpose(box_matrix))
