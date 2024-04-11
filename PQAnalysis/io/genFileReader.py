import numpy as np

from beartype.typing import List

from .exceptions import GenFileReaderError
from PQAnalysis.io import BaseReader
from PQAnalysis.types import PositiveInt, Np2DNumberArray
from PQAnalysis.core import Cell, Atom
from PQAnalysis.traj import Frame
from PQAnalysis.atomicSystem import AtomicSystem


class GenFileReader(BaseReader):
    def __init__(self, filename: str) -> None:
        super().__init__(filename)

    def read(self) -> Frame:

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

    def read_header(self, header: List[str]) -> PositiveInt, bool, List[str]:
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

    def read_coords(self, lines: List[str]) -> Np2DNumberArray:
        coords = np.zeros((self.n_atoms, 3))
        ids = np.zeros(self.n_atoms, dtype=int)

        for i, line in enumerate(lines):
            coords[i] = np.array(line.split()[2:], dtype=float)
            ids[i] = int(line.split()[1])

        return coords, ids

    def read_cell(self, lines: List[str]) -> Cell:
        box_matrix = np.zeros((3, 3))
        for i, line in enumerate(lines):
            box_matrix[i] = np.array(line.split(), dtype=float)

        # NOTE: The box matrix is stored in row-major order in the gen file, so we need to transpose it.
        return Cell(box_matrix=np.transpose(box_matrix))
