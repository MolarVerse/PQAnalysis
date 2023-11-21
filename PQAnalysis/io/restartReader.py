import numpy as np

from beartype.typing import Tuple, List

from .base import BaseReader
from ..core.atomicSystem import AtomicSystem
from ..core.atom import Atom
from ..core.cell import Cell
from ..traj.formats import MDEngineFormat
from ..types import Np1DIntArray, FILE, Np2DNumberArray


class RestartFileReader(BaseReader):
    def __init__(self, filename: str, format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF) -> None:
        super().__init__(filename)

        self.format = MDEngineFormat(format)

    def read(self) -> Tuple[AtomicSystem, Np1DIntArray]:

        cell = Cell()
        atom_lines = []

        with open(self.filename, 'r') as file:
            for line in file.readlines():
                if line.strip().startswith('#'):
                    continue

                if len(line.strip().split()) == 0:
                    continue

                line = line.strip().split()

                if line[0].lower() == "box":
                    cell = self._parse_box(line)
                elif line[0].lower() == "step":
                    continue
                elif line[0].lower() == "chi":
                    continue
                else:
                    atom_lines.append(" ".join(line))

        return self._parse_atoms(atom_lines, cell)

    @classmethod
    def _parse_box(cls, line: List[str]) -> Cell:
        if len(line) == 1:
            return Cell()
        elif len(line) == 4:
            box_lengths = [float(l) for l in line[1:]]
            return Cell(*box_lengths)
        elif len(line) == 7:
            box_lengths = [float(l) for l in line[1:4]]
            box_angles = [float(a) for a in line[4:]]
            return Cell(*box_lengths, *box_angles)
        else:
            raise ValueError(
                f"Invalid number of arguments for box: {len(line)}")

    @classmethod
    def _parse_atoms(cls, lines: List[str], cell: Cell = Cell()) -> Tuple[AtomicSystem, Np1DIntArray]:

        atoms = []
        positions = []
        velocities = []
        forces = []
        mol_types = []

        for line in lines:
            line = line.strip().split()

            if len(line) != 12 and len(line) != 21:
                raise ValueError(
                    f"Invalid number of arguments for atom: {len(line)}")

            atoms.append(Atom(line[0], use_guess_element=False))
            mol_types.append(int(line[2]))
            positions.append(np.array([float(l) for l in line[3:6]]))
            velocities.append(np.array([float(l) for l in line[6:9]]))
            forces.append(np.array([float(l) for l in line[9:12]]))

        if atoms == []:
            raise ValueError("No atoms found in restart file.")

        system = AtomicSystem(atoms=atoms, pos=np.array(positions), vel=np.array(
            velocities), forces=np.array(forces), cell=cell)

        return system, np.array(mol_types)
