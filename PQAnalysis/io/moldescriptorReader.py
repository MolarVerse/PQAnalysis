from beartype.typing import List

from .base import BaseReader
from ..topology import MolType
from ..core.atom import Atom


class MoldescriptorReader(BaseReader):
    def __init__(self, filename: str) -> None:
        super().__init__(filename)

    def read(self) -> List[MolType]:
        with open(self.filename, 'r') as f:
            lines = f.readlines()

            mol_types = []

            counter = 0
            while counter < len(lines):
                line = lines[counter]

                if line.strip().startswith('#'):
                    continue

                if len(line.strip().split()) == 0:
                    continue

                line = line.strip().split()

                if line[0].lower() == "water_type":
                    water_type = int(line[1])
                    counter += 1
                elif line[0].lower() == "ammonia_type":
                    ammonia_type = int(line[1])
                    counter += 1
                else:
                    if len(line) != 3:
                        raise ValueError(
                            "The number of columns in the header of a mol type must be 3.")

                    n_atoms = int(line[1])

                    mol_types.append(self._read_mol_type(
                        lines[counter:counter+n_atoms+1]))

                    counter += n_atoms + 1

        return mol_types

    @classmethod
    def _read_mol_type(cls, lines: List[str]) -> MolType:
        header_line = lines[0].strip().split()
        name = header_line[0]
        total_charge = float(header_line[2])

        elements = []
        atom_types = []
        partial_charges = []

        for line in lines[1:]:
            line = line.strip().split()
            elements.append(Atom(line[0]))
            atom_types.append(int(line[1]))
            partial_charges.append(float(line[2]))

        return MolType(name, total_charge, elements, atom_types, partial_charges)
