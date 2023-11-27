"""
A module for reading moldescriptor files.

Classes
-------
MoldescriptorReader
    A class for reading moldescriptor files.
"""

import numpy as np

from beartype.typing import List

from . import BaseReader
from .exceptions import MoldescriptorReaderError
from ..topology import MolType
from ..core import Atom


class MoldescriptorReader(BaseReader):
    """
    A class for reading moldescriptor files.

    Inherits from the BaseReader class.

    Attributes
    ----------
    filename : str
        The filename of the moldescriptor file.
    """

    def __init__(self, filename: str) -> None:
        """
        Initializes the MoldescriptorReader with the given parameters.

        Parameters
        ----------
        filename : str
            The filename of the moldescriptor file.
        """
        super().__init__(filename)

    def read(self) -> List[MolType]:
        """
        Reads the moldescriptor file and returns the mol types.

        Returns
        -------
        List[MolType]
            The mol types.

        Raises
        ------
        MoldescriptorReaderError
            If the number of columns in the header of a mol type is not 3.
        """
        with open(self.filename, 'r') as f:
            lines = f.readlines()

            mol_types = []

            counter = 0
            while counter < len(lines):
                line = lines[counter]

                if line.strip().startswith('#'):
                    counter += 1
                    continue

                if len(line.strip().split()) == 0:
                    counter += 1
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
                        raise MoldescriptorReaderError(
                            "The number of columns in the header of a mol type must be 3.")

                    n_atoms = int(line[1])

                    mol_types.append(self._read_mol_type(
                        lines[counter:counter+n_atoms+1], len(mol_types) + 1))

                    counter += n_atoms + 1

        return mol_types

    @classmethod
    def _read_mol_type(cls, lines: List[str], mol_type_id: int) -> MolType:
        """
        Parses a mol type from the given lines.

        Parameters
        ----------
        lines : List[str]
            The header line and the body lines of the mol type.
        mol_type_id : int
            The id of the mol type.

        Returns
        -------
        MolType
            The mol type.

        Raises
        ------
        MoldescriptorReaderError
            If the number of columns in the body lines is not 3 or 4.
        """
        header_line = lines[0].strip().split()
        name = header_line[0]
        total_charge = float(header_line[2])

        elements = []
        atom_types = []
        partial_charges = []

        for line in lines[1:]:
            line = line.strip().split()

            if len(line) != 3 and len(line) != 4:
                raise MoldescriptorReaderError(
                    "The number of columns in the body of a mol type must be 3 or 4.")

            elements.append(Atom(line[0]))
            atom_types.append(int(line[1]))
            partial_charges.append(float(line[2]))

        return MolType(name, mol_type_id, total_charge, elements, np.array(atom_types), np.array(partial_charges))
