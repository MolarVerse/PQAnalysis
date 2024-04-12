"""
A module for reading moldescriptor files.
"""

import numpy as np

from beartype.typing import List

from . import BaseReader
from .exceptions import MoldescriptorReaderError
from PQAnalysis.core import Residue, Residues, Element


class MoldescriptorReader(BaseReader):
    """
    This is a class to read moldescriptor files. Moldescriptor files are used by the PIMD-QMCF and QMCFC MD engines to store the mol types of a system. The moldescriptor file is a text file that contains the mol types of a system. The mol types are defined by the elements, atom types and partial charges of the atoms of the mol type. 
    For more information of how a moldescriptor file should look like please see the documentation of the `PIMD-QMCF <https://molarverse.github.io/pimd_qmcf>`_ code.

    Calling the read method returns a list of mol types. Each mol type is a :py:class:`~PQAnalysis.core.residue.Residue` object.
    """

    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the moldescriptor file.
        """
        super().__init__(filename)

    def read(self) -> Residues:
        """
        The read method reads a moldescriptor file and returns a list of mol types. Each mol type is a :py:class:`~PQAnalysis.core.residue.Residue` object.

        Returns
        -------
        Residues
            The residues (mol types) read from the moldescriptor file.

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

                # remove everything from line after first # to avoid comments
                splitted_line = line.split('#')[0]
                splitted_line = splitted_line.strip().split()

                if splitted_line[0].lower() == "water_type":
                    water_type = int(splitted_line[1])  # TODO:
                    counter += 1
                elif splitted_line[0].lower() == "ammonia_type":
                    ammonia_type = int(splitted_line[1])  # TODO:
                    counter += 1
                else:
                    if len(splitted_line) != 3:
                        line = line.replace('\n', '')
                        raise MoldescriptorReaderError(
                            f"The number of columns in the header of a mol type must be 3.\nGot {len(splitted_line)} columns instead in text: '{line}'\n"
                        )

                    n_atoms = int(splitted_line[1])

                    mol_types.append(self._read_mol_type(
                        lines[counter:counter+n_atoms+1], len(mol_types) + 1))

                    counter += n_atoms + 1

        return mol_types

    @classmethod
    def _read_mol_type(cls, lines: List[str], mol_type_id: int) -> Residue:
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
            splitted_line = line.split('#')[0]
            splitted_line = splitted_line.strip().split()

            if len(splitted_line) != 3 and len(splitted_line) != 4:
                line = line.replace('\n', '')
                raise MoldescriptorReaderError(
                    f"The number of columns in the body of a mol type must be 3 or 4.\nGot {len(splitted_line)} columns instead in text: '{line}'\n"
                )

            elements.append(Element(splitted_line[0]))
            atom_types.append(int(splitted_line[1]))
            partial_charges.append(float(splitted_line[2]))

        return Residue(name, mol_type_id, total_charge, elements, np.array(atom_types), np.array(partial_charges))
