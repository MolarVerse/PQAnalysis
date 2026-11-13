"""
A module for reading moldescriptor files.
"""

import logging

import numpy as np

from beartype.typing import List

from PQAnalysis.core import Residue, Residues, Element
from PQAnalysis.io.base import BaseReader
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import MoldescriptorReaderError



class MoldescriptorReader(BaseReader):

    """
    This is a class to read moldescriptor files. Moldescriptor files
    are used by the PQ and QMCFC MD engines to store the mol types of
    a system. The moldescriptor file is a text file that contains the
    mol types of a system. The mol types are defined by the elements,
    atom types and partial charges of the atoms of the mol type.
    For more information of how a moldescriptor file should look like
    please see the documentation of the `PQ <https://molarverse.github.io/PQ>`_ code.

    Calling the read method returns a list of mol types. Each mol
    type is a :py:class:`~PQAnalysis.core.residue.Residue` object.

    The lines ``WATER_TYPE n`` and ``AMMONIA_TYPE n`` of the file name the
    mol type (by its number, i.e. the ``id`` of the residue) of the water and
    the ammonia molecules. They are available after reading as the attributes
    ``water_type`` and ``ammonia_type`` (None if the line is missing). They
    only point to a mol type of the file, the mol type itself still needs its
    own entry.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the moldescriptor file.
        """
        super().__init__(filename)

        #: int | None: The mol type number of the water molecules.
        self.water_type = None

        #: int | None: The mol type number of the ammonia molecules.
        self.ammonia_type = None

    def read(self) -> Residues:
        """
        The read method reads a moldescriptor file and returns a
        list of mol types. Each mol type is a :py:class:`~PQAnalysis.core.residue.Residue`
        object.

        Returns
        -------
        Residues
            The residues (mol types) read from the moldescriptor file.

        Raises
        ------
        MoldescriptorReaderError
            If the number of columns in the header of a mol type is not 3.
        MoldescriptorReaderError
            If a mol type declares more atoms than atom lines follow
            before the end of the file.
        MoldescriptorReaderError
            If the number of a WATER_TYPE or AMMONIA_TYPE line is not an
            integer.

        Notes
        -----
        A warning is logged if the water type does not point to a mol type
        of the file or the mol type is not made of the three atoms O, H and
        H in this order (as PQ requires), and if the ammonia type does not
        point to a mol type of the file.
        """
        self.water_type = None
        self.ammonia_type = None

        with open(self.filename, 'r', encoding='utf-8') as f:
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
                    self.water_type = self._read_type_number(splitted_line)
                    counter += 1
                elif splitted_line[0].lower() == "ammonia_type":
                    self.ammonia_type = self._read_type_number(splitted_line)
                    counter += 1
                else:
                    if len(splitted_line) != 3:
                        line = line.replace('\n', '')
                        self.logger.error(
                            (
                            "The number of columns in the header of a mol "
                            f"type must be 3.\nGot {len(splitted_line)} "
                            f"columns instead in text: '{line}'\n"
                            ),
                            exception=MoldescriptorReaderError
                        )

                    n_atoms = int(splitted_line[1])

                    if counter + n_atoms + 1 > len(lines):
                        self.logger.error(
                            (
                            f"The mol type {splitted_line[0]} declares "
                            f"{n_atoms} atoms, but only "
                            f"{len(lines) - counter - 1} atom lines follow "
                            "before the end of the file.\n"
                            ),
                            exception=MoldescriptorReaderError
                        )

                    mol_types.append(
                        self._read_mol_type(
                        lines[counter:counter + n_atoms + 1],
                        len(mol_types) + 1
                        )
                    )

                    counter += n_atoms + 1

        self._check_special_types(mol_types)

        return mol_types

    def _read_type_number(self, splitted_line: List[str]) -> int:
        """Parses the number of a WATER_TYPE or AMMONIA_TYPE line."""

        try:
            if len(splitted_line) != 2:
                raise ValueError
            return int(splitted_line[1])
        except ValueError:
            self.logger.error(
                (
                    f"The line of {splitted_line[0]} needs exactly one "
                    f"integer, got: '{' '.join(splitted_line)}'"
                ),
                exception=MoldescriptorReaderError
            )

        return None  # pragma: no cover - the error above raises

    def _check_special_types(self, mol_types: Residues) -> None:
        """Warns about water and ammonia types that do not fit the file."""

        numbers = {mol_type.id: mol_type for mol_type in mol_types}

        for key, number in (
            ("WATER_TYPE", self.water_type),
            ("AMMONIA_TYPE", self.ammonia_type),
        ):
            if number is not None and number not in numbers:
                self.logger.warning(
                    (
                        f"{key} {number} points to a mol type that is not in "
                        f"the file, which has the mol types 1 to "
                        f"{len(mol_types)}."
                    )
                )

        water = numbers.get(self.water_type)
        if water is not None and [e.symbol.lower() for e in water.elements] != [
            "o", "h", "h"
        ]:
            self.logger.warning(
                (
                    f"WATER_TYPE {self.water_type} points to the mol type "
                    f"'{water.name}', which is not made of the atoms O, H "
                    "and H in this order, as PQ requires."
                )
            )

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
                cls.logger.error(
                    (
                    "The number of columns in the body of a mol "
                    f"type must be 3 or 4.\nGot {len(splitted_line)} "
                    f"columns instead in text: '{line}'\n"
                    ),
                    exception=MoldescriptorReaderError
                )

            elements.append(Element(splitted_line[0]))
            atom_types.append(int(splitted_line[1]))
            partial_charges.append(float(splitted_line[2]))

        return Residue(
            name,
            mol_type_id,
            total_charge,
            elements,
            np.array(atom_types),
            np.array(partial_charges)
        )
