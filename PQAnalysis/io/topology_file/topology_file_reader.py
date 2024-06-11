"""
A module to read a topology file for the bonded topology of the PQ or QMCFC project 
and return a BondedTopology object. For more information please visit the 
documentation page of PQ https://molarverse.github.io/PQ/
"""

import logging

from beartype.typing import List, Tuple

from PQAnalysis.io.base import BaseReader
from PQAnalysis.topology import Bond, BondedTopology, Angle, Dihedral
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.utils.string import is_comment_line
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis import __package_name__

from .exceptions import TopologyFileError



class TopologyFileReader(BaseReader):

    """
    A class to read a topology file for the bonded topology of the PQ or QMCFC 
    project and return a BondedTopology object. The topology file can have the 
    following blocks in any order and case-insensitive:
    - bonds
    - shake
    - angles
    - dihedrals
    - impropers

    For more information please visit the documentation page of PQ https://molarverse.github.io/PQ/
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the topology file.
        """
        super().__init__(filename)

    def read(self) -> BondedTopology:
        """
        Read the topology file and return a BondedTopology object.

        It gets all defined blocks in the topology file and parses them to 
        create a BondedTopology object.

        Returns
        -------
        BondedTopology
            A BondedTopology object.
        """
        blocks = self._get_definitions()
        return self._parse_blocks(blocks)

    def _get_definitions(self) -> dict[str, List[str]]:
        """
        Read the topology file and return a dictionary with the blocks.

        The blocks are defined by the first word of the first line of each block.
        The value of the dictionary is the rest of the block apart from the last line.

        Returns
        -------
        dict[str, List[str]]
            A dictionary with the blocks of the topology file. The key is the first
            word of the first line in lower case of each block and the value is 
            the rest of the block apart from the last line.

        Raises
        ------
        TopologyFileError
            If not all blocks end with 'END'.
        """
        with open(self.filename, "r", encoding="utf-8") as file:
            lines = file.readlines()

            # remove all comment lines and empty lines
            no_comment_lines = [
                line for line in lines if
                not is_comment_line(line, comment_token="#", empty_line=True)
            ]

            # check if last line is END else raise error
            if no_comment_lines[-1].strip().lower() != "end":
                self.logger.error(
                    "Something went wrong. Each block should end with 'END'",
                    exception=TopologyFileError,
                )

            # split all lines into blocks where each block ends with END
            blocks = []
            block = []
            for line in lines:
                if line.strip().lower() == "end":
                    block.append(line)
                    blocks.append(block)
                    block = []
                elif not is_comment_line(
                    line, comment_token="#", empty_line=True
                ):
                    block.append(line)

            # make a dictionary for each block with the key
            # being the first word of the first line of each block
            # and the value being the rest of the block apart from the last line
            data = {}
            for block in blocks:
                key = block[0].split()[0].lower()
                value = block[1:]
                data[key] = value

            return data

    def _parse_blocks(self, blocks: dict[str, List[str]]) -> BondedTopology:
        """
        Parse the blocks of the topology file and return a BondedTopology object.

        One block should have the following format:
        KEY #can also be in lower case
        ...
        ...
        END #can also be in lower case

        Parameters
        ----------
        blocks : dict[str, List[str]]
            A dictionary with the blocks of the topology file. The key is the first 
            word of the first line in lower case of each block and the value is the
            rest of the block apart from the last line.

        Returns
        -------
        BondedTopology
            A BondedTopology object.

        Raises
        ------
        TopologyFileError
            If an unknown block is found.
        """
        bonds = None
        shake_bonds = None
        angles = None
        dihedrals = None
        impropers = None

        for key, value in blocks.items():
            if key == "bonds":
                bonds = self._parse_bonds(value)
            elif key == "shake":
                shake_bonds = self._parse_shake(value)
            elif key == "angles":
                angles = self._parse_angles(value)
            elif key == "dihedrals":
                dihedrals = self._parse_dihedrals(value)
            elif key == "impropers":
                impropers = self._parse_impropers(value)
            else:
                self.logger.error(
                    f"Unknown block {key}", exception=TopologyFileError
                )

        return BondedTopology(
            bonds=bonds,
            angles=angles,
            dihedrals=dihedrals,
            impropers=impropers,
            shake_bonds=shake_bonds
        )

    def _parse_bonds(self, block: List[str]) -> List[Bond]:
        """
        Parse the bonds block of the topology file and return a list of Bond objects.

        One block should have the following format:

        bond_index target_index bond_type (*)
        ...
        ...
        bond_index target_index bond_type (*)

        Parameters
        ----------
        block : List[str]
            A list with the lines of the bonds block.

        Returns
        -------
        List[Bond]
            A list of Bond objects.

        Raises
        ------
        TopologyFileError
            If the number of columns in the block is not 3 or 4.
        """
        bonds = []
        for line in block[:-1]:  # [-1] to avoid the "END" line of the block

            line, comment = self._get_data_line_comment(line)

            index = None  # to avoid linter warning
            target_index = None  # to avoid linter warning
            bond_type = None  # to avoid linter warning
            is_linker = None  # to avoid linter warning

            if len(line.split()) == 4:
                index, target_index, bond_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 3:
                index, target_index, bond_type = line.split()
                is_linker = False
            else:
                self.logger.error(
                    "Invalid number of columns in bond block. Expected 3 or 4.",
                    exception=TopologyFileError
                )

            bonds.append(
                Bond(
                    index1=int(index),
                    index2=int(target_index),
                    bond_type=int(bond_type),
                    is_linker=is_linker,
                    comment=comment
                )
            )

        return bonds

    def _parse_angles(self, block: List[str]) -> List[Angle]:
        """
        Parse the angles block of the topology file and return a list of Angle objects.

        One block should have the following format:
        angle_index1 angle_index2 angle_index3 angle_type (*)
        ...
        ...
        angle_index1 angle_index2 angle_index3 angle_type (*)

        Parameters
        ----------
        block : List[str]
            A list with the lines of the angles block.

        Returns
        -------
        List[Angle]
            A list of Angle objects.

        Raises
        ------
        TopologyFileError
            If the number of columns in the block is not 4 or 5.
        """
        angles = []
        for line in block[:-1]:  # [-1] to avoid the "END" line of the block

            line, comment = self._get_data_line_comment(line)

            index1 = None  # to avoid linter warning
            index2 = None  # to avoid linter warning
            index3 = None  # to avoid linter warning
            angle_type = None  # to avoid linter warning
            is_linker = None  # to avoid linter warning

            if len(line.split()) == 5:
                index1, index2, index3, angle_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 4:
                index1, index2, index3, angle_type = line.split()
                is_linker = False
            else:
                self.logger.error(
                    "Invalid number of columns in angle block. Expected 4 or 5.",
                    exception=TopologyFileError
                )

            angles.append(
                Angle(
                    index1=int(index1),
                    index2=int(index2),
                    index3=int(index3),
                    angle_type=int(angle_type),
                    is_linker=is_linker,
                    comment=comment
                )
            )

        return angles

    def _parse_dihedrals(self, block: List[str]) -> List[Dihedral]:
        """
        Parse the dihedrals block of the topology file and return a list of Dihedral objects.

        One block should have the following format:
        dihedral_index1 dihedral_index2 dihedral_index3 dihedral_index4 dihedral_type (*)
        ...
        ...
        dihedral_index1 dihedral_index2 dihedral_index3 dihedral_index4 dihedral_type (*)

        Parameters
        ----------
        block : List[str]
            A list with the lines of the dihedrals block.

        Returns
        -------
        List[Dihedral]
            A list of Dihedral objects.

        Raises
        ------
        TopologyFileError
            If the number of columns in the block is not 5 or 6.
        """
        dihedrals = []
        for line in block[:-1]:  # [-1] to avoid the "END" line of the block

            line, comment = self._get_data_line_comment(line)

            index1 = None  # to avoid linter warning
            index2 = None  # to avoid linter warning
            index3 = None  # to avoid linter warning
            index4 = None  # to avoid linter warning
            dihedral_type = None  # to avoid linter warning
            is_linker = None  # to avoid linter warning

            if len(line.split()) == 6:
                index1, index2, index3, index4, dihedral_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 5:
                index1, index2, index3, index4, dihedral_type = line.split()
                is_linker = False
            else:
                self.logger.error(
                    "Invalid number of columns in dihedral block. Expected 5 or 6.",
                    exception=TopologyFileError
                )

            dihedrals.append(
                Dihedral(
                    index1=int(index1),
                    index2=int(index2),
                    index3=int(index3),
                    index4=int(index4),
                    dihedral_type=int(dihedral_type),
                    is_linker=is_linker,
                    comment=comment
                )
            )

        return dihedrals

    def _parse_impropers(self, block: List[str]) -> List[Dihedral]:
        """
        Parse the impropers block of the topology file and return a list of Dihedral objects.

        One block should have the following format:
        improper_index1 improper_index2 improper_index3 improper_index4 improper_type (*)
        ...
        ...
        improper_index1 improper_index2 improper_index3 improper_index4 improper_type (*)

        Parameters
        ----------
        block : List[str]
            A list with the lines of the impropers block.

        Returns
        -------
        List[Dihedral]
            A list of Dihedral objects.

        Raises
        ------
        TopologyFileError
            If the number of columns in the block is not 5 or 6.
        """
        dihedrals = []
        for line in block[:-1]:  # [-1] to avoid the "END" line of the block

            line, comment = self._get_data_line_comment(line)

            index1 = None  # to avoid linter warning
            index2 = None  # to avoid linter warning
            index3 = None  # to avoid linter warning
            index4 = None  # to avoid linter warning
            dihedral_type = None  # to avoid linter warning
            is_linker = None  # to avoid linter warning

            if len(line.split()) == 6:
                index1, index2, index3, index4, dihedral_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 5:
                index1, index2, index3, index4, dihedral_type = line.split()
                is_linker = False
            else:
                self.logger.error(
                    "Invalid number of columns in dihedral block. Expected 5 or 6.",
                    exception=TopologyFileError
                )

            dihedrals.append(
                Dihedral(
                    index1=int(index1),
                    index2=int(index2),
                    index3=int(index3),
                    index4=int(index4),
                    dihedral_type=int(dihedral_type),
                    is_linker=is_linker,
                    is_improper=True,
                    comment=comment
                )
            )

        return dihedrals

    def _parse_shake(self, block: List[str]) -> List[Bond]:
        """
        Parse the shake block of the topology file and return a list of Bond objects.

        One block should have the following format:
        bond_index target_index equilibrium-distance (*)
        ...
        ...
        bond_index target_index equilibrium-distance (*)

        Parameters
        ----------
        block : List[str]
            A list with the lines of the shake block.

        Returns
        -------
        List[Bond]
            A list of Bond objects.

        Raises
        ------
        TopologyFileError
            If the number of columns in the block is not 3 or 4.
        """
        shake_bonds = []
        for line in block[:-1]:  # [-1] to avoid the "END" line of the block

            line, comment = self._get_data_line_comment(line)

            index = None  # to avoid linter warning
            target_index = None  # to avoid linter warning
            distance = None  # to avoid linter warning
            is_linker = None  # to avoid linter warning

            if len(line.split()) == 4:
                index, target_index, distance, _ = line.split()
                is_linker = True
            else:
                index, target_index, distance = line.split()
                is_linker = False

            shake_bonds.append(
                Bond(
                    index1=int(index),
                    index2=int(target_index),
                    equilibrium_distance=float(distance),
                    is_linker=is_linker,
                    is_shake=True,
                    comment=comment
                )
            )

        return shake_bonds

    @staticmethod
    def _get_data_line_comment(line: str) -> Tuple[str, str | None]:
        """
        Get the data and the comment from a line.

        Parameters
        ----------
        line : str
            A line from the topology file.

        Returns
        -------
        tuple[str, str]
            A tuple with the data and the comment.
        """
        splitted_line = line.split("#")
        data = splitted_line[0].strip()

        if len(splitted_line) > 1:
            comment = "#".join(splitted_line[1:]).strip()
        else:
            comment = None

        return data, comment
