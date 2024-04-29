"""
A module to read a topology file for the bonded topology of the PQ or QMCFC project 
and return a BondedTopology object. For more information please visit the 
documentation page of PQ https://molarverse.github.io/PQ/
"""

from beartype.typing import List

from PQAnalysis.io.base import BaseReader
from PQAnalysis.topology import Bond, BondedTopology, Angle, Dihedral
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
        blocks = self.get_definitions()
        return self.parse_blocks(blocks)

    def get_definitions(self) -> dict[str, List[str]]:
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

            # remove all #.... parts from all lines
            lines = [line.split("#")[0] for line in lines]

            # remove all empty lines
            lines = [line for line in lines if line.strip()]

            # check if last line is END else raise error
            if lines[-1].strip() != "END":
                raise TopologyFileError(
                    "Something went wrong. Each block should end with 'END'"
                )

            # split all lines into blocks where each block ends with END
            blocks = []
            block = []
            for line in lines:
                if line.strip().lower() == "END":
                    blocks.append(block)
                    block = []
                else:
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

    def parse_blocks(self, blocks: dict[str, List[str]]) -> BondedTopology:
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
                bonds = self.parse_bonds(value)
            elif key == "shake":
                shake_bonds = self.parse_shake(value)
            elif key == "angles":
                angles = self.parse_angles(value)
            elif key == "dihedrals":
                dihedrals = self.parse_dihedrals(value)
            elif key == "impropers":
                impropers = self.parse_impropers(value)
            else:
                raise TopologyFileError(f"Unknown block {key}")

        return BondedTopology(
            bonds=bonds,
            angles=angles,
            dihedrals=dihedrals,
            impropers=impropers,
            shake_bonds=shake_bonds
        )

    def parse_bonds(self, block: List[str]) -> List[Bond]:
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
        for line in block:
            if len(line.split()) == 4:
                index, target_index, bond_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 3:
                index, target_index, bond_type = line.split()
                is_linker = False
            else:
                raise TopologyFileError(
                    "Invalid number of columns in bond block. Expected 3 or 4."
                )

            bonds.append(
                Bond(
                    index1=int(index),
                    index2=int(target_index),
                    bond_type=int(bond_type),
                    is_linker=is_linker,
                )
            )

        return bonds

    def parse_angles(self, block: List[str]) -> List[Angle]:
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
        for line in block:
            if len(line.split()) == 5:
                index1, index2, index3, angle_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 4:
                index1, index2, index3, angle_type = line.split()
                is_linker = False
            else:
                raise TopologyFileError(
                    "Invalid number of columns in angle block. Expected 4 or 5."
                )

            angles.append(
                Angle(
                    index1=int(index1),
                    index2=int(index2),
                    index3=int(index3),
                    angle_type=int(angle_type),
                    is_linker=is_linker,
                )
            )

        return angles

    def parse_dihedrals(self, block: List[str]) -> List[Dihedral]:
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
        for line in block:
            if len(line.split()) == 6:
                index1, index2, index3, index4, dihedral_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 5:
                index1, index2, index3, index4, dihedral_type = line.split()
                is_linker = False
            else:
                raise TopologyFileError(
                    "Invalid number of columns in dihedral block. Expected 5 or 6."
                )

            dihedrals.append(
                Dihedral(
                    index1=int(index1),
                    index2=int(index2),
                    index3=int(index3),
                    index4=int(index4),
                    dihedral_type=int(dihedral_type),
                    is_linker=is_linker,
                )
            )

        return dihedrals

    def parse_impropers(self, block: List[str]) -> List[Dihedral]:
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
        for line in block:
            if len(line.split()) == 6:
                index1, index2, index3, index4, dihedral_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 5:
                index1, index2, index3, index4, dihedral_type = line.split()
                is_linker = False
            else:
                raise TopologyFileError(
                    "Invalid number of columns in dihedral block. Expected 5 or 6."
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
                )
            )

        return dihedrals

    def parse_shake(self, block: List[str]) -> List[Bond]:
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
        for line in block:
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
                )
            )

        return shake_bonds
