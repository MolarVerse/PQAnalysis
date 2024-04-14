from beartype.typing import List

from PQAnalysis.io import BaseReader
from PQAnalysis.topology import Bond, BondedTopology, Angle, Dihedral


class TopologyFileReader(BaseReader):
    def __init__(self, filename: str) -> None:
        super().__init__(filename)

    def read(self) -> BondedTopology:
        blocks = self.get_definitions()
        return self.parse_blocks(blocks)

    def parse_blocks(self, blocks) -> BondedTopology:
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
                raise ValueError(f"Unknown block {key}")

        return BondedTopology(
            bonds=bonds,
            angles=angles,
            dihedrals=dihedrals,
            impropers=impropers,
            shake_bonds=shake_bonds
        )

    def get_definitions(self):
        with open(self.filename, "r") as file:
            lines = file.readlines()

            # remove all #.... parts from all lines
            lines = [line.split("#")[0] for line in lines]

            # remove all empty lines
            lines = [line for line in lines if line.strip()]

            # check if last line is END else raise error
            if lines[-1].strip() != "END":
                raise ValueError(
                    "Something went wrong. Each block should end with 'END'")

            # split all lines into blocks where each block ends with END
            blocks = []
            block = []
            for line in lines:
                if line.strip() == "END":
                    blocks.append(block)
                    block = []
                else:
                    block.append(line)

            # make a dictionary for each block with the key being the first word of the first line of each block
            # and the value being the rest of the block appart from the last line
            data = {}
            for block in blocks:
                key = block[0].split()[0].lower()
                value = block[1:]
                data[key] = value

            return data

    def parse_bonds(self, block) -> List[Bond]:
        bonds = []
        for line in block:
            if len(line.split()) == 4:
                index, target_index, bond_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 3:
                index, target_index, bond_type = line.split()
                is_linker = False
            else:
                # TODO: make an own exception for this
                raise ValueError(
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

    def parse_angles(self, block):
        angles = []
        for line in block:
            if len(line.split()) == 5:
                index1, index2, index3, angle_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 4:
                index1, index2, index3, angle_type = line.split()
                is_linker = False
            else:
                raise ValueError(
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

    def parse_dihedrals(self, block):
        dihedrals = []
        for line in block:
            if len(line.split()) == 6:
                index1, index2, index3, index4, dihedral_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 5:
                index1, index2, index3, index4, dihedral_type = line.split()
                is_linker = False
            else:
                raise ValueError(
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

    def parse_impropers(self, block):
        dihedrals = []
        for line in block:
            if len(line.split()) == 6:
                index1, index2, index3, index4, dihedral_type, _ = line.split()
                is_linker = True
            elif len(line.split()) == 5:
                index1, index2, index3, index4, dihedral_type = line.split()
                is_linker = False
            else:
                raise ValueError(
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

    def parse_shake(self, block) -> List[Bond]:
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
