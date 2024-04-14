from beartype.typing import List

from PQAnalysis.io import BaseWriter, FileWritingMode
from PQAnalysis.topology import Bond, BondedTopology, Topology


class TopologyFileWriter(BaseWriter):
    def __init__(self,
                 filename: str,
                 mode: str | FileWritingMode = "w"
                 ) -> None:

        super().__init__(filename, mode=mode)

    def write(self, bonded_topology: Topology | BondedTopology) -> None:
        """
        Writes the bonded topology to the file.
        """

        if isinstance(bonded_topology, Topology) and bonded_topology.bonded_topology is not None:
            bonded_topology = bonded_topology.bonded_topology
        elif isinstance(bonded_topology, BondedTopology):
            bonded_topology = bonded_topology
        else:
            raise ValueError("Invalid bonded topology.")

        self.open()

        # TODO: may result in different section order

        if len(bonded_topology.bonds) != 0:
            lines = get_bond_lines(bonded_topology)
            for line in lines:
                print(line, file=self.file)

        if len(bonded_topology.angles) != 0:
            lines = get_angle_lines(bonded_topology)
            for line in lines:
                print(line, file=self.file)

        if len(bonded_topology.dihedrals) != 0:
            lines = get_dihedral_lines(bonded_topology)
            for line in lines:
                print(line, file=self.file)

        if len(bonded_topology.impropers) != 0:
            lines = get_improper_lines(bonded_topology)
            for line in lines:
                print(line, file=self.file)

        if len(bonded_topology.shake_bonds) != 0:
            lines = get_shake_lines(bonded_topology)
            for line in lines:
                print(line, file=self.file)

        self.close()


def get_bond_lines(bonded_topology: BondedTopology) -> List[str]:
    n_unique_indices = len(bonded_topology.unique_bond1_indices)
    n_unique_target_indices = len(bonded_topology.unique_bond2_indices)
    n_linkers = len(bonded_topology.linkers)

    lines = []

    lines.append(
        f"BONDS {n_unique_indices} {n_unique_target_indices} {n_linkers}"
    )

    for bond in bonded_topology.bonds:
        lines.append(
            f"{bond.index1} {bond.index2} {bond.bond_type}"
        )

    lines.append("END")

    return lines


def get_angle_lines(bonded_topology: BondedTopology) -> List[str]:
    n_unique_indices1 = len(bonded_topology.unique_angle1_indices)
    n_unique_indices2 = len(bonded_topology.unique_angle2_indices)
    n_unique_indices3 = len(bonded_topology.unique_angle_target_indices)
    n_linkers = len(bonded_topology.angle_linkers)

    lines = []

    lines.append(
        f"ANGLES {n_unique_indices1} {n_unique_indices2} {
            n_unique_indices3} {n_linkers}"
    )

    for angle in bonded_topology.angles:
        lines.append(
            f"{angle.index1} {angle.index2} {
                angle.index3} {angle.angle_type}"
        )

    lines.append("END")

    return lines


def get_dihedral_lines(bonded_topology: BondedTopology) -> List[str]:
    n_unique_indices1 = len(bonded_topology.unique_dihedral1_indices)
    n_unique_indices2 = len(bonded_topology.unique_dihedral2_indices)
    n_unique_indices3 = len(bonded_topology.unique_dihedral3_indices)
    n_unique_indices4 = len(bonded_topology.unique_dihedral4_indices)

    lines = []

    lines.append(
        f"DIHEDRALS {n_unique_indices1} {n_unique_indices2} {
            n_unique_indices3} {n_unique_indices4}"
    )

    for dihedral in bonded_topology.dihedrals:
        lines.append(
            f"{dihedral.index1} {dihedral.index2} {
                dihedral.index3} {dihedral.index4} {dihedral.dihedral_type}"
        )

    lines.append("END")

    return lines


def get_improper_lines(bonded_topology: BondedTopology) -> List[str]:
    n_unique_indices1 = len(bonded_topology.unique_improper1_indices)
    n_unique_indices2 = len(bonded_topology.unique_improper2_indices)
    n_unique_indices3 = len(bonded_topology.unique_improper3_indices)
    n_unique_indices4 = len(bonded_topology.unique_improper4_indices)

    lines = []

    lines.append(
        f"IMPROPERS {n_unique_indices1} {n_unique_indices2} {
            n_unique_indices3} {n_unique_indices4}"
    )

    for improper in bonded_topology.impropers:
        lines.append(
            f"{improper.index1} {improper.index2} {
                improper.index3} {improper.index4} {improper.improper_type}"
        )

    lines.append("END")

    return lines


def get_shake_lines(bonded_topology: BondedTopology) -> List[str]:
    n_unique_indices = len(bonded_topology.unique_shake_indices)
    n_unique_target_indices = len(bonded_topology.unique_shake_target_indices)
    n_linkers = len(bonded_topology.shake_linkers)

    lines = []

    lines.append(
        f"SHAKE {n_unique_indices} {n_unique_target_indices} {n_linkers}"
    )

    for bond in bonded_topology.shake_bonds:
        lines.append(
            f"{bond.index1} {bond.index2} {bond.equilibrium_distance} {
                "*" if bond.is_linker else ""}"
        )

    lines.append("END")

    return lines
