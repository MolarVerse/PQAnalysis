from beartype.typing import List

from PQAnalysis.io import BaseWriter, FileWritingMode
from PQAnalysis.topology import Bond, BondedTopology, Topology


class TopologyFileWriter(BaseWriter):
    def __init__(self,
                 filename: str,
                 mode: str | FileWritingMode = "w"
                 ) -> None:

        super().__init__(filename, mode=mode)

    def write(self, bonded_topology: Topology | BondedTopology | None) -> None:
        """
        Writes the bonded topology to the file.
        """

        if isinstance(bonded_topology, Topology) and bonded_topology.bonded_topology is not None:
            bonded_topology = bonded_topology.bonded_topology
        elif isinstance(bonded_topology, BondedTopology):
            raise ValueError("No bonded topology provided in Topology object.")
        elif bonded_topology is None:
            raise ValueError("No bonded topology provided.")

        self.open()

        if len(bonded_topology.shake_bonds) != 0:
            lines = get_shake_lines(bonded_topology)
            for line in lines:
                self.file.write(line)

        self.close()


def get_bond_lines(bonded_topology: BondedTopology) -> List[Bond]:
    lines = []

    for bond in bonded_topology.bonds:
        lines.append(
            f"{bond.index1} {bond.index2} {bond.equilibrium_distance}"
        )

    return lines


def get_shake_lines(bonded_topology: BondedTopology) -> List[Bond]:
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
