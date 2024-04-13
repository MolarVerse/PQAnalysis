from beartype.typing import List

from PQAnalysis.io import BaseWriter, FileWritingMode
from PQAnalysis.topology import Bond


class TopologyFileWriter(BaseWriter):
    def __init__(self,
                 filename: str,
                 shake_bonds: List[Bond] = None,
                 mode: str | FileWritingMode = "w"
                 ) -> None:

        super().__init__(filename)
        self.shake_bonds = shake_bonds if shake_bonds is not None else []

    def write(self) -> None:
        """
        Writes the shake bonds to the file.
        """

        self.open()

        n_unique_i, n_unique_j, n_linker = write_shake_header(self.shake_bonds)
        print(f"SHAKE {n_unique_i} {n_unique_j} {n_linker}", file=self.file)

        for bond in self.shake_bonds:
            self.file.write(
                f"{bond.index1} {bond.index2} {bond.equilibrium_distance} {"*" if bond.is_linker else ""}\n")

        self.close()


def write_shake_header(bonds: List[Bond]):  # TODO: Add type hint
    n_unique_indices = len(set([bond.index1 for bond in bonds]))
    n_unique_target_indices = len(set([bond.index2 for bond in bonds]))
    n_linker = len([bond for bond in bonds if bond.is_linker])

    return n_unique_indices, n_unique_target_indices, n_linker
