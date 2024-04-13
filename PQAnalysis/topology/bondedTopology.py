from beartype.typing import List

from .bond import Bond
from .angle import Angle
from .dihedral import Dihedral

from PQAnalysis.types import PositiveInt


class BondedTopology:
    def __init__(self,
                 bonds: List[Bond] | None = None,
                 angles: List[Angle] | None = None,
                 dihedrals: List[Dihedral] | None = None,
                 impropers: List[Dihedral] | None = None,
                 shake_bonds: List[Bond] | None = None,
                 ) -> None:

        self.bonds = bonds or []
        self.angles = angles or []
        self.dihedrals = dihedrals or []
        self.impropers = impropers or []
        self.shake_bonds = shake_bonds or []

    def extend_shake_bonds(self,
                           shake_bonds: List[Bond],
                           n_atoms: PositiveInt,
                           ) -> None:
        for bond in shake_bonds:
            bond.index1 += n_atoms
            bond.index2 += n_atoms
            self.shake_bonds.append(bond)

    @property
    def unique_shake_indices(self):
        return set([bond.index1 for bond in self.shake_bonds])

    @property
    def unique_shake_target_indices(self):
        return set([bond.index2 for bond in self.shake_bonds])

    @property
    def shake_linkers(self):
        return [bond for bond in self.shake_bonds if bond.is_linker]
