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
                           n_extensions: PositiveInt = 1,
                           n_atoms_per_extension: PositiveInt | None = None
                           ) -> None:

        if n_extensions != 1 and n_atoms_per_extension is None:
            raise ValueError(
                "n_atoms_per_extension must be provided if n_extensions is not 1."
            )

        max_index_shake_bonds = max(
            [bond.index1 for bond in shake_bonds] +
            [bond.index2 for bond in shake_bonds]
        )

        if n_atoms_per_extension is None:
            n_atoms_per_extension = 0
        elif n_atoms_per_extension < max_index_shake_bonds:
            raise ValueError(
                "n_atoms_per_extension must be greater or equal than the highest index in the provided shake bonds."
            )

        for i in range(n_extensions):

            _n_atoms = n_atoms_per_extension * i + n_atoms

            for bond in shake_bonds:
                _bond = bond.copy()
                _bond.index1 += _n_atoms
                _bond.index2 += _n_atoms
                self.shake_bonds.append(_bond)

    @property
    def unique_shake_indices(self):
        return set([bond.index1 for bond in self.shake_bonds])

    @property
    def unique_shake_target_indices(self):
        return set([bond.index2 for bond in self.shake_bonds])

    @property
    def shake_linkers(self):
        return [bond for bond in self.shake_bonds if bond.is_linker]
