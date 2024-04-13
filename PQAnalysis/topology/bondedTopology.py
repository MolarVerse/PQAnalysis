from beartype.typing import List

from .bond import Bond
from .angle import Angle
from .dihedral import Dihedral


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
