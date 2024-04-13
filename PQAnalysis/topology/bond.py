from PQAnalysis.types import PositiveInt, PositiveReal


class Bond:
    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 equilibrium_distance: PositiveReal | None = None,
                 bond_index: PositiveInt | None = None,
                 is_linker: bool = False,
                 is_shake: bool = False,
                 ) -> None:

        self.index1 = index1
        self.index2 = index2
        self.equilibrium_distance = equilibrium_distance
        self.bond_index = bond_index
        self.is_linker = is_linker
        self.is_shake = is_shake
