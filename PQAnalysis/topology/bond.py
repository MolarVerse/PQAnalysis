from PQAnalysis.types import PositiveInt, PositiveReal


class Bond:
    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 equilibrium_distance: PositiveReal,
                 is_linker: bool = False
                 ) -> None:

        self.index1 = index1
        self.index2 = index2
        self.equilibrium_distance = equilibrium_distance
        self.is_linker = is_linker
