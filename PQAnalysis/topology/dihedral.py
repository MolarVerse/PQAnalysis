from PQAnalysis.types import PositiveInt, PositiveReal


class Dihedral:
    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 index3: PositiveInt,
                 index4: PositiveInt,
                 equilibrium_angle: PositiveReal | None = None,
                 dihedral_type: int | None = None,
                 is_linker: bool = False,
                 is_improper: bool = False,
                 ) -> None:

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.index4 = index4
        self.equilibrium_angle = equilibrium_angle
        self.angle_type = dihedral_type
        self.is_linker = is_linker
        self.is_improper = is_improper
