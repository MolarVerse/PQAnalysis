from PQAnalysis.types import PositiveInt, PositiveReal


class Angle:
    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 index3: PositiveInt,
                 equilibrium_angle: PositiveReal | None = None,
                 angle_type: PositiveInt | None = None,
                 is_linker: bool = False,
                 ) -> None:

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.equilibrium_angle = equilibrium_angle
        self.angle_type = angle_type
        self.is_linker = is_linker
