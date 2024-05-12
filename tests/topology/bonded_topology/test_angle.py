import numpy as np

from .. import pytestmark

from PQAnalysis.topology import Angle



class TestAngle:

    def test__init__(self):
        angle = Angle(index1=1, index2=2, index3=3)
        assert angle.index1 == 1
        assert angle.index2 == 2
        assert angle.index3 == 3
        assert angle.equilibrium_angle is None
        assert angle.angle_type is None
        assert not angle.is_linker

        angle = Angle(
            index1=1,
            index2=2,
            index3=3,
            equilibrium_angle=1.0,
            angle_type=1,
            is_linker=True,
        )
        assert angle.index1 == 1
        assert angle.index2 == 2
        assert angle.index3 == 3
        assert np.isclose(angle.equilibrium_angle, 1.0)
        assert angle.angle_type == 1
        assert angle.is_linker

    def test_copy(self):
        angle = Angle(index1=1, index2=2, index3=3)
        angle_copy = angle.copy()
        assert angle_copy.index1 == 1
        assert angle_copy.index2 == 2
        assert angle_copy.index3 == 3
        assert angle_copy.equilibrium_angle is None
        assert angle_copy.angle_type is None
        assert not angle_copy.is_linker

        angle = Angle(
            index1=1,
            index2=2,
            index3=3,
            equilibrium_angle=1.0,
            angle_type=1,
            is_linker=True,
        )
        angle_copy = angle.copy()
        assert angle_copy.index1 == 1
        assert angle_copy.index2 == 2
        assert angle_copy.index3 == 3
        assert np.isclose(angle_copy.equilibrium_angle, 1.0)
        assert angle_copy.angle_type == 1
        assert angle_copy.is_linker

        angle_copy.index1 = 2
        angle_copy.index2 = 3
        angle_copy.index3 = 4
        angle_copy.equilibrium_angle = 2.0
        angle_copy.angle_type = 2
        angle_copy.is_linker = False

        assert angle_copy.index1 != angle.index1
        assert angle_copy.index2 != angle.index2
        assert angle_copy.index3 != angle.index3
        assert angle_copy.equilibrium_angle != angle.equilibrium_angle
        assert angle_copy.angle_type != angle.angle_type
        assert angle_copy.is_linker != angle.is_linker

    def test__eq__(self):
        angle = Angle(index1=1, index2=2, index3=3)
        assert angle == angle

        angle_copy = angle.copy()
        assert angle == angle_copy

        angle_copy.index1 = 2
        assert angle != angle_copy

        assert angle != 1
        assert angle != "a"
        assert angle != None
