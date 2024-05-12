import numpy as np

from .. import pytestmark

from PQAnalysis.topology import Dihedral



class TestDihedral:

    def test__init__(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        assert dihedral.index1 == 1
        assert dihedral.index2 == 2
        assert dihedral.index3 == 3
        assert dihedral.index4 == 4
        assert dihedral.equilibrium_angle is None
        assert dihedral.dihedral_type is None
        assert not dihedral.is_linker
        assert not dihedral.is_improper

        dihedral = Dihedral(
            index1=1,
            index2=2,
            index3=3,
            index4=4,
            equilibrium_angle=1.0,
            dihedral_type=1,
            is_linker=True,
            is_improper=True
        )
        assert dihedral.index1 == 1
        assert dihedral.index2 == 2
        assert dihedral.index3 == 3
        assert dihedral.index4 == 4
        assert np.isclose(dihedral.equilibrium_angle, 1.0)
        assert dihedral.dihedral_type == 1
        assert dihedral.is_linker
        assert dihedral.is_improper

    def test_copy(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        dihedral_copy = dihedral.copy()
        assert dihedral_copy.index1 == 1
        assert dihedral_copy.index2 == 2
        assert dihedral_copy.index3 == 3
        assert dihedral_copy.index4 == 4
        assert dihedral_copy.equilibrium_angle is None
        assert dihedral_copy.dihedral_type is None
        assert not dihedral_copy.is_linker
        assert not dihedral_copy.is_improper

        dihedral = Dihedral(
            index1=1,
            index2=2,
            index3=3,
            index4=4,
            equilibrium_angle=1.0,
            dihedral_type=1,
            is_linker=True,
            is_improper=True
        )
        dihedral_copy = dihedral.copy()
        assert dihedral_copy.index1 == 1
        assert dihedral_copy.index2 == 2
        assert dihedral_copy.index3 == 3
        assert dihedral_copy.index4 == 4
        assert np.isclose(dihedral_copy.equilibrium_angle, 1.0)
        assert dihedral_copy.dihedral_type == 1
        assert dihedral_copy.is_linker
        assert dihedral_copy.is_improper

        dihedral_copy.index1 = 2
        dihedral_copy.index2 = 3
        dihedral_copy.index3 = 4
        dihedral_copy.index4 = 5
        dihedral_copy.equilibrium_angle = 2.0
        dihedral_copy.dihedral_type = 2
        dihedral_copy.is_linker = False
        dihedral_copy.is_improper = False

        assert dihedral_copy.index1 != dihedral.index1
        assert dihedral_copy.index2 != dihedral.index2
        assert dihedral_copy.index3 != dihedral.index3
        assert dihedral_copy.index4 != dihedral.index4
        assert dihedral_copy.equilibrium_angle != dihedral.equilibrium_angle
        assert dihedral_copy.dihedral_type != dihedral.dihedral_type
        assert dihedral_copy.is_linker != dihedral.is_linker
        assert dihedral_copy.is_improper != dihedral.is_improper

    def test__eq__(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        assert dihedral == dihedral

        dihedral_copy = dihedral.copy()
        assert dihedral == dihedral_copy

        dihedral_copy.index1 = 2
        assert dihedral != dihedral_copy

        assert dihedral != 1
        assert dihedral != "a"
        assert dihedral != None
