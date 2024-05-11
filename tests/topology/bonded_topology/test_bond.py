import numpy as np

from .. import pytestmark

from PQAnalysis.topology import Bond



class TestBond:

    def test__init__(self):
        bond = Bond(index1=1, index2=2)
        assert bond.index1 == 1
        assert bond.index2 == 2
        assert bond.equilibrium_distance is None
        assert bond.bond_type is None
        assert not bond.is_linker
        assert not bond.is_shake

        bond = Bond(
            index1=1,
            index2=2,
            equilibrium_distance=1.0,
            bond_type=1,
            is_linker=True,
            is_shake=True
        )
        assert bond.index1 == 1
        assert bond.index2 == 2
        assert np.isclose(bond.equilibrium_distance, 1.0)
        assert bond.bond_type == 1
        assert bond.is_linker
        assert bond.is_shake

    def test_copy(self):
        bond = Bond(index1=1, index2=2)
        bond_copy = bond.copy()
        assert bond_copy.index1 == 1
        assert bond_copy.index2 == 2
        assert bond_copy.equilibrium_distance is None
        assert bond_copy.bond_type is None
        assert not bond_copy.is_linker
        assert not bond_copy.is_shake

        bond = Bond(
            index1=1,
            index2=2,
            equilibrium_distance=1.0,
            bond_type=1,
            is_linker=True,
            is_shake=True
        )
        bond_copy = bond.copy()
        assert bond_copy.index1 == 1
        assert bond_copy.index2 == 2
        assert np.isclose(bond_copy.equilibrium_distance, 1.0)
        assert bond_copy.bond_type == 1
        assert bond_copy.is_linker
        assert bond_copy.is_shake

        bond_copy.index1 = 2
        bond_copy.index2 = 3
        bond_copy.equilibrium_distance = 2.0
        bond_copy.bond_type = 2
        bond_copy.is_linker = False
        bond_copy.is_shake = False

        assert bond.index1 != bond_copy.index1
        assert bond.index2 != bond_copy.index2
        assert bond.equilibrium_distance != bond_copy.equilibrium_distance
        assert bond.bond_type != bond_copy.bond_type
        assert bond.is_linker != bond_copy.is_linker
        assert bond.is_shake != bond_copy.is_shake

    def test__eq__(self):
        bond = Bond(index1=1, index2=2)
        bond_copy = bond.copy()
        assert bond == bond_copy

        bond = Bond(
            index1=1,
            index2=2,
            equilibrium_distance=1.0,
            bond_type=1,
            is_linker=True,
            is_shake=True
        )
        bond_copy = bond.copy()
        assert bond == bond_copy

        bond_copy.index1 = 2
        bond_copy.index2 = 3
        bond_copy.equilibrium_distance = 2.0
        bond_copy.bond_type = 2
        bond_copy.is_linker = False
        bond_copy.is_shake = False

        assert bond != bond_copy
        assert bond_copy != 1
        assert bond != 1
        assert bond != "1"
        assert bond != None
