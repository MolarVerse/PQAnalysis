import pytest
import numpy as np

from PQAnalysis.exceptions import PQValueError

from PQAnalysis.topology import (
    BondedTopology,
    Bond,
    Angle,
    Dihedral,
)

from .. import pytestmark



class TestBondedTopology:

    def test__init__(self):
        bonded_topology = BondedTopology()
        assert bonded_topology.bonds == []
        assert bonded_topology.angles == []
        assert bonded_topology.dihedrals == []
        assert bonded_topology.impropers == []
        assert bonded_topology.shake_bonds == []
        assert bonded_topology.ordering_keys is None

        bond = Bond(index1=1, index2=2)
        angle = Angle(index1=1, index2=2, index3=3)
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)

        bonded_topology = BondedTopology(
            bonds=[bond],
            angles=[angle],
            dihedrals=[dihedral],
            impropers=[dihedral],
            shake_bonds=[bond],
            ordering_keys=[
            "bonds",
            "angles",
            "dihedrals",
            "impropers",
            "shake_bonds"
            ]
        )
        assert bonded_topology.bonds == [bond]
        assert bonded_topology.angles == [angle]
        assert bonded_topology.dihedrals == [dihedral]
        assert bonded_topology.impropers == [dihedral]
        assert bonded_topology.shake_bonds == [bond]
        assert bonded_topology.ordering_keys == [
            "bonds",
            "angles",
            "dihedrals",
            "impropers",
            "shake_bonds"
        ]

    def test_extend_shake_bonds(self):
        bond = Bond(index1=1, index2=2)
        bonded_topology = BondedTopology(bonds=[bond])

        shake_bond = Bond(index1=2, index2=3)
        bonded_topology.extend_shake_bonds(shake_bonds=[shake_bond], n_atoms=3)

        assert all(
            [
            a == b for a,
            b in zip(bonded_topology.shake_bonds,
            [Bond(index1=5,
            index2=6)])
            ]
        )

        bonded_topology = BondedTopology(bonds=[bond])

        with pytest.raises(PQValueError) as exception:
            bonded_topology.extend_shake_bonds(
                shake_bonds=[shake_bond],
                n_atoms=3,
                n_extensions=2
            )
        assert str(
            exception.value
        ) == "n_atoms_per_extension must be provided if n_extensions is not 1."

        with pytest.raises(PQValueError) as exception:
            bonded_topology.extend_shake_bonds(
                shake_bonds=[shake_bond],
                n_atoms=3,
                n_extensions=2,
                n_atoms_per_extension=1
            )
        assert str(
            exception.value
        ) == "n_atoms_per_extension must be greater or equal than the highest index in the provided shake bonds."

        bonded_topology.extend_shake_bonds(
            shake_bonds=[shake_bond],
            n_atoms=3,
            n_extensions=2,
            n_atoms_per_extension=3
        )

        assert all(
            [
            a == b for a,
            b in zip(
            bonded_topology.shake_bonds,
            [Bond(index1=5,
            index2=6),
            Bond(index1=8,
            index2=9)]
            )
            ]
        )

    def test_unique_bond1_indices(self):
        bond = Bond(index1=1, index2=2)
        bonded_topology = BondedTopology(bonds=[bond, bond])

        assert bonded_topology.unique_bond1_indices == {1}

    def test_unique_bond2_indices(self):
        bond = Bond(index1=1, index2=2)
        bonded_topology = BondedTopology(bonds=[bond, bond])

        assert bonded_topology.unique_bond2_indices == {2}

    def test_bond_linkers(self):
        bond = Bond(index1=1, index2=2, is_linker=True)
        bonded_topology = BondedTopology(bonds=[bond, bond])

        assert bonded_topology.bond_linkers == [bond, bond]

    def test_unique_angle1_indices(self):
        angle = Angle(index1=1, index2=2, index3=3)
        bonded_topology = BondedTopology(angles=[angle, angle])

        assert bonded_topology.unique_angle1_indices == {1}

    def test_unique_angle2_indices(self):
        angle = Angle(index1=1, index2=2, index3=3)
        bonded_topology = BondedTopology(angles=[angle, angle])

        assert bonded_topology.unique_angle2_indices == {2}

    def test_unique_angle3_indices(self):
        angle = Angle(index1=1, index2=2, index3=3)
        bonded_topology = BondedTopology(angles=[angle, angle])

        assert bonded_topology.unique_angle3_indices == {3}

    def test_angle_linkers(self):
        angle = Angle(index1=1, index2=2, index3=3, is_linker=True)
        bonded_topology = BondedTopology(angles=[angle, angle])

        assert bonded_topology.angle_linkers == [angle, angle]

    def test_unique_dihedral1_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(dihedrals=[dihedral, dihedral])

        assert bonded_topology.unique_dihedral1_indices == {1}

    def test_unique_dihedral2_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(dihedrals=[dihedral, dihedral])

        assert bonded_topology.unique_dihedral2_indices == {2}

    def test_unique_dihedral3_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(dihedrals=[dihedral, dihedral])

        assert bonded_topology.unique_dihedral3_indices == {3}

    def test_unique_dihedral4_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(dihedrals=[dihedral, dihedral])

        assert bonded_topology.unique_dihedral4_indices == {4}

    def test_dihedral_linkers(self):
        dihedral = Dihedral(
            index1=1,
            index2=2,
            index3=3,
            index4=4,
            is_linker=True
        )
        bonded_topology = BondedTopology(dihedrals=[dihedral, dihedral])

        assert bonded_topology.dihedral_linkers == [dihedral, dihedral]

    def test_unique_improper1_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(impropers=[dihedral, dihedral])

        assert bonded_topology.unique_improper1_indices == {1}

    def test_unique_improper2_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(impropers=[dihedral, dihedral])

        assert bonded_topology.unique_improper2_indices == {2}

    def test_unique_improper3_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(impropers=[dihedral, dihedral])

        assert bonded_topology.unique_improper3_indices == {3}

    def test_unique_improper4_indices(self):
        dihedral = Dihedral(index1=1, index2=2, index3=3, index4=4)
        bonded_topology = BondedTopology(impropers=[dihedral, dihedral])

        assert bonded_topology.unique_improper4_indices == {4}

    def test_improper_linkers(self):
        dihedral = Dihedral(
            index1=1,
            index2=2,
            index3=3,
            index4=4,
            is_linker=True
        )
        bonded_topology = BondedTopology(impropers=[dihedral, dihedral])

        assert bonded_topology.improper_linkers == [dihedral, dihedral]

    def test_unique_shake_indices(self):
        bond = Bond(index1=1, index2=2)
        bonded_topology = BondedTopology(shake_bonds=[bond, bond])

        assert bonded_topology.unique_shake_indices == {1}

    def test_unique_shake_target_indices(self):
        bond = Bond(index1=1, index2=2)
        bonded_topology = BondedTopology(shake_bonds=[bond, bond])

        assert bonded_topology.unique_shake_target_indices == {2}

    def test_shake_linkers(self):
        bond = Bond(index1=1, index2=2, is_linker=True)
        bonded_topology = BondedTopology(shake_bonds=[bond, bond])

        assert bonded_topology.shake_linkers == [bond, bond]
