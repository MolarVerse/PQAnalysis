"""
Unit tests for the topologyWriter module.
"""

import pytest

from PQAnalysis.io.topology_file import TopologyFileWriter, TopologyFileError
from PQAnalysis.topology import BondedTopology
from PQAnalysis.topology.bonded_topology import (
    Bond,
    Angle,
    Dihedral,
)

from . import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access



class TestTopologyFileWriter:

    """
    Test cases for the TopologyFileWriter class.
    """

    def test__check_type_given(self):
        """
        Test the _check_type_given method.
        """
        bonds = [Bond(index1=1, index2=2, bond_type=1)]

        TopologyFileWriter._check_type_given(bonds, "bond")

        bonds.append(Bond(index1=2, index2=3))

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileWriter._check_type_given(bonds, "bond")
        assert str(exception.value) == (
            "In order to write the bond information in 'PQ' topology format, "
            "all bonds must have a bond type defined."
        )

        angles = [Angle(index1=1, index2=2, index3=3, angle_type=1)]

        TopologyFileWriter._check_type_given(angles, "angle")

        angles.append(Angle(index1=2, index2=3, index3=4))

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileWriter._check_type_given(angles, "angle")
        assert str(exception.value) == (
            "In order to write the angle information in 'PQ' topology format, "
            "all angles must have a angle type defined."
        )

        dihedrals = [
            Dihedral(index1=1, index2=2, index3=3, index4=4, dihedral_type=1)
        ]

        TopologyFileWriter._check_type_given(dihedrals, "dihedral")

        dihedrals.append(Dihedral(index1=2, index2=3, index3=4, index4=5))

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileWriter._check_type_given(dihedrals, "dihedral")
        assert str(exception.value) == (
            "In order to write the dihedral information in 'PQ' topology format, "
            "all dihedrals must have a dihedral type defined."
        )

    def test__get_bond_lines(self):
        """
        Test the _get_bond_lines method.
        """
        topology = BondedTopology()

        lines = TopologyFileWriter._get_bond_lines(topology)
        assert lines[0] == "BONDS 0 0 0"
        assert lines[1] == "END"

        bond1 = Bond(index1=1, index2=2, bond_type=1)
        bond2 = Bond(index1=2, index2=3, bond_type=1, is_linker=True)
        bond3 = Bond(
            index1=5,
            index2=2,
            bond_type=1,
            comment="This is a comment.",
        )
        bond4 = Bond(
            index1=5,
            index2=2,
            bond_type=1,
            is_linker=True,
            comment="This is a comment."
        )

        topology.bonds = [bond1, bond2, bond3, bond4]

        lines = TopologyFileWriter._get_bond_lines(topology)
        assert lines[0] == "BONDS 3 2 2"
        assert lines[1] == f"{1:>5d} {2:>5d} {1:>5d}"
        assert lines[2] == f"{2:>5d} {3:>5d} {1:>5d} *"
        assert lines[3] == f"{5:>5d} {2:>5d} {1:>5d} # This is a comment."
        assert lines[4] == f"{5:>5d} {2:>5d} {1:>5d} * # This is a comment."
        assert lines[5] == "END"

    def test__get_angle_lines(self):
        """
        Test the _get_angle_lines method.
        """
        topology = BondedTopology()

        lines = TopologyFileWriter._get_angle_lines(topology)
        assert lines[0] == "ANGLES 0 0 0 0"
        assert lines[1] == "END"

        angle1 = Angle(index1=1, index2=2, index3=3, angle_type=1)
        angle2 = Angle(
            index1=2,
            index2=3,
            index3=4,
            angle_type=1,
            is_linker=True,
        )
        angle3 = Angle(
            index1=5,
            index2=2,
            index3=3,
            angle_type=1,
            comment="This is a comment.",
        )
        angle4 = Angle(
            index1=5,
            index2=2,
            index3=3,
            angle_type=1,
            is_linker=True,
            comment="This is a comment."
        )

        topology.angles = [angle1, angle2, angle3, angle4]

        lines = TopologyFileWriter._get_angle_lines(topology)
        assert lines[0] == "ANGLES 3 2 2 2"
        assert lines[1] == f"{1:>5d} {2:>5d} {3:>5d} {1:>5d}"
        assert lines[2] == f"{2:>5d} {3:>5d} {4:>5d} {1:>5d} *"
        assert lines[3] == (
            f"{5:>5d} {2:>5d} {3:>5d} {1:>5d} "
            "# This is a comment."
        )
        assert lines[4] == (
            f"{5:>5d} {2:>5d} {3:>5d} {1:>5d} * # "
            "This is a comment."
        )
        assert lines[5] == "END"

    def test__get_dihedral_lines(self):
        """
        Test the _get_dihedral_lines method.
        """
        topology = BondedTopology()

        lines = TopologyFileWriter._get_dihedral_lines(topology)
        assert lines[0] == "DIHEDRALS 0 0 0 0"
        assert lines[1] == "END"

        dihedral1 = Dihedral(
            index1=1, index2=2, index3=3, index4=4, dihedral_type=1
        )
        dihedral2 = Dihedral(
            index1=2,
            index2=3,
            index3=4,
            index4=5,
            dihedral_type=1,
            is_linker=True,
        )
        dihedral3 = Dihedral(
            index1=5,
            index2=2,
            index3=3,
            index4=4,
            dihedral_type=1,
            comment="This is a comment.",
        )
        dihedral4 = Dihedral(
            index1=5,
            index2=2,
            index3=3,
            index4=4,
            dihedral_type=1,
            is_linker=True,
            comment="This is a comment."
        )

        topology.dihedrals = [dihedral1, dihedral2, dihedral3, dihedral4]

        lines = TopologyFileWriter._get_dihedral_lines(topology)
        assert lines[0] == "DIHEDRALS 3 2 2 2"
        assert lines[1] == f"{1:>5d} {2:>5d} {3:>5d} {4:>5d} {1:>5d}"
        assert lines[2] == f"{2:>5d} {3:>5d} {4:>5d} {5:>5d} {1:>5d} *"
        assert lines[3] == (
            f"{5:>5d} {2:>5d} {3:>5d} {4:>5d} {1:>5d} "
            "# This is a comment."
        )
        assert lines[4] == (
            f"{5:>5d} {2:>5d} {3:>5d} {4:>5d} {1:>5d} * # "
            "This is a comment."
        )
        assert lines[5] == "END"

    def test__get_improper_lines(self):
        """
        Test the _get_improper_lines method.
        """
        topology = BondedTopology()

        lines = TopologyFileWriter._get_improper_lines(topology)
        assert lines[0] == "IMPROPERS 0 0 0 0"
        assert lines[1] == "END"

        improper1 = Dihedral(
            index1=1,
            index2=2,
            index3=3,
            index4=4,
            dihedral_type=1,
            is_improper=True
        )
        improper2 = Dihedral(
            index1=2,
            index2=3,
            index3=4,
            index4=5,
            dihedral_type=1,
            is_linker=True,
            is_improper=True,
        )
        improper3 = Dihedral(
            index1=5,
            index2=2,
            index3=3,
            index4=4,
            dihedral_type=1,
            is_improper=True,
            comment="This is a comment.",
        )
        improper4 = Dihedral(
            index1=5,
            index2=2,
            index3=3,
            index4=4,
            dihedral_type=1,
            is_linker=True,
            is_improper=True,
            comment="This is a comment."
        )

        topology.impropers = [improper1, improper2, improper3, improper4]

        lines = TopologyFileWriter._get_improper_lines(topology)
        assert lines[0] == "IMPROPERS 3 2 2 2"
        assert lines[1] == f"{1:>5d} {2:>5d} {3:>5d} {4:>5d} {1:>5d}"
        assert lines[2] == f"{2:>5d} {3:>5d} {4:>5d} {5:>5d} {1:>5d} *"
        assert lines[3] == (
            f"{5:>5d} {2:>5d} {3:>5d} {4:>5d} {1:>5d} "
            "# This is a comment."
        )
        assert lines[4] == (
            f"{5:>5d} {2:>5d} {3:>5d} {4:>5d} {1:>5d} * # "
            "This is a comment."
        )
        assert lines[5] == "END"

    def test__get_shake_lines(self):
        """
        Test the _get_shake_lines method.
        """

        topology = BondedTopology()

        lines = TopologyFileWriter._get_shake_lines(topology)
        assert lines[0] == "SHAKE 0 0 0"
        assert lines[1] == "END"

        bond1 = Bond(
            index1=1,
            index2=2,
            equilibrium_distance=1.4,
            is_shake=True,
        )
        bond2 = Bond(
            index1=2,
            index2=3,
            equilibrium_distance=1.4,
            is_linker=True,
            is_shake=True
        )
        bond3 = Bond(
            index1=5,
            index2=2,
            equilibrium_distance=1.5,
            is_shake=True,
            comment="This is a comment.",
        )
        bond4 = Bond(
            index1=5,
            index2=2,
            equilibrium_distance=1.5,
            is_linker=True,
            is_shake=True,
            comment="This is a comment."
        )

        topology.shake_bonds = [bond1, bond2, bond3, bond4]

        lines = TopologyFileWriter._get_shake_lines(topology)
        assert lines[0] == "SHAKE 3 2 2"
        assert lines[1] == f"{1:>5d} {2:>5d} {1.4:>16.12f}\t"
        assert lines[2] == f"{2:>5d} {3:>5d} {1.4:>16.12f}\t*"
        assert lines[3] == (
            f"{5:>5d} {2:>5d} {1.5:>16.12f}\t "
            "# This is a comment."
        )
        assert lines[4] == (
            f"{5:>5d} {2:>5d} {1.5:>16.12f}\t* # "
            "This is a comment."
        )
        assert lines[5] == "END"
