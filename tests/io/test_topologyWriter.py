"""
Unit tests for the topologyWriter module.
"""

from PQAnalysis.io import TopologyFileWriter
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
