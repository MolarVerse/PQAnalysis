"""
Unit tests for the topologyReader module.
"""

import pytest

from PQAnalysis.io.topology_file import (
    TopologyFileReader,
    TopologyFileWriter,
    TopologyFileError,
)
from PQAnalysis.topology import BondedTopology
from PQAnalysis.topology.bonded_topology import (
    Bond,
    Angle,
    Dihedral,
    JCoupling,
    DistanceConstraint,
)

from . import pytestmark  # pylint: disable=unused-import



class TestTopologyFileReader:

    """
    Test cases for the TopologyFileReader class.
    """

    def test_read_j_couplings(self, tmp_path):
        """
        Test that a j_couplings block is parsed into JCoupling objects.
        """
        filename = str(tmp_path / "topology.top")
        content = (
            "J_COUPLINGS 2 2 2 2\n"
            "1 2 3 4 1\n"
            "5 6 7 8 2 # a comment\n"
            "END\n"
        )

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        topology = TopologyFileReader(filename).read()

        assert topology.j_couplings == [
            JCoupling(index1=1, index2=2, index3=3, index4=4,
                      j_coupling_type=1),
            JCoupling(index1=5, index2=6, index3=7, index4=8,
                      j_coupling_type=2),
        ]
        assert topology.j_couplings[1].comment == "a comment"

    def test_read_j_couplings_invalid_number_of_columns(self, tmp_path):
        """
        Test that a malformed j-coupling line raises a TopologyFileError.
        """
        filename = str(tmp_path / "topology.top")
        content = "J_COUPLINGS 1 1 1 1\n1 2 3 4\nEND\n"

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileReader(filename).read()
        assert str(exception.value) == (
            "Invalid number of columns in j-coupling block. Expected 5."
        )

    def test_read_j_couplings_duplicate_atoms(self, tmp_path):
        """
        Test that repeated atom indices in a j-coupling raise.
        """
        filename = str(tmp_path / "topology.top")
        content = "J_COUPLINGS 1 1 1 1\n1 2 2 4 1\nEND\n"

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileReader(filename).read()
        assert str(exception.value) == (
            "Atoms in j-coupling block cannot be the same."
        )

    def test_read_distance_constraints(self, tmp_path):
        """
        Test that a dist_constraints block is parsed into
        DistanceConstraint objects.
        """
        filename = str(tmp_path / "topology.top")
        content = (
            "DIST_CONSTRAINTS 2 2\n"
            "1 2 1.0 2.0 100.0 -0.1\n"
            "3 4 1.5 2.5 200.0 0.0 # a comment\n"
            "END\n"
        )

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        topology = TopologyFileReader(filename).read()

        assert topology.distance_constraints == [
            DistanceConstraint(
                index1=1,
                index2=2,
                lower_distance=1.0,
                upper_distance=2.0,
                spring_constant=100.0,
                d_spring_constant_dt=-0.1
            ),
            DistanceConstraint(
                index1=3,
                index2=4,
                lower_distance=1.5,
                upper_distance=2.5,
                spring_constant=200.0,
                d_spring_constant_dt=0.0
            ),
        ]
        assert topology.distance_constraints[1].comment == "a comment"

    def test_read_distance_constraints_invalid_number_of_columns(
        self, tmp_path
    ):
        """
        Test that a malformed distance constraint line raises
        a TopologyFileError.
        """
        filename = str(tmp_path / "topology.top")
        content = "DIST_CONSTRAINTS 1 1\n1 2 1.0 2.0 100.0\nEND\n"

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileReader(filename).read()
        assert str(exception.value) == (
            "Invalid number of columns in distance constraints "
            "block. Expected 6."
        )

    def test_read_distance_constraints_same_atoms(self, tmp_path):
        """
        Test that identical atom indices in a distance constraint raise.
        """
        filename = str(tmp_path / "topology.top")
        content = "DIST_CONSTRAINTS 1 1\n1 1 1.0 2.0 100.0 -0.1\nEND\n"

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileReader(filename).read()
        assert str(exception.value) == (
            "Atoms in distance constraints block cannot be the same."
        )

    def test_read_distance_constraints_lower_greater_than_upper(
        self, tmp_path
    ):
        """
        Test that lower > upper distance raises, matching PQ.
        """
        filename = str(tmp_path / "topology.top")
        content = "DIST_CONSTRAINTS 1 1\n1 2 2.0 1.0 100.0 -0.1\nEND\n"

        with open(filename, "w", encoding="utf-8") as file:
            file.write(content)

        with pytest.raises(TopologyFileError) as exception:
            TopologyFileReader(filename).read()
        assert str(exception.value) == (
            "Lower distance cannot be greater than upper distance "
            "in distance constraints block."
        )

    def test_round_trip_all_sections(self, tmp_path):
        """
        Test that a topology file containing all seven sections reads
        into a BondedTopology and writes back identically.
        """
        filename = str(tmp_path / "topology.top")
        filename_rewritten = str(tmp_path / "topology_rewritten.top")

        topology = BondedTopology(
            bonds=[
                Bond(index1=1, index2=2, bond_type=1),
                Bond(index1=3, index2=4, bond_type=2, is_linker=True),
            ],
            angles=[
                Angle(index1=1, index2=2, index3=3, angle_type=1),
                Angle(
                    index1=4, index2=5, index3=6, angle_type=2, is_linker=True
                ),
            ],
            dihedrals=[
                Dihedral(
                    index1=1, index2=2, index3=3, index4=4, dihedral_type=1
                ),
                Dihedral(
                    index1=5,
                    index2=6,
                    index3=7,
                    index4=8,
                    dihedral_type=2,
                    is_linker=True
                ),
            ],
            impropers=[
                Dihedral(
                    index1=1,
                    index2=2,
                    index3=3,
                    index4=4,
                    dihedral_type=1,
                    is_improper=True
                ),
            ],
            shake_bonds=[
                Bond(index1=1, index2=2, equilibrium_distance=1.0,
                     is_shake=True),
                Bond(
                    index1=3,
                    index2=4,
                    equilibrium_distance=1.2,
                    is_linker=True,
                    is_shake=True
                ),
            ],
            j_couplings=[
                JCoupling(
                    index1=1, index2=2, index3=3, index4=4, j_coupling_type=1
                ),
            ],
            distance_constraints=[
                DistanceConstraint(
                    index1=1,
                    index2=2,
                    lower_distance=1.0,
                    upper_distance=2.0,
                    spring_constant=100.0,
                    d_spring_constant_dt=-0.1
                ),
            ]
        )

        TopologyFileWriter(filename).write(topology)

        read_topology = TopologyFileReader(filename).read()

        assert read_topology.bonds == topology.bonds
        assert read_topology.angles == topology.angles
        assert read_topology.dihedrals == topology.dihedrals
        assert read_topology.impropers == topology.impropers
        assert read_topology.shake_bonds == topology.shake_bonds
        assert read_topology.j_couplings == topology.j_couplings
        assert (
            read_topology.distance_constraints ==
            topology.distance_constraints
        )

        TopologyFileWriter(filename_rewritten).write(read_topology)

        with open(filename, "r", encoding="utf-8") as file:
            content = file.read()

        with open(filename_rewritten, "r", encoding="utf-8") as file:
            content_rewritten = file.read()

        assert content == content_rewritten
