# pylint: disable=protected-access

# 3rd Party library Imports
import pytest
import numpy as np

# Local Imports
from PQAnalysis.topology.topology import _find_residue_by_id, _unique_residues_
from PQAnalysis.topology import Topology, TopologyError, BondedTopology
from PQAnalysis.core import (
    Atom,
    Element,
    Residue,
    ResidueError,
    QMResidue,
    ResidueWarning
)

# import topology marker
from . import pytestmark  # pylint: disable=unused-import
from ..conftest import assert_logging



def test_find_residue_by_id(caplog):
    residue1 = Residue(
        name="name",
        residue_id=0,
        total_charge=0.0,
        elements=[],
        atom_types=np.array([]),
        partial_charges=np.array([])
    )
    residue2 = Residue(
        name="name",
        residue_id=1,
        total_charge=0.0,
        elements=[],
        atom_types=np.array([]),
        partial_charges=np.array([])
    )
    residue3 = Residue(
        name="name",
        residue_id=0,
        total_charge=0.0,
        elements=[],
        atom_types=np.array([]),
        partial_charges=np.array([])
    )

    assert_logging(
        caplog=caplog,
        logging_name=Topology.__qualname__,
        logging_level="ERROR",
        message_to_test=("The residue id 2 was not found."),
        function=_find_residue_by_id,
        res_id=2,
        residues=[residue1,
        residue2,
        residue3]
    )

    assert_logging(
        caplog=caplog,
        logging_name=Topology.__qualname__,
        logging_level="ERROR",
        message_to_test=("The residue id 0 is not unique."),
        function=_find_residue_by_id,
        res_id=0,
        residues=[residue1,
        residue2,
        residue3]
    )

    residue = _find_residue_by_id(1, [residue1, residue2, residue3])
    assert residue == residue2



def test_unique_residues():
    residues = []
    assert _unique_residues_(residues) == []

    residues = [
        Residue(
        name="name",
        residue_id=0,
        total_charge=0.0,
        elements=[Element("C")],
        atom_types=np.array([0]),
        partial_charges=np.array([0.0])
        )
    ]
    assert _unique_residues_(residues) == residues

    residues.append(
        Residue(
        name="name",
        residue_id=1,
        total_charge=0.0,
        elements=[Element("C")],
        atom_types=np.array([0]),
        partial_charges=np.array([0.0])
        )
    )
    assert _unique_residues_(residues) == residues

    residues.append(residues[0])
    assert _unique_residues_(residues) == residues[:2]

    residues = [QMResidue(Element("C")), QMResidue(Element("C"))]
    assert _unique_residues_(residues) == residues[:1]

    residues.append(QMResidue(Element("H")))
    assert _unique_residues_(residues) == [residues[0], residues[2]]



class TestTopology:
    atoms = [Atom('C'), Atom('H'), Atom('H')]

    def test__init__(self, caplog):
        residue_ids = np.array([0, 0, 0, 1])

        with pytest.raises(TopologyError) as exception:
            Topology(atoms=self.atoms, residue_ids=residue_ids)
        assert str(
            exception.value
        ) == "The number of atoms does not match the number of residue ids."

        topology = Topology(atoms=self.atoms)
        assert topology.n_atoms == 3
        assert topology.atoms == self.atoms
        assert np.all(topology.residue_ids == np.array([0, 0, 0]))
        assert topology.atomtype_names == ['C', 'H', 'H']
        assert topology.reference_residues == []
        assert topology.residues == []

        residue_ids = np.array([0, 1, 1, 0])
        atoms = [Atom('C'), Atom('H'), Atom('H'), Atom('H')]
        reference_residues = [
            Residue(
            name="ALA",
            residue_id=1,
            total_charge=0.0,
            elements=[Element("H"),
            Element("H")],
            atom_types=np.array([0,
            1]),
            partial_charges=np.array([0.1,
            0.1])
            )
        ]
        topology = Topology(
            atoms=atoms,
            residue_ids=residue_ids,
            reference_residues=reference_residues
        )
        assert topology.n_atoms == 4
        assert topology.atoms == atoms
        assert np.all(topology.residue_ids == np.array([0, 1, 1, 0]))
        assert topology.atomtype_names == ['C', 'H', 'H', 'H']
        assert topology.reference_residues == reference_residues
        assert topology.residues == [
            QMResidue(Element('C')),
            reference_residues[0],
            QMResidue(Element('H'))
        ]
        assert np.all(topology.residue_numbers == [0, 1, 1, 2])

        assert_logging(
            caplog,
            Topology.__qualname__,
            "WARNING",
            "There is no check yet if the bonded topology is compatible with the topology. Please make sure that the bonded topology is compatible with the topology!",
            Topology,
            bonded_topology=BondedTopology()
        )

    def test_setup_residues(self, caplog):
        residue_ids = np.array([0, 1, 1])

        topology = Topology(atoms=self.atoms)
        residues, new_atoms = topology._setup_residues(residue_ids, self.atoms)
        assert not residues
        assert new_atoms == self.atoms

        # residues is here empty because no reference residues are set
        # otherwise it would throw an error!
        topology = Topology(atoms=self.atoms, residue_ids=residue_ids)
        residues, new_atoms = topology._setup_residues(residue_ids, self.atoms)
        assert not residues
        assert new_atoms == self.atoms

        atoms = [Atom('C', use_guess_element=False), Atom('H'), Atom('H')]
        topology = Topology(atoms=atoms, residue_ids=residue_ids)
        reference_residues = [
            Residue(
            name="ALA",
            residue_id=1,
            total_charge=0.0,
            elements=[Element("H"),
            Element("H")],
            atom_types=np.array([0,
            1]),
            partial_charges=np.array([0.1,
            0.1])
            )
        ]
        topology.reference_residues = reference_residues

        assert_logging(
            caplog=caplog,
            logging_name=Topology.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The element of atom 0 is not set. If any reference residues are given the "
            "program tries to automatically deduce the residues from the residue ids and "
            "the reference residues. This means that any atom with an unknown element "
            "raises an error. To avoid deducing residue information please set 'check_residues' "
            "to False"
            ),
            function=topology._setup_residues,
            residue_ids=residue_ids,
            atoms=atoms
        )

        atoms = [Atom('C'), Atom('C'), Atom('H')]
        topology = Topology(atoms=atoms, residue_ids=residue_ids)
        topology.reference_residues = reference_residues

        residues, new_atoms = assert_logging(
            caplog=caplog,
            logging_name=Topology.__qualname__,
            logging_level="WARNING",
            message_to_test=(
                "The element of atom 1 (Element(c, 6, 12.0107)) does not match "
                "the element of the reference residue ALA (Element(h, 1, 1.00794)). "
                "Therefore the element type of the residue description will be used "
                "within the topology format!"
            ),
            function=topology._setup_residues,
            residue_ids=residue_ids,
            atoms=atoms,
        )

        assert len(residues) == 2
        assert new_atoms == atoms
        assert residues[0] == QMResidue(Element('C'))
        assert residues[1] == reference_residues[0]

        topology.check_residues = False
        residues, new_atoms = topology._setup_residues(residue_ids, atoms)
        assert residues == []

        residue_ids = np.array([0, 2, 2])
        topology = Topology(atoms=atoms, residue_ids=residue_ids)
        topology.reference_residues = reference_residues

        assert_logging(
            caplog,
            Topology.__qualname__,
            "ERROR",
            "Residue ids [2] have no corresponding reference residue.",
            topology._setup_residues,
            residue_ids=residue_ids,
            atoms=atoms
        )

        residue_ids = np.array([0, 1, 0])
        atoms = [Atom('C'), Atom('H'), Atom('H')]
        topology = Topology(atoms=atoms, residue_ids=residue_ids)
        topology.reference_residues = reference_residues

        assert_logging(
            caplog,
            Topology.__qualname__,
            "ERROR",
            "The residue ids are not contiguous. Problems with residue ALA with indices 1-2.",
            topology._setup_residues,
            residue_ids=residue_ids,
            atoms=atoms
        )

        residue_ids = np.array([1, 1, 0])
        atoms = [Atom('H'), Atom('H'), Atom('C')]
        topology = Topology(atoms=atoms, residue_ids=residue_ids)
        topology.reference_residues = reference_residues
        residues, new_atoms = topology._setup_residues(residue_ids, atoms)
        assert new_atoms == atoms
        assert residues[0] == reference_residues[0]
        assert residues[1] == QMResidue(Element('C'))

        topology = Topology(
            atoms=atoms,
            residue_ids=residue_ids,
            reference_residues=reference_residues
        )
        residues, new_atoms = topology._setup_residues(residue_ids, atoms)
        assert new_atoms == atoms
        assert residues[0] == reference_residues[0]
        assert residues[1] == QMResidue(Element('C'))

    def test__eq__(self):
        assert Topology() != Atom('C')
        assert Topology() != Topology(atoms=self.atoms)
        assert Topology(atoms=self.atoms) == Topology(
            atoms=self.atoms,
            residue_ids=np.array([0,
            0,
            0])
        )
        assert Topology(atoms=self.atoms) != Topology(
            atoms=self.atoms,
            residue_ids=np.array([0,
            0,
            1])
        )

        reference_residues = [
            Residue(
            name="ALA",
            residue_id=1,
            total_charge=0.0,
            elements=[Element("H"),
            Element("H")],
            atom_types=np.array([0,
            1]),
            partial_charges=np.array([0.1,
            0.1])
            )
        ]
        assert Topology(atoms=self.atoms) == Topology(
            atoms=self.atoms,
            reference_residues=reference_residues
        )

    def test__getitem__(self):
        topology = Topology()
        assert topology[0] == Topology()

        topology = Topology(atoms=self.atoms, residue_ids=np.array([0, 1, 0]))
        assert topology[0] == Topology(
            atoms=[Atom('C')],
            residue_ids=np.array([0])
        )
        assert topology[1] == Topology(
            atoms=[Atom('H')],
            residue_ids=np.array([1])
        )
        assert topology[np.array([0,
            1])] == Topology(
            atoms=[Atom('C'),
            Atom('H')],
            residue_ids=np.array([0,
            1])
            )

    def test__str__(self):
        reference_residues = [
            Residue(
            name="ALA",
            residue_id=1,
            total_charge=0.0,
            elements=[Element("H"),
            Element("H")],
            atom_types=np.array([0,
            1]),
            partial_charges=np.array([0.1,
            0.1])
            )
        ]

        topology = Topology(atoms=self.atoms, residue_ids=np.array([0, 1, 1]))
        assert str(
            topology
        ) == "Topology with 3 atoms and 0 residues and 0 unique residues."
        assert str(topology) == repr(topology)
        assert topology.n_mm_residues == 0

        topology = Topology(
            atoms=self.atoms,
            residue_ids=np.array([0,
            1,
            1]),
            reference_residues=reference_residues
        )
        assert str(
            topology
        ) == "Topology with 3 atoms and 2 residues (1 QM residues) and 2 unique residues."
        assert str(topology) == repr(topology)
        assert topology.n_mm_residues == 1

        topology = Topology(
            atoms=self.atoms,
            residue_ids=np.array([0,
            0,
            0]),
            reference_residues=reference_residues
        )

        topology_str = (
            "Topology with 3 atoms and 3 residues"
            " (3 QM residues) and 2 unique residues."
        )
        assert str(topology) == topology_str
        assert str(topology) == repr(topology)
        assert topology.n_mm_residues == 0
