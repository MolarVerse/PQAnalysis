# 3rd Party library Imports
import pytest
import numpy as np

# import topology marker
from . import pytestmark

# Local Imports
from PQAnalysis.topology.topology import _find_residue_by_id
from PQAnalysis.topology import Topology, TopologyError, Residue, ResidueError
from PQAnalysis.core import Atom


def test_find_residue_by_id():
    residue1 = Residue(name="name", id=0, total_charge=0.0, elements=[
    ], atom_types=np.array([]), partial_charges=np.array([]))
    residue2 = Residue(name="name", id=1, total_charge=0.0, elements=[
    ], atom_types=np.array([]), partial_charges=np.array([]))
    residue3 = Residue(name="name", id=0, total_charge=0.0, elements=[
    ], atom_types=np.array([]), partial_charges=np.array([]))

    with pytest.raises(ResidueError) as exception:
        _find_residue_by_id(2, [residue1, residue2, residue3])
    assert str(exception.value) == "The residue id 2 was not found."

    with pytest.raises(ResidueError) as exception:
        _find_residue_by_id(0, [residue1, residue2, residue3])
    assert str(exception.value) == "The residue id 0 is not unique."

    residue = _find_residue_by_id(1, [residue1, residue2, residue3])
    assert residue == residue2


class TestTopology:
    atoms = [Atom('C'), Atom('H'), Atom('H')]

    def test__init__(self):
        residue_ids = np.array([0, 0, 0, 1])

        with pytest.raises(TopologyError) as exception:
            Topology(atoms=self.atoms, residue_ids=residue_ids)
        assert str(
            exception.value) == "The number of atoms does not match the number of residue ids."

        topology = Topology(atoms=self.atoms)
        assert topology.n_atoms == 3
        assert topology.atoms == self.atoms
        assert np.all(topology.residue_ids == np.array([0, 0, 0]))
        assert topology.atomtype_names == ['C', 'H', 'H']
        assert topology.reference_residues == []
        assert topology.residues == []

    # def test_setup_residues(self):
    #     residue_ids = np.array([0, 1, 0])

    #     topology = Topology(atoms=self.atoms)
    #     residues = topology.setup_residues(residue_ids, self.atoms)
    #     assert residues == []

    #     # residues is here empty because no reference residues are set - otherwise it would throw an error!
    #     topology = Topology(atoms=self.atoms, residue_ids=residue_ids)
    #     residues = topology.setup_residues(residue_ids, self.atoms)
    #     assert residues == []
