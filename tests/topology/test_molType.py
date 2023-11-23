import pytest
import numpy as np

from PQAnalysis.topology import MolType
from PQAnalysis.core import Atom


class TestMolType:
    def test__init__(self):
        mol_type = MolType(name="name", id=0, total_charge=0.0,
                           elements=[], atom_types=np.array([]), partial_charges=np.array([]))
        assert mol_type.name == "name"
        assert mol_type.id == 0
        assert np.allclose(mol_type.total_charge, 0.0)
        assert mol_type.elements == []
        assert len(mol_type.atom_types) == 0
        assert len(mol_type.partial_charges) == 0
        assert mol_type.n_atoms == 0

        mol_type = MolType(name="name", id=0, total_charge=0.0,
                           elements=[Atom("C"), Atom("H"), Atom("H")], atom_types=np.array([0, 1, 1]),
                           partial_charges=np.array([0.0, 0.1, 0.1]))

        assert mol_type.name == "name"
        assert mol_type.id == 0
        assert np.allclose(mol_type.total_charge, 0.0)
        assert mol_type.elements == [Atom("C"), Atom("H"), Atom("H")]
        assert np.allclose(mol_type.atom_types, np.array([0, 1, 1]))
        assert np.allclose(mol_type.partial_charges,
                           np.array([0.0, 0.1, 0.1]))
        assert mol_type.n_atoms == 3

        mol_type.elements = [Atom("C"), Atom("H"), Atom("O")]
        assert mol_type.elements == [Atom("C"), Atom("H"), Atom("O")]

        mol_type.atom_types = np.array([0, 1, 2])
        assert np.allclose(mol_type.atom_types, np.array([0, 1, 2]))

        mol_type.partial_charges = np.array([0.0, 0.1, 0.2])
        assert np.allclose(mol_type.partial_charges,
                           np.array([0.0, 0.1, 0.2]))

        with pytest.raises(ValueError) as exception:
            mol_type.elements = [Atom("C"), Atom("H")]
        assert str(
            exception.value) == "The number of elements must be the same as the number of atoms."

        with pytest.raises(ValueError) as exception:
            mol_type.atom_types = np.array([0, 1])
        assert str(
            exception.value) == "The number of atom_types must be the same as the number of atoms."

        with pytest.raises(ValueError) as exception:
            mol_type.partial_charges = np.array([0.0, 0.1])
        assert str(
            exception.value) == "The number of partial_charges must be the same as the number of atoms."

        with pytest.raises(ValueError) as exception:
            mol_type = MolType(name="name", id=0, total_charge=0.0,
                               elements=[Atom("C"), Atom("H"), Atom("H")], atom_types=np.array([0, 1]),
                               partial_charges=np.array([0.0, 0.1, 0.1]))
        assert str(
            exception.value) == "The number of elements, atom_types and partial_charges must be the same."
