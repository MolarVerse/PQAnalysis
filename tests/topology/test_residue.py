import pytest
import numpy as np

from PQAnalysis.topology import Residue, ResidueError
from PQAnalysis.core import Element


class TestResidue:
    def test__init__(self):
        residue = Residue(name="name", id=0, total_charge=0.0,
                          elements=[], atom_types=np.array([]), partial_charges=np.array([]))
        assert residue.name == "name"
        assert residue.id == 0
        assert np.allclose(residue.total_charge, 0.0)
        assert residue.elements == []
        assert len(residue.atom_types) == 0
        assert len(residue.partial_charges) == 0
        assert residue.n_atoms == 0

        residue = Residue(name="name", id=0, total_charge=0.0,
                          elements=[Element("C"), Element("H"), Element("H")], atom_types=np.array([0, 1, 1]),
                          partial_charges=np.array([0.0, 0.1, 0.1]))

        assert residue.name == "name"
        assert residue.id == 0
        assert np.allclose(residue.total_charge, 0.0)
        assert residue.elements == [Element("C"), Element("H"), Element("H")]
        assert np.allclose(residue.atom_types, np.array([0, 1, 1]))
        assert np.allclose(residue.partial_charges,
                           np.array([0.0, 0.1, 0.1]))
        assert residue.n_atoms == 3

        residue.elements = [Element("C"), Element("H"), Element("O")]
        assert residue.elements == [Element("C"), Element("H"), Element("O")]

        residue.atom_types = np.array([0, 1, 2])
        assert np.allclose(residue.atom_types, np.array([0, 1, 2]))

        residue.partial_charges = np.array([0.0, 0.1, 0.2])
        assert np.allclose(residue.partial_charges,
                           np.array([0.0, 0.1, 0.2]))

        with pytest.raises(ResidueError) as exception:
            residue.elements = [Element("C"), Element("H")]
        assert str(
            exception.value) == "The number of elements must be the same as the number of atoms."

        with pytest.raises(ResidueError) as exception:
            residue.atom_types = np.array([0, 1])
        assert str(
            exception.value) == "The number of atom_types must be the same as the number of atoms."

        with pytest.raises(ResidueError) as exception:
            residue.partial_charges = np.array([0.0, 0.1])
        assert str(
            exception.value) == "The number of partial_charges must be the same as the number of atoms."

        with pytest.raises(ResidueError) as exception:
            residue = Residue(name="name", id=0, total_charge=0.0,
                              elements=[Element("C"), Element("H"), Element("H")], atom_types=np.array([0, 1]),
                              partial_charges=np.array([0.0, 0.1, 0.1]))
        assert str(
            exception.value) == "The number of elements, atom_types and partial_charges must be the same."