import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.io import MoldescriptorReader
from PQAnalysis.io.exceptions import MoldescriptorReaderError
from PQAnalysis.core import Element
from PQAnalysis.exceptions import PQFileNotFoundError



class TestMoldescriptorReader:

    @pytest.mark.parametrize(
        "example_dir",
        ["readMoldescriptor"],
        indirect=False
    )
    def test__init__(self, test_with_data_dir):
        reader = MoldescriptorReader("moldescriptor.dat")
        assert reader.filename == "moldescriptor.dat"

        with pytest.raises(PQFileNotFoundError) as exception:
            MoldescriptorReader("tmp")
        assert str(exception.value) == "File tmp not found."

    def test_read_mol_type(self):
        mol_type = f"""H2o 3 1.0
H 1 2.3
H 1 2.3
O 2 -4.6
"""

        mol_type = MoldescriptorReader._read_mol_type(mol_type.splitlines(), 1)

        assert mol_type.name == "H2o"
        assert mol_type.n_atoms == 3
        assert mol_type.total_charge == 1.0
        assert mol_type.elements == [Element("H"), Element("H"), Element("O")]
        assert np.allclose(mol_type.atom_types, np.array([1, 1, 2]))
        assert np.allclose(
            mol_type.partial_charges,
            np.array([2.3,
            2.3,
            -4.6])
        )
        assert mol_type.id == 1

        mol_type = f"""H2o 3 1.0
H 1 2.3
H 1
O 2 -4.6
"""

        with pytest.raises(MoldescriptorReaderError) as exception:
            MoldescriptorReader._read_mol_type(mol_type.splitlines(), 1)
        assert str(
            exception.value
        ) == "The number of columns in the body of a mol type must be 3 or 4.\nGot 2 columns instead in text: 'H 1'\n"

    @pytest.mark.parametrize(
        "example_dir",
        ["readMoldescriptor"],
        indirect=False
    )
    def test_read(self, test_with_data_dir):
        reader = MoldescriptorReader("moldescriptor.dat")
        mol_types = reader.read()

        assert len(mol_types) == 2
        water_type = mol_types[0]

        assert water_type.name == "H2O"
        assert water_type.n_atoms == 3
        assert np.allclose(water_type.total_charge, 0.0)
        assert water_type.elements == [
            Element("O"),
            Element("H"),
            Element("H")
        ]
        assert np.allclose(water_type.atom_types, np.array([0, 1, 1]))
        assert np.allclose(
            water_type.partial_charges,
            np.array([-0.65966,
            0.32983,
            0.32983])
        )
        assert water_type.id == 1

        ammonia_type = mol_types[1]

        assert ammonia_type.name == "Ammonia"
        assert ammonia_type.n_atoms == 4
        assert np.allclose(ammonia_type.total_charge, 0.0)
        assert ammonia_type.elements == [
            Element("N"),
            Element("H"),
            Element("H"),
            Element("H")
        ]
        assert np.allclose(ammonia_type.atom_types, np.array([0, 1, 1, 1]))
        assert np.allclose(
            ammonia_type.partial_charges,
            np.array([-0.8022,
            0.2674,
            0.2674,
            0.2674])
        )
        assert ammonia_type.id == 2

        reader = MoldescriptorReader("moldescriptor_withError.dat")
        with pytest.raises(MoldescriptorReaderError) as exception:
            reader.read()
        assert str(
            exception.value
        ) == "The number of columns in the header of a mol type must be 3.\nGot 2 columns instead in text: '  H2O            3'\n"
