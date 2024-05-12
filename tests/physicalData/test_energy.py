from collections import defaultdict

import pytest
import numpy as np

from PQAnalysis.physical_data import Energy, EnergyError



class TestEnergy:

    def test__init__(self):
        data = np.array([1, 2, 3])
        energy = Energy(data)

        assert np.allclose(energy.data, [data])

    def test__setup_info_dictionary(self):
        data = np.array([[1], [2]])
        info = {1: 0, 2: 1}
        units = {1: "a", 2: "b"}

        energy = Energy(data, info=info, units=units)
        assert energy.info == info
        assert energy.units == units
        assert np.allclose(energy.data, data)
        assert energy.info_given == True
        assert energy.units_given == True

        energy = Energy(data)
        assert energy.info == defaultdict(lambda: None)
        assert energy.units == defaultdict(lambda: None)
        assert np.allclose(energy.data, data)
        assert energy.info_given == False
        assert energy.units_given == False

        with pytest.raises(EnergyError) as exception:
            Energy(data, info={1: 0})
        assert str(
            exception.value
        ) == "The length of info dictionary has to be equal to the length of data."

        with pytest.raises(EnergyError) as exception:
            Energy(data, units={1: 0})
        assert str(
            exception.value
        ) == "The length of units dictionary has to be equal to the length of data."

        with pytest.raises(EnergyError) as exception:
            Energy(data, info={1: 0, 2: 0}, units={1: 0, 3: 0})
        assert str(
            exception.value
        ) == "The keys of the info and units dictionary do not match."
