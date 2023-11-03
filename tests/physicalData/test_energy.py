import pytest
import numpy as np

from collections import defaultdict

from PQAnalysis.physicalData.energy import Energy


class TestEnergy:
    def test__init__(self):
        with pytest.raises(TypeError) as exception:
            Energy(1)
        assert str(exception.value) == "data has to be iterable."

        with pytest.raises(ValueError) as exception:
            Energy([[[1]]])
        assert str(exception.value) == "data has to be a 1D or 2D array."

        data = [1, 2, 3]
        energy = Energy(data)

        assert np.allclose(energy.data, [data])

    def test__setup_info_dictionary(self):
        with pytest.raises(TypeError) as exception:
            Energy([1], info=1)
        assert str(exception.value) == "info has to be a dictionary."

        with pytest.raises(TypeError) as exception:
            Energy([1], units=1)
        assert str(exception.value) == "units has to be a dictionary."

        with pytest.raises(ValueError) as exception:
            Energy([1], info={1: 1, 2: 2})
        assert str(
            exception.value) == "The length of info dictionary has to be equal to the length of data."

        with pytest.raises(ValueError) as exception:
            Energy([1], units={1: 1, 2: 2})
        assert str(
            exception.value) == "The length of units dictionary has to be equal to the length of data."

        with pytest.raises(ValueError) as exception:
            Energy([1], info={1: 1}, units={2: 2})
        assert str(
            exception.value) == "The keys of the info and units dictionary do not match."

        data = [[1], [2]]
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
