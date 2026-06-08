import pytest

from PQAnalysis.exceptions import PQValueError
from PQAnalysis.utils import calculate_simulation_time
from PQAnalysis.utils.units import calculate_simulation_time as calculate_time

from . import pytestmark  # pylint: disable=unused-import



def test_calculate_simulation_time():
    assert calculate_simulation_time(10000, 0.5) == 5.0
    assert calculate_simulation_time(10000, 0.5, unit="ps") == 5.0
    assert calculate_simulation_time(10000, 0.5, unit="fs") == 5000.0
    assert calculate_simulation_time(10000, 0.5, unit="PS") == 5.0
    assert calculate_simulation_time(0, 0.5) == 0.0
    assert calculate_time(250, 2.0, unit="fs") == 500.0


def test_calculate_simulation_time_raises_for_invalid_input():
    with pytest.raises(PQValueError) as exception:
        calculate_simulation_time(-1, 0.5)
    assert str(exception.value) == (
        "timesteps has to be greater than or equal to 0."
    )

    with pytest.raises(PQValueError) as exception:
        calculate_simulation_time(1, -0.5)
    assert str(exception.value) == (
        "timestep_fs has to be greater than or equal to 0."
    )

    with pytest.raises(PQValueError) as exception:
        calculate_simulation_time(1, 0.5, unit="ns")
    assert str(exception.value) == (
        "Unsupported time unit ns. Supported units are fs and ps."
    )
