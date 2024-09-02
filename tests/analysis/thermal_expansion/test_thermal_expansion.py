import numpy as np
import pytest

from beartype.typing import List

from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.analysis.thermal_expansion.exceptions import ThermalExpansionError, ThermalExpansionWarning
from PQAnalysis.analysis import ThermalExpansion
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError
from PQAnalysis.core.cell import Cells, Cell

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access

np.random.seed(0)
a = np.random.normal(9.2, 0.5, 1000)
b = np.random.normal(9.2, 0.5, 1000)
c = np.random.normal(11.4, 0.5, 1000)
alpha = np.ones(1000) * 90
beta = np.ones(1000) * 90
gamma = np.ones(1000) * 90
a_avg = np.array([np.average(a)] * 5)
b_avg = np.array([np.average(b)] * 5)
c_avg = np.array([np.average(c)] * 5)
v_avg = np.array([np.average(a * b * c)] * 5)
a_std = np.array([np.std(a)] * 5)
b_std = np.array([np.std(b)] * 5)
c_std = np.array([np.std(c)] * 5)
v_std = np.array([np.std(a * b * c)] * 5)
box_avg = [a_avg, b_avg, c_avg, v_avg]
box_std = [a_std, b_std, c_std, v_std]
temp = np.array([248.15, 273.15, 298.15, 323.15, 348.15])
temp_to_many = np.array([248.15, 273.15, 298.15, 323.15, 348.15, 373.15])
box_not_four_columns = [a_avg, b_avg, c_avg, v_avg, v_avg]
a_avg_not_five_points = np.array([np.average(a)] * 6)
b_avg_not_five_points = np.array([np.average(b)] * 6)
c_avg_not_five_points = np.array([np.average(c)] * 6)
v_avg_not_five_points = np.array([np.average(a * b * c)] * 6)
box_not_five_points = [
    a_avg_not_five_points,
    b_avg_not_five_points,
    c_avg_not_five_points,
    v_avg_not_five_points
]
box_empty_avg = [np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5)]
box_empty_std = [np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5)]



class TestThermalExpansion:

    def test__init__type_checking(self, caplog):

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "temperature_points", 1.0, Np1DNumberArray | None
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=1.0,
            boxes_avg=box_avg,
            boxes_std=box_std,
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "boxes_avg", 1.0, List[Np1DNumberArray] | None
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=temp,
            boxes_avg=1.0,
            boxes_std=box_std,
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "boxes_std", 1.0, List[Np1DNumberArray] | None
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=temp,
            boxes_avg=box_avg,
            boxes_std=1.0,
        )

    def test__init__(self, caplog):

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points must be provided",
            boxes_avg=box_avg,
            boxes_std=box_std,
            temperature_points=None,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Cell data must be provided",
            boxes_avg=None,
            boxes_std=box_std,
            temperature_points=temp
        )
        # assert_logging_with_exception(
        #     caplog=caplog,
        #     logging_name=ThermalExpansion.__qualname__,
        #     logging_level="WARNING",
        #     exception=ThermalExpansionWarning,
        #     function=ThermalExpansion,
        #     message_to_test="The standard deviation of the boxes is set to zero because no data is provided,
        #     boxes_avg=box_avg,
        #     boxes_std=None,
        #     temperature_points=temp
        # )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test=
            "The number of temperature points must be 5. You have provided 6 points. Only 5-point stencil is supported at the moment!",
            boxes_avg=box_avg,
            boxes_std=box_std,
            temperature_points=temp_to_many
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points must have the same step size",
            boxes_avg=box_avg,
            boxes_std=box_std,
            temperature_points=np.array(
                [248.15, 276.15, 298.15, 323.15, 349.15]
            )
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="The boxes data must have 4 columns",
            boxes_avg=box_not_four_columns,
            boxes_std=box_not_four_columns,
            temperature_points=temp
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test=
            "The boxes data must have the same length as the temperature points",
            boxes_avg=box_not_five_points,
            boxes_std=box_not_five_points,
            temperature_points=temp
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="The boxes std data must have 4 columns",
            boxes_avg=box_avg,
            boxes_std=box_not_four_columns,
            temperature_points=temp
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test=
            "The boxes std data must have the same length as the temperature points",
            boxes_avg=box_avg,
            boxes_std=box_not_five_points,
            temperature_points=temp
        )

    def test__init__dummy_implementation(self):

        thermal_expansion = ThermalExpansion(
            temperature_points=temp,
            boxes_avg=box_empty_avg,
            boxes_std=box_empty_std
        )

        assert np.allclose(
            thermal_expansion.temperature_step_size, 25, rtol=1e-4
        )
        assert np.all(thermal_expansion.boxes_avg == [np.zeros(5)] * 4)
        assert np.all(thermal_expansion.boxes_std == [np.zeros(5)] * 4)
        assert np.all(thermal_expansion.thermal_expansions == np.zeros(4))
        assert np.all(thermal_expansion.middle_points == np.zeros(4))

    def test__init__temperature_points(self):
        temperature_points = np.array([248.15, 273.15, 298.15, 323.15, 348.15])
        thermal_expansion = ThermalExpansion(
            temperature_points=temp,
            boxes_avg=box_empty_avg,
            boxes_std=box_empty_std
        )

        assert np.all(
            thermal_expansion._temperature_points == temperature_points
        )

    def test__init__box(self, caplog):

        thermal_expansion = ThermalExpansion(
            temperature_points=temp, boxes_avg=box_avg, boxes_std=box_std
        )
        assert np.all(thermal_expansion.boxes_avg == box_avg)
        assert np.all(thermal_expansion.boxes_std == box_std)

    def test__initialize_run(self):

        thermal_expansion = ThermalExpansion(
            temperature_points=temp, boxes_avg=box_avg, boxes_std=box_std
        )
        thermal_expansion._initialize_run()

        average_ref = [a_avg, b_avg, c_avg, v_avg]
        std_ref = [a_std, b_std, c_std, v_std]

        print("THErmal expansion\n", thermal_expansion.boxes_avg)
        print("THErmal expansion\n", average_ref)
        assert np.allclose(thermal_expansion.boxes_avg, average_ref, rtol=1e-6)
        assert np.allclose(thermal_expansion.boxes_std, std_ref, rtol=1e-6)

        middle_points = np.array([a_avg[2], b_avg[2], c_avg[2], v_avg[2]])

        assert np.allclose(
            thermal_expansion.middle_points, middle_points, rtol=1e-4
        )

    def test_run(self):
        box_avg_ref = np.array([a_avg, b_avg, c_avg, v_avg])
        print(box_avg_ref.shape)
        print(box_avg_ref[0])
        reference_thermal_expansions = (
            (
                box_avg_ref[:, 0] - 8 * box_avg_ref[:, 1] +
                8 * box_avg_ref[:, 3] - box_avg_ref[:, 4]
            ) / (12 * 25)
        ) / box_avg_ref[:, 2]
        thermal_expansion = ThermalExpansion(
            temperature_points=temp, boxes_avg=box_avg, boxes_std=box_std
        )
        thermal_expansion._initialize_run()
        thermal_expansion.run()

        assert np.allclose(
            thermal_expansion.thermal_expansions,
            reference_thermal_expansions,
            rtol=1e-6
        )
