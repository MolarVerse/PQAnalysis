import numpy as np
import pytest


from PQAnalysis.analysis.thermal_expansion.exceptions import ThermalExpansionError
from PQAnalysis.analysis import ThermalExpansion
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError
from PQAnalysis.physical_data.box import Box

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access


class TestThermalExpansion:

    def test__init__type_checking(self, caplog):
        box1 = [Box(), Box(), Box(), Box(), Box()]
        temp1 = np.linspace(248.15, 348.15, 25)
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "temperature_points",
                "Np1DNumberArray",
                "None"
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=None,
            box=box1,
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "box",
                "List[Box]",
                "None"
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=temp1,
            box=None,
        )

    def test__init(self, caplog):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box = Box(a, b, c, alpha, beta, gamma)
        box1 = [box, box, box, box, box]

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points must be provided",
            box=box1,
            temperature_points=None,
        )

        box_not_same_length = [box, box, box, box]
        temperature_points = np.array([248.15, 273.15, 298.15, 323.15, 348.15])
        temperature_points_not_same_step = np.array(
            [248.15, 276.15, 298.15, 323.15, 348.15])
        temperature_points_not_5_points = np.array(
            [248.15, 273.15, 298.15, 323.15, 348.15, 373.15])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Box data must be provided",
            box=None,
            temperature_points=temperature_points
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points must have the same step size",
            box=box1,
            temperature_points=temperature_points_not_same_step
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points must have the same length as the box data",
            box=box_not_same_length,
            temperature_points=temperature_points,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test=(
                "The number of temperature points must be 5",
                "You have provided 6 points",
                "Only 5-point stencil is supported at the moment!"
            ),
            box=box1,
            temperature_points=temperature_points_not_5_points
        )

    def test__init__dummy_implementation(self):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box = Box(a, b, c, alpha, beta, gamma)

        temperature_points = np.linspace(248.15, 348.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, box=[box, box, box, box, box])

        assert thermal_expansion._temperature_step_size == 4
        assert np.all(thermal_expansion._box_av == np.zeros(4))
        assert np.all(thermal_expansion._box_std == np.zeros(4))
        assert np.all(thermal_expansion._thermal_expansions == np.zeros(4))
        assert np.all(thermal_expansion._middle_points == np.zeros(4))

    def test__init__temperature_points(self):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box = Box(a, b, c, alpha, beta, gamma)

        temperature_points = np.linspace(248.15, 348.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, box=[box, box, box, box, box])

        assert np.all(thermal_expansion._temperature_points ==
                      temperature_points)

    def test__init__box(self, caplog):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box = Box(a, b, c, alpha, beta, gamma)

        temperature_points = np.linspace(248.15, 348.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, box=[box, box, box, box, box])

        assert np.all(thermal_expansion._box == [box, box, box, box, box])

    def test__initialize_run(self):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box = Box(a, b, c, alpha, beta, gamma)
        v = a*b*c*np.sqrt(1-np.cos(alpha)**2 -
                          np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

        temperature_points = np.linspace(248.15, 348.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, box=[box, box, box, box, box])

        thermal_expansion._initialize_run()

        a_avg = np.mean(a)
        b_avg = np.mean(b)
        c_avg = np.mean(c)

        a_std = np.std(a)
        b_std = np.std(b)
        c_std = np.std(c)

        v_avg = np.mean(v)
        v_std = np.std(v)

        assert np.all(thermal_expansion.box_avg ==
                      [a_avg, b_avg, c_avg, v_avg])
        assert np.all(thermal_expansion.box_std ==
                      [a_std, b_std, c_std, v_std])

        middle_points = np.array([a_avg[2], b_avg[2], c_avg[2], v_avg[2]])

        assert np.all(thermal_expansion.middle_points == middle_points)
