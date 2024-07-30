import numpy as np
import pytest

from beartype.typing import List

from PQAnalysis.types import Np1DNumberArray
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
                1.0,
                Np1DNumberArray | None
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=1.0,
            boxes=box1,
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "boxes",
                1.0,
                List[Box] | None
            ),
            exception=PQTypeError,
            function=ThermalExpansion,
            temperature_points=temp1,
            boxes=1.0,
        )

    def test__init__(self, caplog):
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
            boxes=box1,
            temperature_points=None,
        )

        box_not_same_length = [box, box, box, box]
        temperature_points = np.arange(248.15, 373.15, 25)
        print(temperature_points)
        temperature_points_not_same_step = np.array(
            [248.15, 276.15, 298.15, 323.15, 349.15])
        temperature_points_not_5_points = np.array(
            [248.15, 273.15, 298.15, 323.15, 348.15, 373.15])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Box data must be provided",
            boxes=None,
            temperature_points=temperature_points
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points must have the same step size",
            boxes=box1,
            temperature_points=temperature_points_not_same_step
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test="Temperature points and boxes data must have the same length",
            boxes=box_not_same_length,
            temperature_points=temperature_points,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ThermalExpansion.__qualname__,
            logging_level="ERROR",
            exception=ThermalExpansionError,
            function=ThermalExpansion,
            message_to_test=(
                (
                    "The number of temperature points must be 5. "
                    "You have provided 6 points. "
                    "Only 5-point stencil is supported at the moment!"
                )

            ),
            boxes=box1,
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

        temperature_points = np.arange(248.15, 373.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, boxes=[box, box, box, box, box])

        assert np.allclose(
            thermal_expansion.temperature_step_size, 25, rtol=1e-4)
        assert np.all(thermal_expansion.boxes_avg == np.zeros(4))
        assert np.all(thermal_expansion.boxes_std == np.zeros(4))
        assert np.all(thermal_expansion.thermal_expansions == np.zeros(4))
        assert np.all(thermal_expansion.middle_points == np.zeros(4))

    def test__init__temperature_points(self):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box = Box(a, b, c, alpha, beta, gamma)

        temperature_points = np.arange(248.15, 373.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, boxes=[box, box, box, box, box])

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

        temperature_points = np.arange(248.15, 373.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, boxes=[box, box, box, box, box])

        assert np.all(thermal_expansion.boxes == [box, box, box, box, box])

    def test__initialize_run(self):
        np.random.seed(0)
        a = np.random.normal(9.2, 0.5, 1000)
        b = np.random.normal(9.2, 0.5, 1000)
        c = np.random.normal(11.4, 0.5, 1000)
        alpha = np.ones(1000)*90
        beta = np.ones(1000)*90
        gamma = np.ones(1000)*90
        box1 = Box(a, b, c, alpha, beta, gamma)
        box2 = Box(a*0.1, b*0.1, c*0.1, alpha, beta, gamma)
        box3 = Box(a*0.2, b*0.2, c*0.2, alpha, beta, gamma)
        box4 = Box(a*0.3, b*0.3, c*0.3, alpha, beta, gamma)
        box5 = Box(a*0.4, b*0.4, c*0.4, alpha, beta, gamma)
        v1 = a*b*c*np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 -
                           np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        v2 = a*0.1*b*0.1*c*0.1*np.sqrt(1 - np.cos(alpha)**2 - np.cos(
            beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        v3 = a*0.2*b*0.2*c*0.2*np.sqrt(1 - np.cos(alpha)**2 - np.cos(
            beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        v4 = a*0.3*b*0.3*c*0.3*np.sqrt(1 - np.cos(alpha)**2 - np.cos(
            beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        v5 = a*0.4*b*0.4*c*0.4*np.sqrt(1 - np.cos(alpha)**2 - np.cos(
            beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

        temperature_points = np.arange(248.15, 373.15, 25)
        thermal_expansion = ThermalExpansion(
            temperature_points=temperature_points, boxes=[box1, box2, box3, box4, box5])

        thermal_expansion._initialize_run()

        a_avg = [np.mean(a), np.mean(a*0.1), np.mean(a*0.2),
                 np.mean(a*0.3), np.mean(a*0.4)]
        b_avg = [np.mean(b), np.mean(b*0.1), np.mean(b*0.2),
                 np.mean(b*0.3), np.mean(b*0.4)]
        c_avg = [np.mean(c), np.mean(c*0.1), np.mean(c*0.2),
                 np.mean(c*0.3), np.mean(c*0.4)]
        v_avg = [np.mean(v1), np.mean(v2), np.mean(v3),
                 np.mean(v4), np.mean(v5)]

        a_std = [np.std(a), np.std(a*0.1), np.std(a*0.2),
                 np.std(a*0.3), np.std(a*0.4)]
        b_std = [np.std(b), np.std(b*0.1), np.std(b*0.2),
                 np.std(b*0.3), np.std(b*0.4)]
        c_std = [np.std(c), np.std(c*0.1), np.std(c*0.2),
                 np.std(c*0.3), np.std(c*0.4)]
        v_std = [np.std(v1), np.std(v2), np.std(v3),
                 np.std(v4), np.std(v5)]
        average_ref = [a_avg, b_avg, c_avg, v_avg]
        std_ref = [a_std, b_std, c_std, v_std]

        print("THErmal expansion\n", thermal_expansion.boxes_avg)
        print("THErmal expansion\n", average_ref)
        assert np.allclose(thermal_expansion.boxes_avg,
                           average_ref, rtol=1e-6)
        assert np.allclose(thermal_expansion.boxes_std,
                           std_ref, rtol=1e-6)

        middle_points = np.array([a_avg[2], b_avg[2], c_avg[2], v_avg[2]])

        assert np.allclose(thermal_expansion.middle_points,
                           middle_points, rtol=1e-4)
