import numpy as np
import pytest
from scipy.optimize import curve_fit

from PQAnalysis.analysis.bulk_modulus.exceptions import BulkModulusError, BulkModulusWarning
from PQAnalysis.analysis import BulkModulus
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.types import PositiveReal, Np1DNumberArray, Np2DNumberArray
from PQAnalysis.analysis.bulk_modulus.exceptions import BulkModulusError
from PQAnalysis.type_checking import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception


# pylint: disable=protected-access


class TestBulkModulus:
    def test__init__type_checking(self, caplog):
        volume_equilibrium = 1.0
        volumes_perturbation = np.array([0.8, 1.2])
        pressures_perturbation = np.array([[1000, 1300], [10, 20]])
        mode = "simple"

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                arg_name="volume_equilibrium",
                value=-1,
                expected_type=PositiveReal | None
            ),
            exception=PQTypeError,
            function=BulkModulus,
            volume_equilibrium=-1,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation=pressures_perturbation,
            mode=mode,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                arg_name="pressures_perturbation",
                value="0.8",
                expected_type=Np1DNumberArray | Np2DNumberArray | None

            ),
            exception=PQTypeError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation="0.8",
            mode=mode,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                arg_name="volumes_perturbation",
                value=-1,
                expected_type=Np1DNumberArray | None
            ),
            exception=PQTypeError,
            function=BulkModulus,
            volumes_perturbation=-1,
            volume_equilibrium=volume_equilibrium,
            pressures_perturbation=pressures_perturbation,
            mode=mode,
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                arg_name="mode",
                value=1,
                expected_type=str | None

            ),
            exception=PQTypeError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation=pressures_perturbation,
            mode=1,
        )

    def test__init__(self, caplog):
        volume_equilibrium = 1.0
        volumes_perturbation = np.array([0.8, 1.2])
        pressures_perturbation = np.array([[1000, 1300], [10, 20]])
        mode = "simple"
        print("here")
        print(len(np.array([[1000, 1300], [10, 20]])))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="The equilibrium volume must be provided.",
            exception=BulkModulusError,
            function=BulkModulus,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation=pressures_perturbation,
            mode=mode,
            volume_equilibrium=None
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="The volume perturbation must be provided.",
            exception=BulkModulusError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=None,
            pressures_perturbation=pressures_perturbation,
            mode=mode
        )
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="The pressure perturbation must be provided.",
            exception=BulkModulusError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation=None,
            mode=mode,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="The number of perturbation pressure and volume must be the same.",
            exception=BulkModulusError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation=np.array(
                [[1000, 1300], [10, 20], [1, 2]]),
            mode=mode,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="The pressure perturbation must be a 1D or 2D array.",
            exception=BulkModulusError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=volumes_perturbation,
            pressures_perturbation=np.array(
                [[1000, 1300, 100], [100, 100, 100]]),
            mode=mode,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="At least two perturbation pressures are required.",
            exception=BulkModulusError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=np.array([100]),
            pressures_perturbation=np.array([[1000, 10]]),
            mode=mode,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=BulkModulus.__qualname__,
            logging_level="ERROR",
            message_to_test="The volume perturbation must be sorted from smallest to largest.",
            exception=BulkModulusError,
            function=BulkModulus,
            volume_equilibrium=volume_equilibrium,
            volumes_perturbation=np.array([1.2, 0.8]),
            pressures_perturbation=pressures_perturbation,
            mode=mode,
        )

        # assert_logging_with_exception(
        #     caplog=caplog,
        #     logging_name="BulkModulusWarning",
        #     logging_level="WARNING",
        #     message_to_test="The pressure perturbation is a 1D array. The standard deviation will be set to 0.",
        #     exception=BulkModulusWarning,
        #     function=BulkModulus,
        #     volume_equilibrium=volume_equilibrium,
        #     volumes_perturbation=np.array([0.8, 1.2]),
        #     pressures_perturbation=np.array([1000, 1300]),
        #     mode=mode,
        # )

        # assert_logging_with_exception(
        #     caplog=caplog,
        #     logging_name="BulkModulusWarning",
        #     logging_level="WARNING",
        #     message_to_test="The pressure perturbation is a 2D array. The first column will be used as the average and the second column as the standard deviation.",
        #     exception=BulkModulusWarning,
        #     function=BulkModulus,
        #     volume_equilibrium=volume_equilibrium,
        #     volumes_perturbation=volumes_perturbation,
        #     pressures_perturbation=pressures_perturbation,
        #     mode=mode,
        # )

        def test__initialize_run_simple(self):
            volume_equilibrium = 1.0
            volumes_perturbation = np.array([0.8, 1.2, 1.4])
            pressures_perturbation = np.array([[1000, 1300], [10, 20], [1, 2]])
            mode = "simple"
            _bulk_modulus = BulkModulus(
                pressures_perturbation=pressures_perturbation,
                volumes_perturbation=volumes_perturbation,
                volume_equilibrium=volume_equilibrium,
                mode=mode)

            # assert logging with exception of the _initialize_run_simple

        def test_two_point_stencel(self):
            volume_equilibrium = 1.0
            volumes_perturbation = np.array([0.8, 1.2, 1.4])
            pressures_perturbation = np.array([[1000, 1300], [10, 20], [1, 2]])
            mode = "two_point_stencel"
            _bulk_modulus = BulkModulus(
                pressures_perturbation=pressures_perturbation,
                volumes_perturbation=volumes_perturbation,
                volume_equilibrium=volume_equilibrium,
                mode=mode)

            finite_difference_ref = np.diff(
                pressures_perturbation, axis=0) / np.diff(volumes_perturbation)

            finite_difference = _bulk_modulus._two_point_stencel()

            assert np.allclose(finite_difference, finite_difference_ref)

        def test__calculate_bulk_modulus(self):
            volume_equilibrium = 1.0
            volumes_perturbation = np.array([0.8, 1.2, 1.4])
            pressures_perturbation = np.array([[1000, 1300], [10, 20], [1, 2]])
            mode = "two_point_stencel"
            _bulk_modulus = BulkModulus(
                pressures_perturbation=pressures_perturbation,
                volumes_perturbation=volumes_perturbation,
                volume_equilibrium=volume_equilibrium,
                mode=mode)

            _bulk_modulus._calculate_bulk_modulus_simple()

            bm_ref = -volume_equilibrium * np.diff(
                pressures_perturbation, axis=0) / np.diff(volumes_perturbation)

            assert np.allclose(_bulk_modulus.bulk_modulus, bm_ref)

        def test__run(self):
            mode = "simple"
            volume_equilibrium = 1.0
            volumes_perturbation = np.array([0.8, 1.2, 1.4])
            pressures_perturbation = np.array([[1000, 1300], [10, 20], [1, 2]])
            _bulk_modulus = BulkModulus(
                pressures_perturbation=pressures_perturbation,
                volumes_perturbation=volumes_perturbation,
                volume_equilibrium=volume_equilibrium,
                mode=mode)

            assert _bulk_modulus._mode == mode

            assert _bulk_modulus.run() == _bulk_modulus._run_simple()

            mode = "MEOS"

            _bulk_modulus = BulkModulus(
                pressures_perturbation=pressures_perturbation,
                volumes_perturbation=volumes_perturbation,
                volume_equilibrium=volume_equilibrium,
                mode=mode)

            assert _bulk_modulus._mode == mode

            # call error
            with pytest.raises(BulkModulusError):
                _bulk_modulus.run()

            mode = "BMEOS"

            _bulk_modulus = BulkModulus(
                pressures_perturbation=pressures_perturbation,
                volumes_perturbation=volumes_perturbation,
                volume_equilibrium=volume_equilibrium,
                mode=mode)

            assert _bulk_modulus._mode == mode

            assert _bulk_modulus.run() == _bulk_modulus._run_bmeos()

    def test__fit_bmeos(self):
        def _third_order_bmeos(v, b0, b0_prime, v0):
            return 3/2 * b0 * ((v0/v)**(7/3) - (v0/v)**(5/3)) * (1 + 0.75 * (b0_prime - 4) * ((v0/v)**(2/3) - 1))
        volume_equilibrium = 1.0
        volumes_perturbation = np.array([0.8, 1.2, 1.4])
        pressures_perturbation = np.array([[1000, 20], [10, 20], [1, 2]])
        mode = "BMEOS"

        initial_guesses = [200, 100, volume_equilibrium]
        boundary = ([0, 0, volume_equilibrium-0.00000001], [
                    np.inf, np.inf, volume_equilibrium+0.00000001])

        popt_ref, pcov_ref = curve_fit(
            f=_third_order_bmeos,
            xdata=volumes_perturbation,
            ydata=pressures_perturbation[:, 0],
            sigma=pressures_perturbation[:, 1],
            absolute_sigma=True,
            p0=initial_guesses,
            bounds=boundary
        )
        _bulk_modulus = BulkModulus(
            pressures_perturbation=pressures_perturbation, volumes_perturbation=volumes_perturbation, volume_equilibrium=volume_equilibrium,
            mode=mode)

        b0_opt_ref, b0_prime_opt_ref, v0_opt_ref = popt_ref
        b0_err_ref, b0_prime_err_ref, v0_err_ref = np.sqrt(np.diag(pcov_ref))

        b0_opt, b0_prime_opt, b0_err, b0_prime_err = _bulk_modulus._fit_bmeos()

        assert np.allclose(v0_opt_ref, volume_equilibrium)
        assert np.allclose(b0_opt, b0_opt_ref)
        assert np.allclose(b0_prime_opt, b0_prime_opt_ref)

        assert np.allclose(b0_err, b0_err_ref)
        assert np.allclose(b0_prime_err, b0_prime_err_ref)
