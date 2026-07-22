"""
Tests for the RawChargeTrajectoryReader fast-path reader and the
ModuleNotFoundError fallback import of its scalar line parser.
"""

import importlib
import sys

from PQAnalysis.analysis.vacf._raw_charge_reader import (
    RawChargeTrajectoryReader,
)

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

#: The compiled Cython kernel that provides ``parse_charge_lines``. The
#: fallback import is only taken when this module cannot be imported.
_KERNEL = "PQAnalysis.analysis.vacf._vacf_kernel"

#: The module object of the reader under test, used for the in-process
#: reload that exercises the ModuleNotFoundError fallback import.
_rcr_module = sys.modules[RawChargeTrajectoryReader.__module__]



class _KernelBlocker:

    """
    A meta_path finder that hides the compiled vacf kernel so that the
    ``except ModuleNotFoundError`` fallback import is taken on reload.
    """

    def find_spec(self, name, path=None, target=None):  # pylint: disable=unused-argument
        """Raise ModuleNotFoundError for the kernel, delegate otherwise."""
        if name == _KERNEL:
            raise ModuleNotFoundError(name)



class TestRawChargeTrajectoryReader:

    """
    Tests for the RawChargeTrajectoryReader construction.
    """

    def test_single_filename_is_wrapped_into_a_list(self, tmp_path):
        """
        A single filename string is wrapped into the ``filenames`` list.
        """
        charge_file = tmp_path / "single.chrg"
        charge_file.write_text("1\ncomment\nO 0.5\n", encoding="utf-8")

        reader = RawChargeTrajectoryReader(str(charge_file))

        assert reader.multiple_files is False
        assert reader.filenames == [str(charge_file)]

    def test_module_not_found_fallback_import(self):
        """
        When the compiled kernel is absent, ``parse_charge_lines`` is
        imported from the pure-python ``_vacf_kernel_py`` fallback.

        The reader module is reloaded in-process with the compiled
        kernel hidden by a meta_path blocker and restored to the
        compiled version afterwards, so that the test does not affect
        any other test regardless of the execution order.
        """
        assert _rcr_module.parse_charge_lines.__module__.endswith(
            "_vacf_kernel"
        )

        sys.modules.pop(_KERNEL, None)
        blocker = _KernelBlocker()
        sys.meta_path.insert(0, blocker)

        try:
            importlib.reload(_rcr_module)
            assert _rcr_module.parse_charge_lines.__module__.endswith(
                "_vacf_kernel_py"
            )
        finally:
            sys.meta_path.remove(blocker)
            sys.modules.pop(_KERNEL, None)
            importlib.reload(_rcr_module)

        assert _rcr_module.parse_charge_lines.__module__.endswith(
            "_vacf_kernel"
        )
