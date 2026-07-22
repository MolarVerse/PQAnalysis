"""
Tests for the pure-Python fallback import of the ADF frame kernel.

The ``adf`` module imports the compiled Cython kernel
:py:mod:`PQAnalysis.analysis.adf._adf_kernel` and falls back to the pure
numpy implementation :py:mod:`PQAnalysis.analysis.adf._adf_kernel_py`
only when the compiled extension is not importable
(``except ModuleNotFoundError``). As the compiled extension is present
in the test environment, that fallback branch is exercised by importing
a *fresh* copy of the ``adf`` module source (from the same file, under a
throwaway name) while a meta_path finder blocks the compiled kernel. A
fresh copy is used instead of reloading the shared module object on
purpose: reloading would rebind the ``ADF`` class to a new object and
break the beartype runtime type checks of other tests that captured the
original class identity. The throwaway module leaves the shared ``adf``
module completely untouched.
"""

import importlib
import importlib.machinery
import importlib.util
import sys

from PQAnalysis.analysis.adf import ADF

from .. import pytestmark  # pylint: disable=unused-import

#: The dotted name of the compiled kernel module blocked below.
_KERNEL = "PQAnalysis.analysis.adf._adf_kernel"

#: The dotted name of the pure-Python fallback kernel module.
_KERNEL_PY = "PQAnalysis.analysis.adf._adf_kernel_py"

#: The throwaway name under which the fresh adf module copy is loaded.
_PROBE = "PQAnalysis.analysis.adf._adf_fallback_probe"

# pylint: disable=protected-access



class _BlockKernel:

    """
    A meta_path finder that raises ModuleNotFoundError for the compiled
    ADF kernel module and defers to the remaining finders for every
    other module.
    """

    def find_spec(self, name, path=None, target=None):
        """
        Blocks the compiled ADF kernel and defers otherwise.

        Parameters
        ----------
        name : str
            The dotted name of the module being imported.
        path : optional
            The parent package search path (unused).
        target : optional
            The reload target module (unused).

        Returns
        -------
        None
            Always None so that the remaining finders resolve the module.

        Raises
        ------
        ModuleNotFoundError
            If the compiled ADF kernel module is requested.
        """
        if name == _KERNEL:
            raise ModuleNotFoundError(name)
        return None



def test_adf_frame_histogram_falls_back_to_python_kernel():
    """
    Exercises the ``except ModuleNotFoundError`` fallback import of the
    ADF frame kernel without disturbing the shared ``adf`` module.
    """
    real_module = sys.modules[ADF.__module__]
    real_kernel = real_module.adf_frame_histogram

    assert real_kernel.__module__.endswith("_adf_kernel")

    # pre-import the pure-Python fallback kernel so the probe's fallback
    # ``from ._adf_kernel_py import ...`` is served from the module cache
    # instead of re-running the import machinery (which, under the
    # beartype claw import hook in DEBUG mode, recurses on a re-entrant
    # first import of a hooked package module)
    importlib.import_module(_KERNEL_PY)

    blocker = _BlockKernel()
    had_kernel = _KERNEL in sys.modules

    # drop the cached compiled kernel so the fresh import re-runs the
    # import machinery (and thus the blocking finder) for the relative
    # ``from ._adf_kernel import ...`` line; already bound references
    # (e.g. real_module.adf_frame_histogram) are unaffected by this
    sys.modules.pop(_KERNEL, None)
    sys.meta_path.insert(0, blocker)

    try:
        # exec the adf.py source through a plain SourceFileLoader so the
        # same file lines are executed (and counted by coverage) without
        # going through the beartype claw import hook, which recurses on
        # a re-entrant fresh load of a hooked package module in DEBUG
        # mode. __package__ is set so the relative fallback import
        # ``from ._adf_kernel_py import ...`` resolves correctly.
        loader = importlib.machinery.SourceFileLoader(
            _PROBE, real_module.__file__
        )
        spec = importlib.util.spec_from_loader(_PROBE, loader)
        probe = importlib.util.module_from_spec(spec)
        probe.__package__ = "PQAnalysis.analysis.adf"
        sys.modules[_PROBE] = probe
        loader.exec_module(probe)

        # the ModuleNotFoundError of the compiled kernel selected the
        # pure-Python fallback import (the except branch of adf.py)
        assert probe.adf_frame_histogram.__module__.endswith(
            "_adf_kernel_py"
        )
    finally:
        sys.meta_path.remove(blocker)
        sys.modules.pop(_PROBE, None)
        if had_kernel:
            # restore the compiled kernel in sys.modules for any later
            # test that imports it freshly
            importlib.import_module(_KERNEL)

    # the shared adf module and its class identity are untouched
    assert sys.modules[ADF.__module__] is real_module
    assert real_module.adf_frame_histogram is real_kernel
    assert real_module.adf_frame_histogram.__module__.endswith("_adf_kernel")
