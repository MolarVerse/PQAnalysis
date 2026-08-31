"""Coverage of the momentum module's pure-Python kernel fallback."""

import importlib
import importlib.machinery
import importlib.util
import sys

from PQAnalysis.analysis.momentum import Momentum

from .. import pytestmark  # pylint: disable=unused-import

_KERNEL = "PQAnalysis.analysis.momentum._momentum_kernel"
_KERNEL_PY = "PQAnalysis.analysis.momentum._momentum_kernel_py"
_PROBE = "PQAnalysis.analysis.momentum._momentum_fallback_probe"



class _BlockKernel:

    """Hide the compiled kernel while the throwaway module is imported."""

    def find_spec(self, name, path=None, target=None):
        del path, target
        if name == _KERNEL:
            raise ModuleNotFoundError(name)
        return None



def test_momentum_module_falls_back_to_python_kernel():
    """The module remains usable when its compiled extension is absent."""
    real_module = sys.modules[Momentum.__module__]
    real_kernel = real_module.legacy_momentum_norm

    assert real_kernel.__module__.endswith("_momentum_kernel")

    importlib.import_module(_KERNEL_PY)
    blocker = _BlockKernel()
    had_kernel = _KERNEL in sys.modules
    sys.modules.pop(_KERNEL, None)
    sys.meta_path.insert(0, blocker)

    try:
        loader = importlib.machinery.SourceFileLoader(
            _PROBE,
            real_module.__file__,
        )
        spec = importlib.util.spec_from_loader(_PROBE, loader)
        probe = importlib.util.module_from_spec(spec)
        probe.__package__ = "PQAnalysis.analysis.momentum"
        sys.modules[_PROBE] = probe
        loader.exec_module(probe)

        assert probe.legacy_momentum_file is None
        assert probe.legacy_momentum_norm.__module__.endswith(
            "_momentum_kernel_py"
        )
    finally:
        sys.meta_path.remove(blocker)
        sys.modules.pop(_PROBE, None)
        if had_kernel:
            importlib.import_module(_KERNEL)

    assert sys.modules[Momentum.__module__] is real_module
    assert real_module.legacy_momentum_norm is real_kernel
