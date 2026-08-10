"""Bitwise parity tests against the legacy equipartition.jl tool."""

import importlib

import numpy as np
import pytest

from PQAnalysis.analysis.momentum import Momentum
from PQAnalysis.analysis.momentum._momentum_kernel_py import (
    legacy_momentum_norm as legacy_momentum_norm_py,
)
from PQAnalysis.io import TrajectoryReader

from .. import pytestmark  # pylint: disable=unused-import
from .test_momentum import GAS3_LEGACY_NORMS

momentum_module = importlib.import_module(
    "PQAnalysis.analysis.momentum.momentum"
)

try:
    from PQAnalysis.analysis.momentum._momentum_kernel import (  # pylint: disable=import-error
        legacy_momentum_norm as legacy_momentum_norm_compiled,
    )
except ModuleNotFoundError:
    legacy_momentum_norm_compiled = None

KERNELS = [legacy_momentum_norm_py]

if legacy_momentum_norm_compiled is not None:
    KERNELS.append(legacy_momentum_norm_compiled)



@pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
@pytest.mark.parametrize("kernel", KERNELS)
def test_full_precision_equipartition_parity(
    test_with_data_dir,
    monkeypatch,
    kernel,
):
    """Every legacy gas-frame result matches at the float64 bit level."""
    monkeypatch.setattr(momentum_module, "legacy_momentum_norm", kernel)
    monkeypatch.setattr(Momentum, "_batch_max_bytes", 0)

    analysis = Momentum(TrajectoryReader("gas3.vel", md_format="qmcfc"))
    actual = analysis.run()

    assert analysis._raw_reader.dtype == "float64"  # pylint: disable=protected-access
    assert np.array_equal(actual, GAS3_LEGACY_NORMS)



@pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
def test_fused_file_path_matches_streaming(
    test_with_data_dir,
    monkeypatch,
):
    """Fused parsing preserves every exact streaming result bit."""
    if momentum_module.legacy_momentum_file is None:
        pytest.skip("compiled fused parser not available")

    fused = Momentum(TrajectoryReader("gas3.vel", md_format="qmcfc")).run()

    monkeypatch.setattr(Momentum, "_batch_max_bytes", 0)
    streamed = Momentum(TrajectoryReader("gas3.vel", md_format="qmcfc")).run()

    assert np.array_equal(fused, streamed)
