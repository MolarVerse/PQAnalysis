"""
Fixtures of the end-to-end tests of the analysis commands.
"""

import importlib

import pytest

#: The compiled kernels and their pure Python fallbacks, wired like the
#: ``except ModuleNotFoundError`` blocks of the importing modules do.
KERNELS = [
    {
        "module": "PQAnalysis.analysis.rdf.rdf",
        "compiled": "PQAnalysis.analysis.rdf._rdf_kernel",
        "fallback": "PQAnalysis.analysis.rdf._rdf_kernel_py",
        "names": [
            "legacy_rdf_batch_histogram",
            "legacy_rdf_frame_histogram",
            "rdf_frame_histogram",
        ],
        "none": [],
    },
    {
        "module": "PQAnalysis.analysis.msd.msd",
        "compiled": "PQAnalysis.analysis.msd._msd_kernel",
        "fallback": "PQAnalysis.analysis.msd._msd_kernel_py",
        "names": ["msd_frame_update"],
        "none": [
            "direct_msd_boundary_states",
            "direct_msd_image_steps",
            "direct_msd_lag_range",
            "legacy_msd_frame_update",
        ],
    },
    {
        "module": "PQAnalysis.analysis.vacf.vacf",
        "compiled": "PQAnalysis.analysis.vacf._vacf_kernel",
        "fallback": "PQAnalysis.analysis.vacf._vacf_kernel_py",
        "names": ["accumulate_frame", "weight_frame"],
        "none": ["direct_correlate_lag_range", "direct_origin_norms"],
    },
    {
        "module": "PQAnalysis.analysis.momentum.momentum",
        "compiled": "PQAnalysis.analysis.momentum._momentum_kernel",
        "fallback": "PQAnalysis.analysis.momentum._momentum_kernel_py",
        "names": ["legacy_momentum_norm"],
        "none": ["legacy_momentum_file"],
    },
    {
        "module": "PQAnalysis.io.traj_file.raw_frame_reader",
        "compiled": "PQAnalysis.io.traj_file._slab_parser",
        "fallback": "PQAnalysis.io.traj_file._slab_parser_py",
        "names": ["scan_header", "parse_body", "parse_xyz_frames"],
        "none": [],
    },
]


def _compiled_available():
    try:
        for kernel in KERNELS:
            importlib.import_module(kernel["compiled"])
    except ModuleNotFoundError:
        return False

    return True


@pytest.fixture(
    params=[
        pytest.param(
            "compiled",
            marks=pytest.mark.skipif(
                not _compiled_available(),
                reason="compiled extensions not available",
            ),
        ),
        "fallback",
    ]
)
def kernel_impl(request, monkeypatch):
    """
    Runs a test once with the compiled kernels and once with the pure Python
    fallbacks wired in.
    """

    if request.param == "fallback":
        for kernel in KERNELS:
            module = importlib.import_module(kernel["module"])
            fallback = importlib.import_module(kernel["fallback"])

            for name in kernel["names"]:
                monkeypatch.setattr(module, name, getattr(fallback, name))
            for name in kernel["none"]:
                monkeypatch.setattr(module, name, None)

    return request.param
