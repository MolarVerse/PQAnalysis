"""Build the Cython extensions used by PQAnalysis."""

import sys

from Cython.Build import cythonize
import numpy as np
from setuptools import Extension, setup

# MSVC does not contract floating-point operations under its default /fp:precise
NO_FP_CONTRACT = [] if sys.platform == 'win32' else ['-ffp-contract=off']

extensions = cythonize(
    [
        Extension(
            'PQAnalysis.io.traj_file.process_lines',
            sources=['PQAnalysis/io/traj_file/process_lines.pyx'],
            include_dirs=[np.get_include()]
        ),
        Extension(
            'PQAnalysis.io.traj_file._slab_parser',
            sources=['PQAnalysis/io/traj_file/_slab_parser.pyx'],
            include_dirs=[np.get_include()]
        ),
        Extension(
            'PQAnalysis.analysis.msd._msd_kernel',
            sources=['PQAnalysis/analysis/msd/_msd_kernel.pyx'],
            include_dirs=[np.get_include()]
        ),
        Extension(
            'PQAnalysis.analysis.momentum._momentum_kernel',
            sources=['PQAnalysis/analysis/momentum/_momentum_kernel.pyx'],
            include_dirs=[np.get_include()],
            # Julia rounds every multiply and add separately in the
            # equipartition loop.
            extra_compile_args=NO_FP_CONTRACT
        ),
        Extension(
            'PQAnalysis.analysis.vacf._vacf_kernel',
            sources=['PQAnalysis/analysis/vacf/_vacf_kernel.pyx'],
            include_dirs=[np.get_include()],
            # Only the explicit fma calls may fuse operations. This
            # keeps the optimized FreqCalc arithmetic reproducible.
            extra_compile_args=NO_FP_CONTRACT
        ),
        Extension(
            'PQAnalysis.analysis.rdf._rdf_kernel',
            sources=['PQAnalysis/analysis/rdf/_rdf_kernel.pyx'],
            include_dirs=[np.get_include()],
            # no FMA contraction: the kernel must round exactly like
            # the separate numpy operations it replicates
            extra_compile_args=NO_FP_CONTRACT
        ),
    ],
    language_level=3,
)

setup(ext_modules=extensions)
