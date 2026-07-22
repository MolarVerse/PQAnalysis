"""Build the Cython extensions used by PQAnalysis."""

from Cython.Build import cythonize
import numpy as np
from setuptools import Extension, setup

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
            'PQAnalysis.analysis.vacf._vacf_kernel',
            sources=['PQAnalysis/analysis/vacf/_vacf_kernel.pyx'],
            include_dirs=[np.get_include()]
        ),
        Extension(
            'PQAnalysis.analysis.rdf._rdf_kernel',
            sources=['PQAnalysis/analysis/rdf/_rdf_kernel.pyx'],
            include_dirs=[np.get_include()],
            # no FMA contraction: the kernel must round exactly like
            # the separate numpy operations it replicates
            extra_compile_args=['-ffp-contract=off']
        ),
        Extension(
            'PQAnalysis.analysis.adf._adf_kernel',
            sources=['PQAnalysis/analysis/adf/_adf_kernel.pyx'],
            include_dirs=[np.get_include()],
            # no FMA contraction: the kernel must round exactly like
            # the separate numpy operations it replicates
            extra_compile_args=['-ffp-contract=off']
        ),
    ],
    language_level=3,
)

setup(ext_modules=extensions)
