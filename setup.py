"""Build the Cython extension used by PQAnalysis."""

from Cython.Build import cythonize
import numpy as np
from setuptools import Extension, setup

extensions = cythonize(
    Extension(
        'PQAnalysis.io.traj_file.process_lines',
        sources=['PQAnalysis/io/traj_file/process_lines.pyx'],
        include_dirs=[np.get_include()]
    ),
    language_level=3,
)

setup(ext_modules=extensions)
