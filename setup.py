"""
This file is used to build the package. It is used to compile the Cython code
"""

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

module = cythonize(
    Extension(
        'PQAnalysis.io.traj_file.process_lines',
        sources=['PQAnalysis/io/traj_file/process_lines.pyx'],
        include_dirs=[np.get_include()]
    ),
    compiler_directives={'language_level': 3}
)

setup(
    name="PQAnalysis",
    packages=["PQAnalysis"],
    ext_modules=module,
)
