from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

module = cythonize(
    Extension(
        'PQAnalysis.io.traj_file.process_lines',
        sources=['PQAnalysis/io/traj_file/process_lines.pyx'],
        include_dirs=[np.get_include()]
    )
)

setup(
    name="PQAnalysis",
    packages=["PQAnalysis"],
    ext_modules=module,
)
