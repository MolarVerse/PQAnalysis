from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

module = cythonize(
    Extension(
        'PQAnalysis.io.traj_file.mytest',
        sources=['PQAnalysis/io/traj_file/mytest.pyx'],
        include_dirs=[np.get_include()]
    )
)

setup(
    name="PQAnalysis",
    packages=["PQAnalysis"],
    ext_modules=module,
)
