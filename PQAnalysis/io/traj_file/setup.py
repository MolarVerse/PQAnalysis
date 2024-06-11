from setuptools import setup, Extension

module = Extension('mytest', sources=['mytest.c'])

setup(
    name='MyTest',
    version='1.0',
    description='Python interface for the FrameReader C library',
    ext_modules=[module],
)
