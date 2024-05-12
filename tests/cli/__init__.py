import argparse

from PQAnalysis._version import __version__



class ArgparseNamespace(argparse.Namespace):

    def __init__(self, **kwargs):
        self.__dict__.update(progress=False, version=__version__)
        self.__dict__.update(kwargs)
