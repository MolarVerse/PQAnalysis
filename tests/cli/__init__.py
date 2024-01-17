from PQAnalysis._version import __version__


class ArgparseNamespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def add(self, **kwargs):
        self.__dict__.update(kwargs)
        return self


argparse_namespace = ArgparseNamespace(progress=False, version=__version__)
