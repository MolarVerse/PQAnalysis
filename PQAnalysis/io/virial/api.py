"""
A module containing the API for reading virial and stress files.
"""

from beartype.typing import List

from .virial_reader import VirialFileReader, StressFileReader


def read_virial_file(filename: str) -> List:
    """
    Read a virial file.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Returns
    -------
    List[np.ndarray]
        The virial data.
    """
    reader = VirialFileReader(filename)
    return reader.read()


def read_stress_file(filename: str) -> List:
    """
    Read a stress file.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Returns
    -------
    List[np.ndarray]
        The stress data.
    """
    reader = StressFileReader(filename)
    return reader.read()