"""
A module containing the API for reading virial and stress files.
"""

from beartype.typing import List

from PQAnalysis.type_checking import runtime_type_checking

from .virial_reader import VirialFileReader, StressFileReader



@runtime_type_checking
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



@runtime_type_checking
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
