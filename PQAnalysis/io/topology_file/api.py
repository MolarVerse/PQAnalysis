"""
A module containing the API for reading and writing topology files.
"""

from PQAnalysis.topology import Topology, BondedTopology
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.type_checking import runtime_type_checking

from .topology_file_writer import TopologyFileWriter
from .topology_file_reader import TopologyFileReader



@runtime_type_checking
def write_topology_file(
    bonded_topology: Topology | BondedTopology,
    filename: str | None = None,
    mode: FileWritingMode | str = "w"
) -> None:
    """
    Wrapper function to write a topology file.
    """

    TopologyFileWriter(filename, mode=mode).write(bonded_topology)



@runtime_type_checking
def read_topology_file(filename: str) -> BondedTopology:
    """
    Function to read a topology file.

    Parameters
    ----------
    filename : str
        The filename of the topology file.

    Returns
    -------
    BondedTopology
        The bonded topology.
    """

    return TopologyFileReader(filename).read()
