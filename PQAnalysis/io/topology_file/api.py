from .topologyFileWriter import TopologyFileWriter
from PQAnalysis.topology import Topology, BondedTopology
from PQAnalysis.io import FileWritingMode


def write_topology_file(bonded_topology: Topology | BondedTopology | None,
                        filename: str | None = None,
                        mode: FileWritingMode | str = "w"
                        ) -> None:
    """
    Wrapper function to write a topology file.
    """

    TopologyFileWriter(filename, mode=mode).write(bonded_topology)


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

    return TopologyFileWriter(filename).read()
