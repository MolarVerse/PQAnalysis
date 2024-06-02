"""
A subpackage to handle trajectory files.
"""

from .exceptions import TopologyFileError

from .topology_file_reader import TopologyFileReader
from .topology_file_writer import TopologyFileWriter
