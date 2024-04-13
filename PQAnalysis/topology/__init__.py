"""
A package containing classes and functions to handle molecular/atomic topologies.
"""

from .exceptions import TopologyError

from .bond import Bond
from .angle import Angle
from .dihedral import Dihedral
from .bondedTopology import BondedTopology
from .selection import Selection, SelectionCompatible
from .topology import Topology

# TODO: partially circular --- from .api import generate_shake_topology_file
# TODO: partially circular --- from .shakeTopology import ShakeTopologyGenerator
