"""
A package containing classes and functions to handle molecular/atomic topologies.
"""

from .exceptions import TopologyError

from .bonded_topology.bond import Bond
from .bonded_topology.angle import Angle
from .bonded_topology.dihedral import Dihedral
from .bonded_topology.bonded_topology import BondedTopology
from .selection import Selection, SelectionCompatible
from .topology import Topology

# TODO: partially circular --- from .api import generate_shake_topology_file
# TODO: partially circular --- from .shakeTopology import ShakeTopologyGenerator
