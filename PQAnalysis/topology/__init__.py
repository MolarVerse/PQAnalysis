"""
A package containing classes and functions to handle molecular/atomic topologies.
"""

from .exceptions import TopologyError

from .selection import Selection, SelectionCompatible
from .topology import Topology
# TODO: partially circular --- from .shakeTopology import ShakeTopologyGenerator
