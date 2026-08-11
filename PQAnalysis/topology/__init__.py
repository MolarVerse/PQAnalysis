"""
A package containing classes and functions to handle molecular/atomic topologies.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:
    from .bonded_topology.angle import Angle
    from .bonded_topology.bond import Bond
    from .bonded_topology.bonded_topology import BondedTopology
    from .bonded_topology.dihedral import Dihedral
    from .exceptions import TopologyError
    from .selection import Selection, SelectionCompatible
    from .topology import Topology

_EXPORTS = {
    "TopologyError": ".exceptions",
    "Bond": ".bonded_topology.bond",
    "Angle": ".bonded_topology.angle",
    "Dihedral": ".bonded_topology.dihedral",
    "BondedTopology": ".bonded_topology.bonded_topology",
    "Selection": ".selection",
    "SelectionCompatible": ".selection",
    "Topology": ".topology",
}

__all__ = list(_EXPORTS)



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
