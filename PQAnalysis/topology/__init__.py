"""
A package containing classes and functions to handle molecular topologies.

The topology package contains the following submodules:
    
        - residue
        - selection
        - topology
        
The topology package contains the following classes:
        
        - Residue
        - QMResidue
        - Selection
        - Topology
                
The topology package contains the following type hints:

        - SelectionCompatible
        - Residues
                
The topology package contains the following exceptions:
    
        - ResidueError
"""

from .exceptions import ResidueError, TopologyError

from .selection import Selection, SelectionCompatible
from .topology import Topology
from .residue import Residue, Residues, QMResidue
