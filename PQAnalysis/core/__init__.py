"""
This package contains the core classes and functions of the package.

The core package contains the following submodules and packages:
    
        - atom
        - cell
        - atomicSystem
        - common
        - exceptions
        
The core package contains the following classes:
        
            - Atom
            - Element
            - Cell
            - AtomicSystem
            
The core package contains the following functions:
        
            - distance

The core package contains the following exceptions:

            - ElementNotFoundError
            - AtomicSystemPositionsError
            - AtomicSystemMassError
            
The core package contains the following type hints:
        
                  - Atoms
                  - Elements
"""

from .common import distance
from .exceptions import ElementNotFoundError, AtomicSystemPositionsError, AtomicSystemMassError

from .atomicSystem import AtomicSystem
from .cell.cell import Cell
from .atom import Atom, Atoms, Element, Elements
