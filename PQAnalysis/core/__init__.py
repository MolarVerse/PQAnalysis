"""
This module contains the core classes and functions of the package.

The core module contains the following submodules:
    
        - atom
        - cell
        - atomicSystem
        - common
        - exceptions
        
The core module contains the following classes:
        
            - Atom
            - Atoms
            - Element
            - Elements
            - Cell
            - AtomicSystem
            
The core module contains the following functions:
        
            - distance

The core module contains the following exceptions:

            - ElementNotFoundError
            - AtomicSystemPositionsError
            - AtomicSystemMassError
"""

from .common import distance
from .exceptions import ElementNotFoundError, AtomicSystemPositionsError, AtomicSystemMassError
from .atom import Atom, Atoms, Element, Elements
from .cell import Cell
from .atomicSystem import AtomicSystem
