"""
This is the core package of PQAnalysis.

It consists of all classes and functions that are used to handle atomic systems, atoms, elements and cells. Additionally it contains some functions that are used to compute properties of atomic systems within the api module, but can be used directly from the core package.
"""

from .api import distance
from .exceptions import ElementNotFoundError, AtomicSystemPositionsError, AtomicSystemMassError

from .cell.cell import Cell, Cells
from .atom import Atom, Atoms, Element, Elements
from .atomicSystem import AtomicSystem
