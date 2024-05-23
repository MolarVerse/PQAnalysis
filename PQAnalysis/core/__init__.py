"""
This is the core package of PQAnalysis.

It consists of all classes and functions that are used to handle 
atomic systems, atoms, elements and cells. Additionally it contains
some functions that are used to compute properties of atomic systems
within the api module, but can be used directly from the core package.
"""

from .exceptions import *

from .cell import Cell, Cells
from .atom import Atom, Atoms, Element, Elements, CustomElement
from .residue import Residue, Residues, QMResidue

from .api import distance
