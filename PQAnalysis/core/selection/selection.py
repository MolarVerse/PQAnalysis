import numpy as np

from lark import Transformer, Visitor, Tree, Lark
from beartype.typing import List, Tuple

from ... import __base_path__
from ...types import Np1DIntArray
from ..atom import Atom
from ..atomicSystem import AtomicSystem

def selection(string: str, atomic_system: AtomicSystem) -> Np1DIntArray:
    grammar_file = "selection.lark"
    grammar_path = __base_path__ / "grammar"

    parser = Lark.open(grammar_path / grammar_file,
                           propagate_positions=True)
    
    tree = parser.parse(string)
    visitor = SelectionVisitor(atomic_system)
    
    return visitor.visit(tree)
    

class SelectionVisitor(Visitor):
    def __init__(self, atomic_system: AtomicSystem):
        self.atomic_system = atomic_system
        self.selection = []
        self.super().__init__()
        
    def atomtype(self, items) -> Np1DIntArray:
        for index in self.atomic_system._indices_by_atom_type_names([items[0]]):
            self.selection.append(index)
            
        return np.array(self.selection)
    
    def atom(self, items) -> Np1DIntArray:
        for index in self.atomic_system._indices_by_atom([items[0]]):
            self.selection.append(index) 
            
        return np.array(self.selection)
        
    def index(self, items) -> Np1DIntArray:
        self.selection.append(items[0])
        
        return np.array(self.selection)
        
    def indices(self, items) -> Np1DIntArray:
        if len(items) == 2
            for index in range(items[0], items[1]):
                self.selection.append(index)
        elif len(items) == 3:
            for index in range(items[0], items[2], items[1]):
                self.selection.append(index)
                
        return np.array(self.selection)

    def visit(self, tree: Tree) -> Np1DIntArray:
        super().visit(tree)

        return np.array(self.selection)
