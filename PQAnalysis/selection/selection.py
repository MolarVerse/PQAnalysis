from PQAnalysis.utils.instance import isinstance_of_list
from collections.abc import Iterable
import numpy as np


class Selection:
    def __init__(self, selection):
        self.selection = selection
        self.n_atoms = len(selection)


class AtomSelection(Selection):
    def __init__(self, selection, frame=None):

        if isinstance(selection, Iterable):
            selection = np.array(selection)
        else:
            selection = np.array([selection])

        if isinstance(selection[0], int):
            super().__init__(selection)
        elif isinstance(selection[0], str):
            if frame is None:
                raise ValueError(
                    'Frame must be provided when selection is a string.')
            else:
                selection = np.nonzero(np.in1d(frame.atoms, selection))
                super().__init__(selection)
        else:
            raise ValueError('Invalid selection type.')
