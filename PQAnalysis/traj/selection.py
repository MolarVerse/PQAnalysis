"""
A module containing the selection class.

...

Classes
-------
Selection
    A class for atom selections.
"""

import numpy as np

from collections.abc import Iterable
from typing import List, Union


class Selection:
    """
    A class for atom selections.
    """

    def __init__(self, selection_param: Union[List[int], List[str], int, str], frame=None):
        """
        Initializes the AtomSelection with the given selection.

        Parameters
        ----------
        selection : np.array or list of int or list of str
            The selection. Can be either a list of int or str, or a np.array.
        frame : Frame, optional
            The frame to select from. Required if selection is a list of str.

        Raises
        ------
        TypeError
            If the selection is not a list of int or str or a single int or str.
        ValueError
            If the selection is a list of str and no frame is provided.
        """

        if isinstance(selection_param, Iterable) and not isinstance(selection_param, str):
            self.selection = np.array(selection_param)
        else:
            self.selection = np.array([selection_param])

        if isinstance(self.selection[0], np.int64):
            self.selection = self.selection
        elif isinstance(self.selection[0], str):
            if frame is None:
                raise ValueError(
                    'Frame must be provided when selection is a string.')
            else:
                self.selection = np.nonzero(
                    np.in1d(frame.atoms, self.selection))[0]
        else:
            raise TypeError('Invalid selection type.')
