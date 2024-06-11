"""
A module containing different decorators 
associated to methods of the AtomicSystem class.

...

Functions
---------
check_atoms_pos
    Decorator which checks that the number 
    of atoms is equal to the number of positions.
check_atoms_has_mass
    Decorator which checks that all atoms 
    have mass information.
check_atom_number_setters
    Decorator which checks that the number 
    of atoms is equal to the number of positions.
"""

import numpy as np

from beartype.typing import Any
from decorator import decorator

from .exceptions import (
    AtomicSystemPositionsError, AtomicSystemMassError, AtomicSystemError
)



@decorator
def check_atom_number_setters(func, self, arg_to_set: Any) -> None:
    """
    Decorator which checks that the number of
    atoms is equal to the number of positions.

    Parameters
    ----------
    func : function
        Function to be decorated.
    arg_to_set : Any
        The argument to be set.

    Returns
    -------
    Any
        The result of the function.

    Raises
    ------
    AtomicSystemError
        If the number of atoms is not equal the number of positions.
    """
    if self.n_atoms != np.shape(arg_to_set)[0]:
        self.logger.error(
            (
                "The number of atoms in the AtomicSystem object have "
                "to be equal to the number of atoms in the new array "
                "in order to set the property."
            ),
            exception=AtomicSystemError
        )

    return func(self, arg_to_set)



@decorator
def check_atoms_pos(func, *args, **kwargs):
    """
    Decorator which checks that the number
    of atoms is equal to the number of positions.

    Parameters
    ----------
    func : function
        Function to be decorated.
    *args : Any
        The arguments of the function.
    **kwargs : Any
        The keyword arguments of the function.

    Raises
    ------
    AtomicSystemPositionsError
        If the number of atoms is not equal the number of positions.
    """
    self = args[0]

    if self.pos.shape[0] != len(self.atoms):
        self.logger.error(
            AtomicSystemPositionsError.message,
            exception=AtomicSystemPositionsError
        )

    return func(*args, **kwargs)



@decorator
def check_atoms_has_mass(func, *args, **kwargs):
    """
    Decorator which checks that all atoms have mass information.

    Parameters
    ----------
    func : function
        Function to be decorated.
    *args : Any
        The arguments of the function.
    **kwargs : Any
        The keyword arguments of the function.    

    Raises
    ------
    AtomicSystemMassError
        If any atom does not have mass information.
    """

    self = args[0]

    if not all(atom.mass is not None for atom in self.atoms):
        self.logger.error(
            AtomicSystemMassError.message, exception=AtomicSystemMassError
        )

    return func(*args, **kwargs)
