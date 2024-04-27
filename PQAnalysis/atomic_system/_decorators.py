"""
A module containing different decorators associated to methods of the AtomicSystem class.

...

Functions
---------
check_atoms_pos
    Decorator which checks that the number of atoms is equal to the number of positions.
check_atoms_has_mass
    Decorator which checks that all atoms have mass information.
check_atom_number_setters
    Decorator which checks that the number of atoms is equal to the number of positions.
"""

import numpy as np

from beartype.typing import Any
from decorator import decorator

from .exceptions import AtomicSystemPositionsError, AtomicSystemMassError, AtomicSystemError


@decorator
def check_atom_number_setters(func, self, arg_to_set: Any) -> Any:
    """
    Decorator which checks that the number of atoms is equal to the number of positions.

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
    ValueError
        If the number of atoms is not equal the number of positions.
    """
    if self.n_atoms != np.shape(arg_to_set)[0]:
        raise AtomicSystemError(
            "The number of atoms in the AtomicSystem object have "
            "to be equal to the number of atoms in the new array "
            "in order to set the property."
        )

    return func(self, arg_to_set)


@decorator
def check_atoms_pos(func, *args, **kwargs):
    """
    Decorator which checks that the number of atoms is equal to the number of positions.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Raises
    ------
    AtomicSystemPositionsError
        If the number of atoms is not equal the number of positions.
    """
    self = args[0]

    if self.pos.shape[0] != len(self.atoms):
        raise AtomicSystemPositionsError()

    return func(*args, **kwargs)


@decorator
def check_atoms_has_mass(func, *args, **kwargs):
    """
    Decorator which checks that all atoms have mass information.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Raises
    ------
    ValueError
        If any atom does not have mass information.
    """

    self = args[0]

    if not all([atom.mass is not None for atom in self.atoms]):
        raise AtomicSystemMassError()

    return func(*args, **kwargs)
